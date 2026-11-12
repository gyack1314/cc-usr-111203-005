"""分位增量订正（加法口径）与同期秩恢复的纯函数实现。

设计口径（与 docs/data-contract.md 一致）：

* 分组按 (station, season) 独立；OH/MH 只取历史期样本、MF 只取未来期样本，
  校准路径上不存在未来观测；
* 每个未来值按包含干日的未来全体样本取 midrank 概率 p，参考分位用公开的
  type-7 插值规则，端点外不外推；
* 变化信号固定为加法增量 ``Q_OH(p) + MF - Q_MH(p)``，不提供比例缩放口径，
  零雨量不参与除法，因此不会产生无限大；
* 干日频率是显式约束：置零与湿日抬升分别记账，频率校准与总雨量守恒不是
  同一目标；
* 模板恢复只做每站边际值的纯置换，跨站依赖仅跟随事先锁定的模板秩结构。

本模块不读写任何输入文件，输入文件保持不变。
"""

import math
from collections import defaultdict

from .io import date, number

HISTORIC_PERIODS = ("OH", "MH")
ALL_PERIODS = ("OH", "MH", "MF", "template")

# 静态保证声明：同时进入机器可读输出，避免口径只停留在注释里。
CORRECTION_GUARANTEES = {
    "preserves": [
        "逐站变化信号以加法增量 Q_OH(p)+MF-Q_MH(p) 保留，不替换成比例缩放，零雨量不会产生无限大",
        "每个站季节只用本列所属时期样本：OH/MH 用历史期，MF 用未来期，未来观测不进入校准",
        "midrank 概率包含干日，参考分位使用公开的 type-7 插值规则且不在端点外外推",
        "干日频率按 clip(f_OH+f_MF-f_MH,0,1) 显式约束，置零个数与目标干日数一致",
    ],
    "does_not_guarantee": [
        "不保证总雨量守恒：置零损失与湿日抬升分别记账，净水量变化可能非零",
        "不保证跨站同期相关：订正在每个站季节内部独立进行，原始未来日期的跨站对齐未被建模",
        "不保证模板之外的任何同期事件结构，未使用模板时输出只反映逐站边际分布",
    ],
}

RESTORATION_GUARANTEES = {
    "preserves": [
        "恢复前后每站边际值多重集完全相同：恢复只是同一组值的置换",
        "恢复前后每站零日数完全相同，每站窗内总量完全相同，全窗汇水区总量完全相同",
        "输出时间身份仍是未来窗日期，锁定模板日期只提供秩位置，不覆盖未来日期",
    ],
    "does_not_guarantee": [
        "不承诺精确跨站相关系数：订正值并列时按 (值, 未来日期) 稳定决胜，并列如何落位影响具体相关",
        "不重建真实同期降雨事件，跨站依赖只保证跟随该锁定模板的秩结构",
        "不验证模板本身的代表性；模板必须事先锁定，流程不允许按站挑年份或按未来观测挑相似段",
    ],
}


class TemplateRestorationError(ValueError):
    """锁定模板不满足窗长、站日完整性或历史锁定要求，整次恢复拒绝。"""


def _nz(x):
    """归一化负零，保证 JSON 输出不出现 -0.0；非有限值在此之前已被拒绝。"""
    return 0.0 if x == 0.0 else float(x)


def quantile_type7(sorted_values, p):
    """公开的 type-7 分位插值：``h=(n-1)*p`` 线性插值，n=1 取唯一值。

    `sorted_values` 必须已升序排列且非空。端点外不向外外推，越界概率钉在
    端点值上（midrank 概率实际落在 (0,1) 内，不会触发该分支）。
    """
    if not sorted_values:
        raise ValueError("分位计算需要非空样本")
    if not 0.0 <= p <= 1.0:
        raise ValueError("分位概率越界")
    n = len(sorted_values)
    if n == 1:
        return float(sorted_values[0])
    h = (n - 1) * p
    lo = math.floor(h)
    if lo <= 0:
        lo = 0
    if lo >= n - 1:
        return float(sorted_values[-1])
    frac = h - lo
    low = sorted_values[lo]
    high = sorted_values[lo + 1]
    return float(low + (high - low) * frac)


def midrank_probabilities(values):
    """含干日的 midrank：``p=(小于x的个数+0.5*等于x的个数)/n``。"""
    n = len(values)
    if n == 0:
        raise ValueError("midrank 需要非空未来样本")
    ordered = sorted(values)
    probs = []
    for x in values:
        less = equal = 0
        for y in ordered:
            if y < x:
                less += 1
            elif y == x:
                equal += 1
        probs.append((less + 0.5 * equal) / n)
    return probs


def half_up_count(ratio, n):
    """目标干日个数：``floor(ratio*n+0.5)``，四舍五入、恰半向上，并钉在 [0,n]。"""
    if n <= 0:
        raise ValueError("样本数必须为正")
    return max(0, min(n, math.floor(ratio * n + 0.5)))


def _parse_config(v):
    if not isinstance(v, dict):
        raise ValueError("顶层资料应为对象")
    for key in ("OH", "MH", "MF"):
        if not isinstance(v.get(key), list) or not v[key]:
            raise ValueError(f"{key} 需要非空记录列表")
    dry_threshold = number(v["dry_threshold"])
    if dry_threshold < 0:
        raise ValueError("干日阈值不得为负")
    wet_floor = number(v["wet_floor"])
    if wet_floor <= dry_threshold:
        raise ValueError("湿日下限应严格高于干日阈值")
    rule = v.get("quantile_rule")
    if rule != "type7":
        raise ValueError("仅支持公开的 type7 分位规则")
    return dry_threshold, wet_floor, rule


def _parse_rows(v):
    """解析四个时期的行，做站日唯一性与基本数值核对，返回 (规范化行, 分组)。"""
    parsed = {p: [] for p in ALL_PERIODS}
    seen = set()
    for period in ALL_PERIODS:
        rows = v.get(period)
        if period == "template":
            # 订正值本身不依赖模板；模板只在恢复阶段强制要求。
            if rows is None:
                continue
            if not isinstance(rows, list) or not rows:
                raise TemplateRestorationError("模板需要非空记录列表")
        elif not isinstance(rows, list) or not rows:
            raise ValueError(f"{period} 需要非空记录列表")
        for row in rows:
            if not isinstance(row, dict):
                raise ValueError(f"{period} 行应为对象")
            station = row["station"]
            season = row["season"]
            if not isinstance(station, str) or not station:
                raise ValueError("站点身份缺失")
            if not isinstance(season, str) or not season:
                raise ValueError("季节身份缺失")
            day = date(row["day"])
            mm = number(row["mm"])
            if mm < 0:
                raise ValueError("降水不得为负")
            key = (period, station, day)
            if key in seen:
                raise ValueError("同列站日重复")
            seen.add(key)
            parsed[period].append({"station": station, "season": season,
                                   "day": day, "mm": mm})
    groups = defaultdict(lambda: {p: [] for p in ALL_PERIODS})
    for period in ALL_PERIODS:
        for row in parsed[period]:
            groups[(row["station"], row["season"])][period].append(row)
    return parsed, groups


def _validate_correction_grid(groups):
    """订正网格核对：每站季节三列非空，同季节各站未来日期集合一致。"""
    stations = set()
    seasons_by_station = defaultdict(set)
    mf_days_by_season = {}
    for (station, season), buckets in groups.items():
        stations.add(station)
        seasons_by_station[station].add(season)
        for period in HISTORIC_PERIODS + ("MF",):
            if not buckets[period]:
                raise ValueError(f"站点 {station} 季节 {season} 缺少 {period} 样本")
        mf_days = frozenset(r["day"] for r in buckets["MF"])
        if season in mf_days_by_season and mf_days_by_season[season] != mf_days:
            raise ValueError(f"季节 {season} 各站未来日期集合不一致")
        mf_days_by_season[season] = mf_days
    base_seasons = None
    for station in stations:
        seasons = seasons_by_station[station]
        if base_seasons is None:
            base_seasons = seasons
        elif seasons != base_seasons:
            raise ValueError(f"站点 {station} 的季节集合与其他站不一致")
    return sorted(stations), sorted(base_seasons)


def _dry_frequency(values, dry_threshold):
    return sum(1 for x in values if x <= dry_threshold) / len(values)


def _correct_group(station, season, buckets, dry_threshold, wet_floor):
    oh = sorted(r["mm"] for r in buckets["OH"])
    mh = sorted(r["mm"] for r in buckets["MH"])
    mf_rows = sorted(buckets["MF"], key=lambda r: r["day"])
    n = len(mf_rows)
    mf_values = [r["mm"] for r in mf_rows]
    probabilities = midrank_probabilities(mf_values)

    f_oh = _dry_frequency(oh, dry_threshold)
    f_mh = _dry_frequency(mh, dry_threshold)
    f_mf = _dry_frequency(mf_values, dry_threshold)
    target_ratio = min(1.0, max(0.0, f_oh + f_mf - f_mh))
    target_dry_count = half_up_count(target_ratio, n)

    records = []
    for row, p in zip(mf_rows, probabilities):
        q_oh = quantile_type7(oh, p)
        q_mh = quantile_type7(mh, p)
        # 加法增量是选定的变化信号口径；显式不做 MF*Q_OH/Q_MH 一类比例缩放。
        raw = q_oh + row["mm"] - q_mh
        if not math.isfinite(raw):
            raise ValueError("订正值非有限，拒绝输出")
        clipped = max(raw, 0.0)
        records.append({
            "day": row["day"].isoformat(),
            "station": station,
            "season": season,
            "original_mm": float(row["mm"]),
            "midrank_p": p,
            "quantile_oh_mm": q_oh,
            "quantile_mh_mm": q_mh,
            "raw_additive_mm": float(raw),
            "clipped_mm": float(clipped),
            "truncation_loss_mm": float(clipped - raw),
            "zeroed": False,
            "lifted": False,
            "zeroing_change_mm": 0.0,
            "lifting_change_mm": 0.0,
            "final_mm": float(clipped),
        })

    # 置零：按 (原始订正值, 未来日期) 取最小的目标个数；日期决胜保证确定性。
    order = sorted(range(n),
                   key=lambda i: (records[i]["raw_additive_mm"], records[i]["day"]))
    for i in order[:target_dry_count]:
        records[i]["zeroed"] = True
        records[i]["zeroing_change_mm"] = _nz(-records[i]["clipped_mm"])

    # 湿日抬升：未置零且不高于干日阈值的值抬到 wet_floor（wet_floor 已验证严格更大）。
    for rec in records:
        if rec["zeroed"]:
            continue
        if rec["clipped_mm"] <= dry_threshold:
            rec["lifted"] = True
            rec["lifting_change_mm"] = _nz(wet_floor - rec["clipped_mm"])

    for rec in records:
        if rec["zeroed"]:
            rec["final_mm"] = 0.0
        elif rec["lifted"]:
            rec["final_mm"] = _nz(wet_floor)

    ledger = {
        "days": n,
        "original_total_mm": _nz(sum(mf_values)),
        "clipped_total_mm": _nz(sum(r["clipped_mm"] for r in records)),
        "final_total_mm": _nz(sum(r["final_mm"] for r in records)),
        "truncation_loss_mm": _nz(sum(r["truncation_loss_mm"] for r in records)),
        "zeroing_change_mm": _nz(sum(r["zeroing_change_mm"] for r in records)),
        "lifting_change_mm": _nz(sum(r["lifting_change_mm"] for r in records)),
    }
    ledger["net_frequency_adjustment_mm"] = _nz(
        ledger["zeroing_change_mm"] + ledger["lifting_change_mm"])
    return {
        "station": station,
        "season": season,
        "dry_frequency_oh": f_oh,
        "dry_frequency_mh": f_mh,
        "dry_frequency_mf": f_mf,
        "target_dry_ratio": target_ratio,
        "target_dry_count": target_dry_count,
        "final_dry_count": sum(1 for r in records if r["final_mm"] <= dry_threshold),
        "ledger": ledger,
        "records": records,
    }


def correct_document(v):
    """对整份资料执行逐站季节分位增量订正，不读取也不修改模板。"""
    dry_threshold, wet_floor, rule = _parse_config(v)
    _parsed, groups = _parse_rows(v)
    stations, seasons = _validate_correction_grid(groups)
    out_groups = []
    for station in stations:
        for season in seasons:
            out_groups.append(_correct_group(
                station, season, groups[(station, season)],
                dry_threshold, wet_floor))

    station_ledger = []
    for station in stations:
        mine = [g for g in out_groups if g["station"] == station]
        merged = defaultdict(float)
        days = 0
        for g in mine:
            days += g["ledger"]["days"]
            for key, value in g["ledger"].items():
                if key != "days":
                    merged[key] += value
        station_ledger.append({
            "station": station,
            "days": days,
            "dry_count_target": sum(g["target_dry_count"] for g in mine),
            "dry_count_final": sum(g["final_dry_count"] for g in mine),
            "ledger": dict(merged),
        })

    return {
        "method": "additive_quantile_delta",
        "quantile_rule": rule,
        "dry_threshold": dry_threshold,
        "wet_floor": wet_floor,
        "groups": out_groups,
        "station_ledger": station_ledger,
        "guarantees": CORRECTION_GUARANTEES,
    }


def _validate_template(v, groups, corrected_groups):
    """模板锁定规则核对，任何不满足都整次拒绝。"""
    window = v.get("template_window")
    if not isinstance(window, list) or len(window) != 2:
        raise TemplateRestorationError("template_window 应为起止两个日期")
    start, end = date(window[0]), date(window[1])
    if start > end:
        raise TemplateRestorationError("模板窗起止颠倒")

    template_days = {r["day"] for r in groups[next(iter(groups))]["template"]}
    if not template_days:
        raise TemplateRestorationError("模板为空")
    if any(not (start <= d <= end) for d in template_days):
        raise TemplateRestorationError("模板日期超出锁定窗")

    future_days_all = set()
    template_days_by_season = {}
    future_days_by_season = {}
    for g in corrected_groups:
        key = (g["station"], g["season"])
        buckets = groups[key]
        trows = buckets["template"]
        if not trows:
            raise TemplateRestorationError(
                f"站点 {g['station']} 季节 {g['season']} 缺模板站日，整次恢复拒绝")
        if len(trows) != len(buckets["MF"]):
            raise TemplateRestorationError(
                f"站点 {g['station']} 季节 {g['season']} 模板窗长度与未来窗不一致")
        t_days = frozenset(r["day"] for r in trows)
        f_days = frozenset(r["day"] for r in buckets["MF"])
        season = g["season"]
        if season in template_days_by_season and template_days_by_season[season] != t_days:
            raise TemplateRestorationError("模板日期集合必须各站相同，不允许按站挑年份")
        if season in future_days_by_season and future_days_by_season[season] != f_days:
            raise TemplateRestorationError("未来日期集合必须各站相同")
        template_days_by_season[season] = t_days
        future_days_by_season[season] = f_days
        future_days_all |= f_days

    # 历史锁定：模板只能取自未来窗之前，杜绝按未来观测挑相似段。
    if max(template_days) >= min(future_days_all):
        raise TemplateRestorationError("模板日期必须整体早于未来窗")


def apply_restoration(v, corrected):
    """按锁定模板的秩把订正值放回未来日期；不满足锁定规则即整次拒绝。"""
    if not v.get("template"):
        raise TemplateRestorationError("缺少锁定模板，整次恢复拒绝")
    _parsed, groups = _parse_rows(v)
    _validate_template(v, groups, corrected["groups"])

    restored_records = []
    station_rows = []
    daily_before = defaultdict(float)
    daily_after = defaultdict(float)

    for g in corrected["groups"]:
        station, season = g["station"], g["season"]
        trows = sorted(groups[(station, season)]["template"], key=lambda r: r["day"])
        # 模板秩：并列按模板日期排序。
        rank_order = sorted(trows, key=lambda r: (r["mm"], r["day"]))
        rank_of_day = {r["day"]: idx for idx, r in enumerate(rank_order)}

        records = sorted(g["records"], key=lambda r: r["day"])
        future_days = [r["day"] for r in records]
        # 订正有序值：并列按原未来日期稳定决胜。
        value_order = sorted(records, key=lambda r: (r["final_mm"], r["day"]))

        before_by_day = {r["day"]: r["final_mm"] for r in records}
        after_by_day = {}
        for position, trow in enumerate(trows):
            source = value_order[rank_of_day[trow["day"]]]
            after_by_day[future_days[position]] = source["final_mm"]

        for day in future_days:
            before_mm = before_by_day[day]
            after_mm = after_by_day[day]
            restored_records.append({
                "day": day,
                "station": station,
                "season": season,
                "before_restore_mm": float(before_mm),
                "restored_mm": float(after_mm),
            })
            daily_before[day] += before_mm
            daily_after[day] += after_mm

        before_values = sorted(before_by_day.values())
        after_values = sorted(after_by_day.values())
        station_rows.append({
            "station": station,
            "days": len(future_days),
            "multiset_unchanged": before_values == after_values,
            "dry_count_before": sum(1 for x in before_values if x <= corrected["dry_threshold"]),
            "dry_count_after": sum(1 for x in after_values if x <= corrected["dry_threshold"]),
            "total_before_restore_mm": _nz(sum(before_values)),
            "total_after_restore_mm": _nz(sum(after_values)),
            "total_change_mm": _nz(sum(after_values) - sum(before_values)),
        })

    daily = []
    for day in sorted(daily_before):
        daily.append({
            "day": day,
            "total_before_restore_mm": _nz(daily_before[day]),
            "total_after_restore_mm": _nz(daily_after[day]),
            "difference_mm": _nz(daily_after[day] - daily_before[day]),
        })

    return {
        "accepted": True,
        "template_window": v["template_window"],
        "station_ledger": station_rows,
        "watershed_daily_totals": daily,
        "window_total_before_restore_mm": _nz(sum(daily_before.values())),
        "window_total_after_restore_mm": _nz(sum(daily_after.values())),
        "records": sorted(restored_records, key=lambda r: (r["day"], r["station"])),
        "guarantees": RESTORATION_GUARANTEES,
    }
