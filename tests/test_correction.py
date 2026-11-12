"""分位增量订正与模板秩恢复的口径、不变量与拒绝路径测试。"""
import io
import json
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from climate_data.io import load
from climate_data.correction import (
    CORRECTION_GUARANTEES,
    RESTORATION_GUARANTEES,
    TemplateRestorationError,
    apply_restoration,
    correct_document,
    half_up_count,
    midrank_probabilities,
    quantile_type7,
)
from climate_data.cli import main

DATA = Path("data")


def group(corrected, station, season):
    hits = [g for g in corrected["groups"]
            if g["station"] == station and g["season"] == season]
    assert len(hits) == 1
    return hits[0]


def records_by_day(g):
    return {r["day"]: r for r in g["records"]}


class QuantileAndRankTests(unittest.TestCase):
    def test_type7_matches_h_definition(self):
        vals = [0.0, 1.0, 3.0, 4.0]
        # p=0.75 -> h=2.25 -> 3 + 0.25*(4-3)
        self.assertAlmostEqual(quantile_type7(vals, 0.75), 3.25)
        self.assertEqual(quantile_type7([42.0], 0.5), 42.0)
        # 端点不外推
        self.assertEqual(quantile_type7(vals, 0.0), 0.0)
        self.assertEqual(quantile_type7(vals, 1.0), 4.0)
        with self.assertRaises(ValueError):
            quantile_type7([], 0.5)

    def test_midrank_includes_ties_and_dry_days(self):
        # (小于x个数 + 0.5*等于x个数)/n
        self.assertEqual(midrank_probabilities([0, 0, 1]), [1 / 3, 1 / 3, 5 / 6])
        self.assertEqual(midrank_probabilities([5] * 4), [0.5] * 4)

    def test_half_up_half_rounds_up(self):
        self.assertEqual(half_up_count(0.5, 7), 4)   # 3.5 -> 4
        self.assertEqual(half_up_count(0.125, 8), 1)  # 1.0
        self.assertEqual(half_up_count(0.0, 8), 0)
        self.assertEqual(half_up_count(1.0, 8), 8)


class NormalScenarioTests(unittest.TestCase):
    def setUp(self):
        self.value, self.digest = load(DATA / "normal.json")
        self.corrected = correct_document(self.value)

    def test_method_and_output_fields(self):
        self.assertEqual(self.corrected["method"], "additive_quantile_delta")
        self.assertEqual(self.corrected["quantile_rule"], "type7")
        g = group(self.corrected, "coast", "summer")
        r = records_by_day(g)["2030-06-02"]
        for key in ("original_mm", "midrank_p", "quantile_oh_mm",
                    "quantile_mh_mm", "raw_additive_mm", "clipped_mm",
                    "truncation_loss_mm", "final_mm"):
            self.assertIn(key, r)

    def test_hand_computed_coast_values(self):
        g = group(self.corrected, "coast", "summer")
        by_day = records_by_day(g)
        r02 = by_day["2030-06-02"]
        self.assertEqual(r02["midrank_p"], 0.5625)
        self.assertAlmostEqual(r02["quantile_oh_mm"], 1.9375)
        self.assertAlmostEqual(r02["quantile_mh_mm"], 0.9375)
        self.assertEqual(r02["raw_additive_mm"], 4.0)
        r01 = by_day["2030-06-01"]
        self.assertEqual(r01["midrank_p"], 0.125)
        self.assertTrue(r01["zeroed"])
        self.assertEqual(r01["final_mm"], 0.0)

    def test_target_dry_frequency_clip_formula(self):
        coast = group(self.corrected, "coast", "summer")
        # f_oh=.375 f_mf=.25 f_mh=.5 -> .125 -> 1 of 8
        self.assertEqual(coast["target_dry_ratio"], 0.125)
        self.assertEqual(coast["target_dry_count"], 1)
        self.assertEqual(coast["final_dry_count"], 1)
        inland = group(self.corrected, "inland", "summer")
        # 0.25+0.125-0.5 = -0.125 -> clip 0
        self.assertEqual(inland["target_dry_ratio"], 0.0)
        self.assertEqual(inland["target_dry_count"], 0)

    def test_wet_floor_lift_bookkeeping(self):
        # 目标干日为 0 时不置零，但低于阈值的湿日仍要抬到 wet_floor：
        # inland 的 2030-06-05 原始订正为 0.0 但保留湿日身份，抬到 0.2。
        inland = group(self.corrected, "inland", "summer")
        self.assertEqual(inland["ledger"]["zeroing_change_mm"], 0.0)
        self.assertEqual(inland["ledger"]["lifting_change_mm"], 0.2)
        g = group(self.corrected, "coast", "summer")
        self.assertEqual(g["ledger"]["lifting_change_mm"], 0.2)
        lifted = [r for r in g["records"] if r["lifted"]]
        self.assertEqual(len(lifted), 1)
        self.assertEqual(lifted[0]["final_mm"], 0.2)
        self.assertEqual(
            len([r for r in inland["records"] if r["lifted"]]), 1)

    def test_frequency_calibration_does_not_claim_mass_conservation(self):
        g = group(self.corrected, "coast", "summer")
        # 抬升增加水量，净频率调整非零即证明不守恒
        self.assertEqual(g["ledger"]["net_frequency_adjustment_mm"], 0.2)
        self.assertNotEqual(g["ledger"]["original_total_mm"],
                            g["ledger"]["final_total_mm"])
        text = "".join(self.corrected["guarantees"]["does_not_guarantee"])
        self.assertIn("不保证总雨量守恒", text)

    def test_station_ledger_reconciles_groups(self):
        coast = next(s for s in self.corrected["station_ledger"]
                     if s["station"] == "coast")
        g = group(self.corrected, "coast", "summer")
        self.assertEqual(coast["ledger"]["final_total_mm"],
                         g["ledger"]["final_total_mm"])
        self.assertEqual(coast["dry_count_final"], g["final_dry_count"])


class RestorationTests(unittest.TestCase):
    def setUp(self):
        self.value, _ = load(DATA / "normal.json")
        self.corrected = correct_document(self.value)
        self.restored = apply_restoration(self.value, self.corrected)

    def test_marginal_multiset_and_zero_count_invariant(self):
        for row in self.restored["station_ledger"]:
            self.assertTrue(row["multiset_unchanged"])
            self.assertEqual(row["dry_count_before"], row["dry_count_after"])
            self.assertEqual(row["total_change_mm"], 0.0)
        self.assertAlmostEqual(self.restored["window_total_before_restore_mm"],
                               self.restored["window_total_after_restore_mm"])

    def test_daily_totals_are_side_by_side_and_reorder_per_station(self):
        # 汇水区日总量只保证窗总量守恒；逐站置换不保留跨站加总的多重集，
        # 这正是“不承诺精确跨站相关系数”的含义。
        before = [d["total_before_restore_mm"]
                  for d in self.restored["watershed_daily_totals"]]
        after = [d["total_after_restore_mm"]
                 for d in self.restored["watershed_daily_totals"]]
        self.assertAlmostEqual(sum(before), sum(after))
        self.assertNotEqual(before, after)
        grand = sum(d["difference_mm"] for d in self.restored["watershed_daily_totals"])
        self.assertAlmostEqual(grand, 0.0)

    def test_output_keeps_future_dates_not_template_dates(self):
        future_days = {r["day"] for r in self.corrected["groups"][0]["records"]}
        out_days = {r["day"] for r in self.restored["records"]
                    if r["station"] == "coast"}
        self.assertEqual(future_days, out_days)
        self.assertNotIn("2016-06-01", out_days)

    def test_tied_template_uses_date_break_tie_deterministically(self):
        value, _ = load(DATA / "tied-template.json")
        corrected = correct_document(value)
        r1 = apply_restoration(value, corrected)
        r2 = apply_restoration(value, corrected)
        self.assertEqual(
            [(x["day"], x["station"], x["restored_mm"]) for x in r1["records"]],
            [(x["day"], x["station"], x["restored_mm"]) for x in r2["records"]],
        )
        self.assertTrue(all(s["multiset_unchanged"] for s in r1["station_ledger"]))
        # 全并列模板：秩完全由模板日期决定，结构可计算但不承诺相关系数
        self.assertIn("不承诺精确跨站相关系数",
                      "".join(r1["guarantees"]["does_not_guarantee"]))

    def test_incomplete_template_is_whole_run_rejected(self):
        value, _ = load(DATA / "incomplete-template.json")
        with self.assertRaises(TemplateRestorationError):
            apply_restoration(value, correct_document(value))

    def test_no_future_observation_in_calibration(self):
        # 复制后扰动 MF 之外再扰动 MF，只有后者应改变结果
        value, _ = load(DATA / "normal.json")
        baseline = json.dumps(correct_document(value), sort_keys=True)
        perturbed = json.loads(json.dumps(value))
        perturbed["OH"][0]["mm"] = 5
        self.assertNotEqual(baseline,
                            json.dumps(correct_document(perturbed), sort_keys=True))
        mf_perturbed = json.loads(json.dumps(value))
        mf_perturbed["MF"][0]["mm"] = 9
        self.assertNotEqual(baseline,
                            json.dumps(correct_document(mf_perturbed), sort_keys=True))


class EdgeScenarioTests(unittest.TestCase):
    def test_all_dry_constant_midrank_and_no_infinity(self):
        value, _ = load(DATA / "all-dry.json")
        corrected = correct_document(value)
        for g in corrected["groups"]:
            self.assertTrue(all(r["midrank_p"] == 0.5 for r in g["records"]))
            self.assertTrue(all(r["final_mm"] == 0.0 for r in g["records"]))
            self.assertEqual(g["target_dry_count"], g["ledger"]["days"])
        restored = apply_restoration(value, corrected)
        self.assertTrue(all(s["multiset_unchanged"]
                            for s in restored["station_ledger"]))

    def test_constant_series_keeps_additive_signal(self):
        value, _ = load(DATA / "constant.json")
        corrected = correct_document(value)
        coast = group(corrected, "coast", "summer")
        # Q_OH(0.5)=2, Q_MH(0.5)=1, MF=4 -> 2+4-1=5
        self.assertTrue(all(r["final_mm"] == 5.0 for r in coast["records"]))
        inland = group(corrected, "inland", "summer")
        # 3+5-2=6
        self.assertTrue(all(r["final_mm"] == 6.0 for r in inland["records"]))
        # 全程无干日：频率约束不动作
        self.assertEqual(coast["target_dry_count"], 0)

    def test_short_season_n1_path_and_multiple_seasons(self):
        value, _ = load(DATA / "short-season.json")
        corrected = correct_document(value)
        winter_coast = group(corrected, "coast", "winter")
        self.assertEqual(winter_coast["ledger"]["days"], 1)
        self.assertEqual(winter_coast["records"][0]["midrank_p"], 0.5)
        self.assertEqual(winter_coast["records"][0]["final_mm"], 0.0)
        winter_inland = group(corrected, "inland", "winter")
        # n=1 唯一值：1 + 2 - 0 = 3
        self.assertEqual(winter_inland["records"][0]["final_mm"], 3.0)
        restored = apply_restoration(value, corrected)
        self.assertTrue(all(s["multiset_unchanged"]
                            for s in restored["station_ledger"]))


    def test_zeroing_tie_breaks_by_original_future_date(self):
        # OH 1 干 / MH 2 干 / MF 2 干：ratio=.25 -> 目标 1 个干日；
        # 两个并列最小原始订正值中，只置零日期更早的一个，确定性可复算。
        value = {
            "dry_threshold": 0.1, "wet_floor": 0.2, "quantile_rule": "type7",
            "OH": [{"day": f"2017-06-0{d}", "station": "s", "season": "summer",
                    "mm": m} for d, m in zip(range(1, 5), [0, 2, 2, 3])],
            "MH": [{"day": f"2018-06-0{d}", "station": "s", "season": "summer",
                    "mm": m} for d, m in zip(range(1, 5), [0, 0, 1, 2])],
            "MF": [{"day": f"2030-06-0{d}", "station": "s", "season": "summer",
                    "mm": m} for d, m in zip(range(1, 5), [0, 0, 4, 5])],
        }
        g = correct_document(value)["groups"][0]
        self.assertEqual(g["target_dry_count"], 1)
        zeroed = [r["day"] for r in g["records"] if r["zeroed"]]
        self.assertEqual(zeroed, ["2030-06-01"])
        tied = next(r for r in g["records"] if r["day"] == "2030-06-02")
        # 并列另一个保持湿日身份，值不变（高于阈值，无需抬升）
        self.assertFalse(tied["zeroed"])
        self.assertFalse(tied["lifted"])
        self.assertEqual(tied["final_mm"], tied["raw_additive_mm"])
        self.assertEqual(g["final_dry_count"], 1)

    def test_correction_does_not_require_template(self):
        value, _ = load(DATA / "normal.json")
        del value["template"]
        corrected = correct_document(value)
        self.assertTrue(corrected["groups"])


class ValidationAndCliTests(unittest.TestCase):
    def test_invalid_domain_rejected(self):
        value, _ = load(DATA / "invalid-domain.json")
        with self.assertRaises((ValueError, KeyError, TypeError)):
            correct_document(value)

    def test_wet_floor_must_exceed_dry_threshold(self):
        value, _ = load(DATA / "invalid-domain.json")
        # 该文件 wet_floor=0.05 < dry_threshold=0.1
        self.assertLess(value["wet_floor"], value["dry_threshold"])

    def test_unsupported_quantile_rule_rejected(self):
        value, _ = load(DATA / "normal.json")
        value["quantile_rule"] = "type1"
        with self.assertRaises(ValueError):
            correct_document(value)

    def test_restoration_without_template_rejected(self):
        value, _ = load(DATA / "normal.json")
        del value["template"]
        with self.assertRaises(TemplateRestorationError):
            apply_restoration(value, correct_document(value))

    def test_cli_inspect_unchanged_and_correct_emits_json(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = main(["inspect", str(DATA / "normal.json")])
        self.assertEqual(rc, 0)
        payload = json.loads(buf.getvalue())
        self.assertIn("source_sha256", payload)
        self.assertIn("summary", payload)

        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = main(["correct", str(DATA / "normal.json")])
        self.assertEqual(rc, 0)
        payload = json.loads(buf.getvalue())
        self.assertEqual(payload["correction"]["method"], "additive_quantile_delta")
        self.assertEqual(len(payload["source_sha256"]), 64)

    def test_cli_restore_refusal_exit_code(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = main(["restore", str(DATA / "incomplete-template.json")])
        self.assertEqual(rc, 3)
        payload = json.loads(buf.getvalue())
        self.assertFalse(payload["restoration"]["accepted"])
        self.assertIn("reason", payload["restoration"])

    def test_source_files_remain_byte_identical(self):
        before = {p.name: p.read_bytes() for p in DATA.glob("*.json")}
        for p in DATA.glob("*.json"):
            if p.name.startswith("invalid-"):
                try:
                    main(["correct", str(p)])
                except (ValueError, KeyError, TypeError):
                    pass
            else:
                buf = io.StringIO()
                with redirect_stdout(buf):
                    main(["restore", str(p)])
        after = {p.name: p.read_bytes() for p in DATA.glob("*.json")}
        self.assertEqual(before, after)

    def test_guarantee_blocks_present(self):
        value, _ = load(DATA / "normal.json")
        corrected = correct_document(value)
        restored = apply_restoration(value, corrected)
        self.assertEqual(corrected["guarantees"], CORRECTION_GUARANTEES)
        self.assertEqual(restored["guarantees"], RESTORATION_GUARANTEES)
        self.assertTrue(corrected["guarantees"]["preserves"])
        self.assertTrue(corrected["guarantees"]["does_not_guarantee"])


if __name__ == "__main__":
    unittest.main()
