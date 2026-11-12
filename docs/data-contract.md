# 降水情景口径

OH、MH、MF、template各行为station/day/season/mm，mm非负。每站季节样本非空，无隐式缺测填补；输入缺少所需站季节时拒绝订正。midrank对某未来值x用p=(小于x的个数+0.5*等于x的个数)/n，干日（mm<=dry_threshold）留在样本内参与排名。参考分位type-7以(n-1)*p插值，n=1为唯一值，不在端点外外推。

## 分位增量订正（climate_data.correction）

* 分组按(station, season)独立；OH与MH只用历史期列，MF只用未来期列，未来观测不进入校准路径。对MF每个值取含干日的未来全体样本midrank概率p，再按同一type-7规则取Q_OH(p)、Q_MH(p)。
* 变化信号固定为加法增量：raw=Q_OH(p)+MF-Q_MH(p)，先截到0（clipped=max(raw,0)，truncation_loss=clipped-raw）。项目不提供比例缩放口径，零雨量不参与除法，因此不会产生无限大。
* 输出逐行原值original_mm、midrank_p、quantile_oh_mm、quantile_mh_mm、raw_additive_mm、clipped_mm、truncation_loss_mm与final_mm，并标记zeroed/lifted。

## 干日频率约束（显式，不与总量守恒混用）

* 目标干日比例为clip(f_OH+f_MF-f_MH,0,1)，频率取各自列实际样本频率；目标干日个数为floor(比例*未来天数+0.5)，即四舍五入、恰半向上。
* 置零选择按(原始订正值raw, MF日期)升序取目标个数，日期是确定性决胜；未置零且clipped<=dry_threshold的湿日抬到wet_floor，wet_floor在输入层验证严格大于dry_threshold。
* 置零水量变化zeroing_change与抬升水量变化lifting_change分别记账，另有截断损失truncation_loss。频率校准只承诺干日个数达标，不承诺总雨量守恒；净水量变化可以非零。

## 锁定模板的同期秩恢复

* 模板必须事先锁定：一个窗（template_window起止）、各站同一模板日期集合、窗长与该季节未来天数相同、模板日期整体早于未来窗。流程不允许为各站分别挑年份，也不允许根据未来观测挑最相似的一段。模板缺任何站日、窗长不符或日期越界，整次恢复拒绝（CLI返回码3），订正版仍随拒绝结果给出。
* 模板恢复在每站季节独立进行；该季节的模板窗、日期集合与所有站相同。模板按(mm,day)排序得到秩，并列按模板日期决胜；未来最终值按(value,MF日期)稳定排序，第i个模板秩位置取回第i个有序值，再按未来窗日期顺序输出。模板日期只提供秩位置，不覆盖未来日期。
* 恢复是纯置换：恢复前后每站边际值多重集、零日数、每站窗内总量完全相同，全窗汇水区总量相同。跨站依赖只保证跟随该模板的秩结构；订正值并列时按(值,未来日期)决胜，不承诺精确跨站相关系数，也不保留汇水区日总量的多重集。模板本身的代表性不在本流程验证范围内。

现有domain不求分位、不订正、不恢复排序，仅核对原始资料。
