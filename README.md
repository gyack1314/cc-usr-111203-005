# 降水情景的分位增量订正与同期秩恢复

这是独立的气候资料Python项目。`data/`保存合成资料，`climate_data.io`负责严格读取、公历日期和来源字节摘要，`climate_data.domain`提供现有描述统计，`climate_data.correction`实现加法口径的分位增量订正、干日频率约束与锁定模板秩恢复，`climate_data.cli`将摘要与订正结果输出为JSON。所有观测值、站点、人群、情景与金额都是虚构。字段和研究约定见[资料约定](docs/data-contract.md)。

当前Python 3.11及以上可直接使用，没有第三方Python依赖。项目根目录运行：

```sh
python -m unittest discover -s tests -v
python -m compileall -q climate_data tests
python -m climate_data.cli inspect data/normal.json
python -m climate_data.cli correct data/normal.json
python -m climate_data.cli restore data/normal.json
```

`inspect`显示当前原始摘要，`correct`输出逐行订正（原值、midrank概率p、Q_OH/Q_MH两条参考分位、加法原始订正、截断损失、置零与抬升标记），`restore`并列给出恢复前后汇水区同期总量与各站变更账。所有命令只读取来源文件，不覆盖资料。`invalid-`前缀文件为不合法资料；模板缺站日时`restore`整次拒绝（返回码3）并说明原因。来源字节SHA-256随输出给出，结果同时标明输入身份、量纲和适用范围。

## 口径要点：保留了什么，没有保证什么

* **保留**：逐站变化信号以加法增量`Q_OH(p)+MF-Q_MH(p)`保留，不替换成比例缩放；零雨量不参与除法，不会出现旧乘法脚本的无限大。OH/MH只用历史期、MF只用未来期，未来观测不进入校准。
* **干日频率**：按`clip(f_OH+f_MF-f_MH,0,1)`显式约束，目标个数四舍五入、恰半向上；置零与湿日抬升分别记账。频率校准**不保证总雨量守恒**。
* **模板恢复**：恢复前后每站边际值多重集与零日数完全相同，全窗汇水区总量相同；跨站依赖只保证跟随事先锁定的模板秩结构，**不承诺精确相关系数**，模板必须锁定，不能按站挑年份或按未来观测挑相似段。

## 可运行例子

| 文件 | 覆盖的边界 |
| --- | --- |
| `data/normal.json` | 常规多站单季节，含置零与湿日抬升 |
| `data/all-dry.json` | 全干样本，midrank恒为0.5，全程零值 |
| `data/constant.json` | 各期常量，midrank并列、模板全并列 |
| `data/short-season.json` | 两个季节，winter仅1天（n=1分位路径） |
| `data/tied-template.json` | 模板值全并列，验证按模板日期决胜的确定性 |
| `data/incomplete-template.json` | 模板缺站日，恢复整次拒绝 |
| `data/invalid-domain.json` | wet_floor不高于dry_threshold等非法输入 |

容器使用固定基础镜像，可运行`docker build -t climate-material .`，随后`docker run --rm climate-material`检查资料。`docker run --rm climate-material python -m unittest discover -s tests -v`执行全部检查。运行不依赖远端账号或实时服务，没有气候预测服务或其他附件。
