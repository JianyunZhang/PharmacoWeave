# 贡献

PharmacoWeave 是研究用软件。拉取请求需要让科学约定保持可测试。

- 新的数值方法先增加基准：真值、容差、失败输入，以及 BLQ 或单位错误时的期望行为。
- 不要把 BLQ 写成 0，不要把 CL/F 写成静脉清除率，不要把预测 ADMET 写进浓度表。
- 不要导入 nlmixr2、PK-Sim、MoBi、NONMEM 或其他 Weave 项目。外部引擎走进程外适配器。
- `core` 不依赖分析包。
- Schema 以 Pydantic 模型为准。改模型后运行 `python -m pharmacoweave.schema_export`，并提交 `schemas/`。
- 提交前运行 `ruff check src tests examples`、`mypy src/pharmacoweave` 和 `pytest`。

例子数据必须是合成数据或可合法再分发的数据。不要加入可识别的人类受试者记录。
