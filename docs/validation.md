# 验证

0.0.2 验证数据契约、单次静脉推注 NCA，以及一室静脉推注的对数线性拟合。

## 已自动执行的检查

- 单位维度和换算。
- 单条记录规则：时区、BLQ、零浓度、输注时长、CL 与 CL/F。
- 关系规则：外键、时钟、LLOQ、基质、未审核身份链接、临床数据的研究用提示。
- BLQ 策略不修改原观测；有偏替换必须 `acknowledge_bias=true`。
- 合成研究在 strict 模式下零错误、零警告。
- CSV、TSV、Parquet 往返后数值和校验结果保持一致。
- 结果包的文件哈希、`research_only`、以及 NCA 产物缺席。
- 证据夹具符合证据交换 Schema；错误升级会被拒绝。
- 一室静脉推注解析解和未来 NCA 的靶值。

本地与 CI 运行：

```bash
pip install -e ".[dev]"
ruff check src tests examples
mypy src/pharmacoweave
pytest
```

CI 在 Python 3.11 和 3.12 上执行 lint、类型检查和测试。

## 合成基准

`examples/data/iv_bolus_single_dose/`：

| 量 | 值 |
|---|---|
| 给药 | 1 mg/kg 静脉推注 |
| CL | 0.20 L/h/kg |
| V | 1.0 L/kg |
| 受试者体重 | 0.25 kg 与 0.30 kg |
| C(0) | 1 µg/mL |
| t½ | ln(2) / 0.2 h |
| AUC∞ | 5 µg·h/mL |
| LLOQ | 0.01 µg/mL |
| 24 h | 低于 LLOQ，状态 BLQ，无数值 |

`tests/numerical_benchmarks/test_iv_bolus_truth.py` 锁定这些靶值。将来的 NCA 应在量化点上、并把 24 h 保持为删失之后，接近 Dose/CL 和 ln(2)·V/CL。允许误差要在实现该计算器时单独写明，不能用本文件的解析解代替梯形积分测试。

## 0.0.2 已经核对的数值

无噪声一室曲线（CL = 0.05 L/h，V = 0.25 L，剂量 0.25 mg）上：

- AUC∞ 等于 Dose/CL = 5 mg·h/L
- λz 等于 0.2 1/h，t½ 等于 ln(2)/0.2 h
- CL 和 Vz 回到 0.05 L/h 和 0.25 L
- 只有两个定量点时不报告 t½ 和 AUC∞
- 采样很早时 AUC∞ 仍可计算，外推比例高于 0.20，`auc_inf_is_reportable` 为 false
- 合成例子里的 24 h BLQ 被省略，AUC_last 对应积分到 12 h
- 对数线性拟合在同一曲线上回收 CL 和 V；浓度随时间上升时不给出清除率

这些检查在 `tests/numerical_benchmarks/test_nca_iv_bolus.py`。它们用解析真值，而不是把 PKNCA 或 Phoenix 引进依赖。和某一外部实现逐行对照仍是后续工作；终末相选窗规则已经固定为 `pharmacoweave.lambda_z.adjusted_r2.v1`。

## 明确不做的符合性声明

字段设计参考 ICH M10 的定量限和 QC 概念，以及 ICH M15 的 question of interest 与 context of use。测试通过不等于生物分析方法学验证，也不等于 MIDD 证据已经被监管机构接受。CDISC SEND 导出、21 CFR Part 11 审计追踪和人类受试者访问控制不在本版本。
