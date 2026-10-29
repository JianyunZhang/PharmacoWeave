# 验证

0.1.0 验证数据契约、锁定设计的 NCA，以及一室静脉推注拟合。通过这些测试不是生物分析验证，也不是监管接受。

## 已自动执行的检查

- 单位维度和换算。
- 单条记录规则：时区、BLQ、零浓度、输注时长、CL 与 CL/F。
- 关系规则：外键、时钟、LLOQ、基质、未审核身份链接、临床数据的研究用提示。
- BLQ 策略不修改原观测；有偏替换必须 `acknowledge_bias=true`。
- 合成研究在 strict 模式下零错误、零警告。
- CSV、TSV、Parquet 往返后数值和校验结果保持一致。
- 只做校验的结果包：文件哈希、`research_only`，以及 `analysis_executed=false`。
- 分析包在 NCA 与拟合写出后 `analysis_executed=true`，并包含浓度–时间图。
- 证据夹具符合证据交换 Schema；错误升级会被拒绝。
- 一室静脉推注解析解，以及下面的 NCA 与拟合靶值。

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

`tests/numerical_benchmarks/test_iv_bolus_truth.py` 锁定这些解析靶值。NCA 使用量化点上的梯形积分，不调用 Dose/CL 闭式解。24 h 保持 BLQ。

## NCA 与拟合的数值容差

规则和冻结数字写在 `tests/numerical_benchmarks/README.md`。`reference_linlog.py` 是独立梯形转录，不导入本包，也不是 PKNCA 或 Phoenix。PKNCA 与 Phoenix 的选点规则彼此不同，因此对照的是这份事先冻结的参考，而不是把其中一个库装进依赖。

| 比较 | 量 | 接受范围 |
|---|---|---|
| 库与冻结梯形参考 | AUC_last、AUC∞、外推比例、Cmax、Tmax、λz、t½ | 绝对 `1e-9` |
| 库与解析靶值 | AUC∞ 相对 Dose/CL = 5 | 绝对 `1e-4` |
| 库与解析靶值 | λz 相对 0.2 h⁻¹ | 绝对 `1e-5` |
| 库与解析靶值 | t½ 相对 ln(2)/0.2 | 绝对 `1e-4` |
| NCA 的 CL、Vz 与生成值 | 受试者 1：0.05 L/h、0.25 L；受试者 2：0.06 L/h、0.30 L | 相对 `1e-4` |
| 无噪声拟合 | 同上的 CL 与 V | 相对 `1e-4` |
| 比例误差 σ = 0.02，种子 20261029 | CL 与 V | 相对 `0.10` |

六点小数舍入使解析窗口宽于梯形窗口。这份数据上的 AUC∞ 约为 5.0000026，测试要求它不等于 5。把 12 h 到 24 h 用零浓度补一条梯形会使 AUC_last 大约增加 0.54，计算器不得得到这个数。

终末相使用 Tmax 之后的定量点，对 ln(浓度) 做普通最小二乘。从最少 3 点到全部可用点中，取调整 R² 最大的窗口。调整 R² 相差不超过 `1e-6` 时取点数更多的窗口。调整 R² 低于 0.90 或点数不足时，省略 λz、t½ 和 AUC∞。外推比例高于 0.20 时保留已合格的 λz 和 t½，省略 AUC∞、CL 和 Vz。

失败或省略还覆盖：内部 BLQ 在 `keep_as_censored` 下不插补、重复时间、高于 ULOQ、剂量与浓度量纲不一致，以及口服、多次给药和稀疏设计。

拟合在少于两个定量时间、优化器不收敛、参数贴上边界（含停在边界内侧 0.1% 以内）或雅可比秩亏时不得标成 pass，并且不写出估计。BLQ 不进入目标函数；这不是删失似然。NCA 的 Vz 与房室模型的 V 是两个量。

## 明确不做的符合性声明

字段设计参考 ICH M10 的定量限和 QC 概念，以及 ICH M15 的 question of interest 与 context of use。测试通过不等于生物分析方法学验证，也不等于 MIDD 证据已经被监管机构接受。CDISC SEND 导出、21 CFR Part 11 审计追踪和人类受试者访问控制不在本版本。
