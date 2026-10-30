# 验证

0.0.5 在静脉推注 NCA、一室拟合、直接效应 Emax、体外 ADME、毒理记录和证据往返之上，补充了原始效应尺度阻尼非线性 Emax / Sigmoid Hill 拟合、参数渐近标准误与 95% 置信区间、可审计平行样聚合与派生观测表，以及半对数/加权残差/浓度–效应诊断 SVG 与 Phase 2 结果包全链路落盘验证。

## 已自动执行的检查

- 单位维度和换算。
- 单条记录规则：时区、BLQ、零浓度、输注时长、CL 与 CL/F、置信区间上下界顺序。
- 关系规则：外键、时钟、LLOQ、基质、未审核身份链接、临床数据的研究用提示。
- BLQ 策略与派生观测表不修改原观测；有偏替换必须 `acknowledge_bias=true`；平行样聚合计算 CV% 且不把 BLQ 平均为 0。
- 合成研究在 strict 模式下零错误、零警告。
- CSV、TSV、Parquet 往返后数值和校验结果保持一致。
- 结果包的文件哈希、`research_only`、分析缺席声明、Phase 2 落盘表及自动生成的 `PharmacologyEvidence`。
- 全部 8 份 JSON Schema 与 Pydantic 模型逐字段同步；证据夹具符合证据交换 Schema，越权升级会被拒绝。

本地与 CI 运行：

```bash
pip install -e ".[dev]"
python -m pharmacoweave.schema_export
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

`tests/numerical_benchmarks/test_iv_bolus_truth.py` 锁定这些靶值。

## 0.0.2 已经核对的数值

无噪声一室曲线（CL = 0.05 L/h，V = 0.25 L，剂量 0.25 mg）上：

- AUC∞ 等于 Dose/CL = 5 mg·h/L
- λz 等于 0.2 1/h，t½ 等于 ln(2)/0.2 h
- CL 和 Vz 回到 0.05 L/h 和 0.25 L
- 只有两个定量点时不报告 t½ 和 AUC∞
- 采样很早时 AUC∞ 仍可计算，外推比例高于 0.20，`auc_inf_is_reportable` 为 false
- 合成例子里的 24 h BLQ 被省略，AUC_last 对应积分到 12 h
- 对数线性拟合在同一曲线上回收 CL 和 V；浓度随时间上升时不给出清除率

这些检查在 `tests/numerical_benchmarks/test_nca_iv_bolus.py`。终末相选窗规则固定为 `pharmacoweave.lambda_z.adjusted_r2.v1`。

## 0.0.3 冻结参考和加权拟合

`tests/regression/iv_bolus_quiet_nca.json` 用 50 位十进制运算独立计算无噪声例子的梯形 NCA。`tests/regression/test_iv_bolus_quiet_reference.py` 要求相对误差不超过下表：

| 量 | 误差上限 |
|---|---|
| AUC_last、AUC∞、外推比例、λz、t½、CL、Vz、终末相 adjusted R² | 相对 1e-9 |
| Cmax | 相对 1e-12 |
| Tmax | 绝对 1e-12 h |

`examples/data/iv_bolus_noisy/`（剂量 1 mg，CL = 0.40 L/h，V = 2.0 L）要求 `1/y` 和 `1/y^2` 加权非线性最小二乘都在 **相对误差 10%** 以内回收 CL 和 V。

## 0.0.4 直接效应、体外 ADME、毒理和证据

`examples/direct_effect_and_evidence.py`（E0 = 10 percent，Emax = 80 percent，EC50 = 8 mg/L，Hill = 1）在 `tests/numerical_benchmarks/test_direct_emax.py` 中要求无噪声 Emax 和 EC50 的相对误差不超过 **1e-8**，预测点 `evidence_grade` 为 `predicted`。体外记录拒绝血浆清除率单位 `L/h`，`scale_to_in_vivo_clearance` 与 `dosing_recommendation` 始终拒绝执行。

## 0.0.5 原始尺度带噪 Emax、Sigmoid Hill、标准误与诊断图

`examples/noisy_emax_case.py` 与 `tests/numerical_benchmarks/test_noisy_emax.py` 锁定原始效应尺度阻尼 Levenberg–Marquardt 求解器的验收标准：

- **带噪直接效应 Emax（含低浓度跨越基线波动与 BLQ 删失点）**：低浓度点 `9.45 percent < E0 = 10.0 percent` 时，传统 `method="linearized"` 因双倒数符号翻转而失败；原始尺度 `nonlinear_ls` 正常收敛，在 **15%** 相对误差内回收真值（Emax = 80 percent，EC50 = 8 mg/L），且 95% 置信区间覆盖真值。
- **无噪声与带噪 Sigmoid Emax（`fit_sigmoid_emax` / `estimate_hill`）**：无噪声曲线（Hill = 2.0）在 **1e-6** 相对误差内精确回收 `Emax`、`EC50` 和 `hill`；带噪曲线在 **20%–25%** 相对误差内稳健回收三参数并给出对数尺度正值 95% 置信区间。
- **退化与不可辨识拒算保护**：仅覆盖线性区（`EC50` 超出最高测试浓度 1.25 倍）、全部落在高浓度平台、阶跃断层缺少 10%–90% 过渡区观测、海森曲率比退化、或 Sigmoid 定量点少于 4 个时，诊断强制为 `fail`，不输出 `EC50` 或 `hill`。
- **PK WNLS 渐近不确定度**：`tests/numerical_benchmarks/test_wnls_iv_bolus.py` 校验 `ParameterEstimate` 的 `standard_error`、`ci_lower`、`ci_upper` 和 `ci_level == 0.95`，在带噪合成曲线上 95% CI 覆盖真值 `CL = 0.40 L/h` 与 `V = 2.0 L`。
- **原生 SVG 诊断图**：`tests/unit/test_fit_figures.py` 通过 XML 解析器校验 `concentration_time_semilog.svg`、`weighted_residuals.svg` 和 `concentration_effect.svg`，确认 LLOQ 线、终末相标记和 Emax 曲线渲染正常且绝不把 BLQ 画成 0。

## 明确不做的符合性声明

字段设计参考 ICH M10 的定量限和 QC 概念，以及 ICH M15 的 question of interest 与 context of use。测试通过不等于生物分析方法学验证，也不等于 MIDD 证据已经被监管机构接受。CDISC SEND 导出、21 CFR Part 11 审计追踪和人类受试者访问控制不在本版本。
