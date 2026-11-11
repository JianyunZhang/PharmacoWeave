# 验证

0.0.6 在 0.0.5 统计加固与 Phase 2 全链路落盘之上，补齐了单次血管外给药与静脉恒速输注 NCA 的 50 位十进制独立手算冻结参考、一室血管外吸收与二室静脉推注加权非线性拟合基准（含不可辨识拒算门禁），以及显式分子量支撑的质量–摩尔量纲互转验证。

## 已自动执行的检查

- 单位维度和换算，含显式分子量下的质量–摩尔量纲互转与无分子量拒算。
- 单条记录规则：时区、BLQ、零浓度、输注时长、`CL` 与 `CL/F` 及 `Vz` 与 `Vz/F` 语义隔离、置信区间上下界顺序。
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

## 0.0.6 血管外/输注冻结参考、房室扩展与显式分子量换算

- **单次口服与静脉恒速输注冻结 NCA 参考**：`tests/regression/oral_quiet_nca.json`（`examples/data/oral_single_dose/`，两只大鼠，表观参数 `CL/F`、`Vz/F`）与 `tests/regression/iv_infusion_quiet_nca.json`（`examples/data/iv_infusion_single_dose/`，0.5 h 恒速静脉输注）均通过 50 位十进制独立手算生成。`test_oral_quiet_reference.py` 与 `test_iv_infusion_quiet_reference.py` 按与 `iv_bolus_quiet_nca.json` 相同的 1e-9 / 1e-12 误差上限严格锁死梯形积分、终末相选窗及表观/静脉语义字段。
- **一室血管外吸收拟合（`fit_one_compartment_extravascular`）**：`tests/numerical_benchmarks/test_extravascular_and_two_compartment.py` 核验无噪声与带噪合成曲线在相对误差上限内回收 `ka`、`CL/F`、`V/F` 及正值对数尺度 95% CI；缺少达峰前上升点、达峰后点不足或 `ka ≈ ke` 发生 flip-flop 混淆时强制诊断为 `fail`，不输出参数。
- **二室静脉推注双指数拟合（`fit_two_compartment_iv_bolus`）**：核验无噪声（`examples/data/two_compartment_iv_bolus/`，`CL = 0.25 L/h`、`V1 = 0.50 L`、`Q = 0.50 L/h`、`V2 = 1.0 L`）与带噪曲线回收四参数及 95% CI；单指数数据传入二室模型发生快慢相坍缩（`alpha / beta < 1.25` 或曲率退化）时强制诊断为 `fail`。
- **显式分子量换算**：`tests/unit/test_units.py` 校验 `convert_with_molecular_weight` 与 `Compound.convert_quantity` 在给定正分子量（`g/mol`、`Da`、`kDa`）时精确互转质量与摩尔浓度/剂量/AUC，并验证 `convert()` 或缺失分子量时持续拒绝跨量纲互转。

## 0.0.7 显式生理假设探索性 IVIVE 与预测校准审计

- **Well-Stirred IVIVE 手算解析基准**：`tests/numerical_benchmarks/test_ivive_and_prediction_audit.py` 针对微粒体与肝细胞 CLint，基于 Yang et al. (2007) 经典方程（包含培养非特异性结合校正 $\text{CL}_{\text{int,u,vitro}} = \text{CL}_{\text{int,app}} / f_{u,\text{inc}}$、全血游离分数 $f_{u,b} = f_{u,p} / R_{b/p}$、全血肝清除率 $\text{CL}_{H,b} = Q_H f_{u,b} \text{CL}_{\text{int,u,H}} / (Q_H + f_{u,b} \text{CL}_{\text{int,u,H}})$ 与血浆肝清除率 $\text{CL}_{H,p} = R_{b/p} \cdot \text{CL}_{H,b}$），手算验证相对误差不超过 **1e-6**。
- **高低清除率渐近极限验证**：在极低清除率极限下核验 $\text{CL}_{H,b} \to f_{u,b} \cdot \text{CL}_{\text{int,u,H}}$；在极高清除率极限下核验 $\text{CL}_{H,b} \to Q_H$ 且肝提取率 $E_H \to 1.0$；验证生理参数缺失 MPPGL/HPGL 或跨化合物借用蛋白结合记录时抛出 `AnalysisError`，调用 `human_starting_dose_from_ivive` 坚决拒绝推导人体起始剂量。
- **先验预测对体内估算折叠误差审计**：`audit_study_predictions` / `audit_evidence_predictions` 严格核验比值 $P/O$、对称折叠误差 $\max(P/O, O/P)$、$\log_{10}$ 误差、2 倍/3 倍阈值标记以及按 `measurement_type` 分组的 AFE（$10^{\frac{1}{N}\sum \log_{10}(P/O)}$）与 AAFE（$10^{\frac{1}{N}\sum |\log_{10}(P/O)|}$），验证在 2 倍高估与 2 倍低估并存时 AFE 正负偏差抵消为 1.0 而 AAFE 保持 2.0，防止跨参数偏差相消。
- **严格语义隔离与量纲自适应**：核验基于受试者体重的单位缩放（$\text{L/h/kg} \to \text{L/h}$）与基于化合物分子量的摩尔–质量换算；核验绝对与表观清除率（$CL$ vs $CL/F$）、通用与二室参数（排除 $Q$、$V2$）、母药与代谢物完全隔离不混同；核验输入证据只读不可变，非预测性先验或非估算性后验输入抛出 `ContractViolation`。
- **Schema 与结果包落盘**：`schemas/prediction_audit.schema.json` 导出并严格冻结，校验 `analysis/ivive_predictions.parquet`、`analysis/prediction_audit.json` 与 HTML 报告内嵌表格。

## 0.1.0 跨尺度时序对齐、声明式映射与结果包核验基准

- **跨尺度端到端可追溯基准**：`tests/numerical_benchmarks/test_multiscale_traceability.py` 核验给药事件（`DOSE-SUBJ-01`）到体内暴露（Cmax / AUC_inf）、再到 EphysWeave 神经表型（海马 CA1 theta 频段功率，$\text{sync\_error} = 8.45\text{ ms} \le 50.0\text{ ms}$）与 OmicsWeave 组学终点（Cyp3a4 表达 log2FC，$\text{sync\_error} = 32.10\text{ ms} \le 50.0\text{ ms}$）、最终溯源至原始文件 SHA-256 摘要的全链条精确对应（`traceability_path` 与 `derived_from` 四跳完整性）。
- **认识论与物理时钟硬性边界**：核验未审核或被拒绝的受试者身份链接（`review_status != ReviewStatus.ACCEPTED`）直接抛出 `AlignmentError`；时钟同步抖动超过声明阈值（$\text{sync\_error\_ms} > \text{max\_allowed\_sync\_error\_ms}$）直接抛出 `AlignmentError`；严守**代谢通量绝非浓度**底线，Omics 反应通量声明为浓度或附带浓度量纲时直接触发 `ContractViolation` 拒算。
- **声明式原始异构表映射基准**：`tests/unit/test_ingest_mapping.py` 校验由 `animals.csv`、`dosing.csv`、`pk_concentrations.csv` 与 `ingest_mapping.yaml` 构建的完整研究，严格核验原始实验文件 SHA-256 哈希完整保存在 `study.source_files`（`role="raw_<table_name>"`），且转换后的规范数据通过 `study.check()` 校验。
- **结果包离线哈希核验与差异对比基准**：`tests/unit/test_bundle_verify_and_diff.py` 校验 `verify_bundle` 对健康包返回 `is_valid=True`、对单字节篡改文件或丢失文件精准返回 `is_valid=False` 并指明破坏文件；核验 `diff_bundles` 对相同包返回 `identical=True`，对参数被修改 50% 的分析表返回 `identical=False` 并计算出精确相对误差（0.3333）。

## 0.2.0 探索性暴露区间与进程外参考引擎

- **闭式解析基准**：`tests/numerical_benchmarks/test_exposure_simulation.py` 锁定剂量 0.25 mg、CL = 0.05 L/h、V = 0.05 L 的一室静脉推注：Cmax = 5 mg/L，AUC∞ = 5 mg·h/L，t½ = ln(2) h，1 mg/L 以上时间为 ln(5) h。输注峰浓度对照 `iv_infusion_concentration`，相对误差不超过 **1e-12**。血管外曲线的 Cmax 对照 Bateman 函数，同样不超过 **1e-12**，并且参数名只有 `CL_over_F` 与 `V_over_F`。
- **终末相命名**：`ka < kel` 时 `terminal_rate_name="ka"`；`ka = kel` 时使用极限形式并记为 `ka_equals_kel`。
- **不确定度语义**：只扰动 V 时 AUC∞ 的 5/50/95 分位与点估计相同；只扰动 CL 时 Cmax 分位与点估计相同。分位带的说明文字写明它不是置信区间或治疗窗。
- **引擎边界**：`tests/unit/test_engine_adapter.py` 要求 `supported_backends()` 只有 `pharmacoweave-analytical`，子进程结果与进程内结果相对误差不超过 **1e-12**，并且 `nlmixr2` / `pk-sim` 在启动子进程之前被拒绝。源码不导入这些引擎。
- **给药建议**：`dosing_recommendation` 与 `human_starting_dose_from_simulation` 抛出 `ContractViolation`。

## 0.2.1 重复给药与稳态

- **重复静脉推注**：两剂、间隔 1 h、CL = 0.05 L/h、V = 0.05 L、剂量 0.25 mg 时，第二次给药瞬间浓度等于 `5 (1 + e^{-1})` mg/L，AUC∞ = 10 mg·h/L。阈值 3 mg/L 以上的时间等于两段 `ln(Cpeak / 3)` 之和。
- **稳态静脉推注**：间隔 2 h 时，`auc_tau` = 5 mg·h/L，AUC∞ 为空，Cmax = `5 / (1 - e^{-2})` mg/L，蓄积比等于 `1 / (1 - e^{-2})`。
- **稳态输注与二室**：峰浓度对照闭式解，相对误差不超过 **1e-12**。`auc_tau` 分别为 Dose/CL。
- **稳态血管外**：`auc_tau` = Dose/(CL/F)。t = 1 h 的浓度与 40 次单剂量 Bateman 叠加的相对误差不超过 **1e-8**。参数名保持 `CL_over_F`，不出现 `CL`。
- **拒绝**：缺少给药间隔、稳态同时声明有限次数、输注时长不短于给药间隔、稳态时间超出一个给药间隔，都会失败。

## 0.2.2 等剂量间隔 NCA 与对数参数相关

- **稳态静脉推注**：kel = 1 h⁻¹、τ = 2 h、剂量 0.25 mg、V = 0.05 L 的指数间隔上，`linear-up/log-down` 的 `AUC_tau` 相对误差不超过 **1e-12** 地等于 5 mg·h/L，清除率回到 0.05 L/h。AUC∞ 不报告。
- **未声明稳态的重复间隔**：同一单次曲线在 0 到 2 h 的 `AUC_tau` 等于 `5(1-e^{-2})`，清除率保持为空。
- **缺端点**：没有 τ 时刻的定量点时不报告 `AUC_tau`，也不外推。
- **剂量或间隔不等**：`equal_interval_profiles` 直接拒绝。
- **相关矩阵**：`clearance|volume = 0.8`、1500 次抽样时，对数清除率与对数体积的样本相关与 0.8 的绝对差不超过 0.05。相关系数取到 ±1，或 3 个参数的相关矩阵不正定，都会拒绝。

## 0.2.3 一室恒速输注拟合

- **闭式对照**：`tests/numerical_benchmarks/test_iv_infusion_fit.py` 用 `iv_infusion_concentration` 生成剂量 1 mg、输注 0.5 h、CL = 0.2 L/h、V = 1 L 的无噪声曲线。拟合回收的 CL 与 V 相对误差不超过 **1e-8**。报告的预测浓度对照同一函数在估计参数下的取值，相对误差不超过 **1e-12**。对数参数灵敏度对照该函数的中心差分，相对误差不超过 **1e-6**（绝对差不超过 **1e-8**）。
- **标准误**：无噪声曲线的渐近标准误绝对值不超过 **1e-8**，95% 置信区间与点估计的相对差不超过 **5e-8**。
- **失败时不输出清除率**：输注时长缺失或非正、时长不短于最后定量时间、定量点少于三个、浓度没有足够曲率，或输注结束后浓度不下降时，诊断为 `fail`，`estimates` 为空。参数语义是 `CL` 与 `V`，给药途径是 `iv_infusion`。
- **BLQ**：示例研究 `examples/data/iv_infusion_single_dose` 的给药前和 24 h BLQ 不进入拟合，也不被写成 0。该示例浓度以六位小数落盘，CL 与 V 的回收相对误差不超过 **1e-4**。

## 0.2.4 尿药排泄

- **手算对照**：`tests/regression/urine_excretion_quiet.json` 用十进制乘法写明 0–4 h、10 mg/L、0.25 L 的排泄量是 2.5 mg，4–8 h、4 mg/L、0.5 L 的排泄量是 2.0 mg，总量 4.5 mg。血浆 AUC 为 9 mg·h/L 时肾清除是 0.5 L/h。`run_urine_excretion` 与这些数值的相对误差不超过 **1e-12**。
- **不补零**：第二段为 BLQ、缺失体积或超定量上限时，该段排泄量为空，总量和肾清除也为空。理由写明该值没有被替换成 0。区间 0–4 h 与 6–8 h 之间有缺口时，两段排泄量保留，总量和肾清除不汇总。区间重叠时拒绝。
- **与血浆分开**：只含尿药的研究运行血浆 `run_nca` 得到空结果。尿药表没有 `auc_inf`。尿药模块源码不调用血浆梯形积分。

## 0.2.5 调用方进程外引擎槽

- **缺凭证不启动**：`tests/unit/test_caller_engine_slot.py` 把 `subprocess.run` 换成会失败的探针。可执行文件不存在、文件为空、版本为空白，或许可证文本为空白时，都抛出 `AnalysisError`，且探针不被调用。
- **完整凭证才记录**：可执行文件、版本和许可证文本都存在时，结果保存这三项，以及可执行文件、请求和响应的 SHA-256。响应文本与子进程写出的字节一致，不解释成 PharmacoWeave 的估计。`constitutes_dosing_recommendation` 为 false。子进程退出码非 0 时不保存这条运行。
- **内置后端不变**：同一测试要求 `supported_backends()` 仍是 `pharmacoweave-analytical`。`run_backend("nlmixr2")` 与 `run_backend("pk-sim")` 在启动子进程之前拒绝。既有 MIT 参考引擎对照见上文 0.2.0。

## 0.2.6 群体汇总、生理表、周转模型和单组织稳态

- **两阶段汇总**：`tests/regression/two_stage_log_parameters.json` 给出 0.05 与 0.20 L/h。几何均值是 0.10 L/h。对数方差按 `[(ln 0.05 − ln 0.10)^2 + (ln 0.20 − ln 0.10)^2] / 1` 手算。`summarize_log_parameters` 与该手算的相对误差不超过 **1e-12**。混合 `CL` 与 `CL/F`、混合 `V` 与 `V/F`、少于两名受试者、同一参数内受试者重复、非正估计、缺估计或错误量纲时拒绝。其中一组失败时，整次调用不返回任何均值。标准误可以随行保存；1e-6 与 50 的标准误不改变几何均值或对数方差。
- **单位**：`tests/regression/two_stage_canonical_units.json` 把 1 与 4 mL/min 写成 0.06 与 0.24 L/h，几何均值 0.12 L/h；把 500 与 2000 mL 写成 0.5 与 2 L，几何均值 1 L。对数方差是 `2·(ln 2)^2`。相对误差不超过 **1e-12**。
- **生理表**：流量总和不等于心输出量时拒绝。体积为 0、分配系数缺失，或来源说明短于 3 个字符时，记录无法建成。通过的 `PhysiologyTableCheck` 没有浓度、AUC 或清除率字段。重复组织名拒绝。`tests/regression/physiology_flow_tolerance.json` 把平衡容差定为相对 1e-6、绝对 1e-9：落在容差内的流量通过，超出容差的流量拒绝。相对容差用心输出量 2 L/h 检查，绝对容差用心输出量 1e-6 L/h 检查。
- **周转模型**：`tests/regression/linear_turnover.json` 取 kin = 2、kout = 1、R0 = 0。稳态是 2。`R(0) = 0`。`R(1)` 对照 `(kin/kout)·(1−e^{−kout})`，相对误差不超过 **1e-12**。非正速率、负的 R0、时间不递增，或 `estimate_parameters=True` 时拒绝。
- **单组织稳态**：`tests/regression/steady_state_one_tissue.json` 取肝血流 1 L/h、另一组织血流 1 L/h、心输出量 2 L/h、CLint = 1 L/h、输入速率 1 mg/h、Kp_liver = 1、Kp_tissue = 2。肝清除率 0.5 L/h，动脉 2 mg/L，肝 1 mg/L，另一组织 4 mg/L。相对误差不超过 **1e-12**。`time_course=True`、组织不是正好两个，或没有名为 `liver` 的组织时拒绝。
- **群体结果表**：缺列、把 `CL/F` 标成 `CL`、证据等级为 `observed`、把 `predicted` 改标成 `estimated`，或出现 `dosing_recommendation` / `constitutes_dosing_recommendation` 时拒绝。读取过程把 `subprocess.run` 换成会失败的探针，探针不被调用。
- **重型算法仍拒绝**：`population_pk`、`pbpk`、`qsp` 的 `run_in_process` 仍抛出 `PlannedCapabilityError`。这三个包的源码不导入 nlmixr2、rpy2、PK-Sim、NONMEM 或 MoBi。

## 0.2.7 加权两阶段汇总

- **相对标准误相同**：`tests/regression/weighted_two_stage.json` 取 0.05 与 0.20 L/h，标准误为 0.005 与 0.02 L/h。`s / theta` 相同，几何均值是 0.10 L/h，对数方差是 `2·(ln 2)^2`。与 `summarize_log_parameters` 的相对误差不超过 **1e-12**。两名标准误再乘同一个正数，结果不变。绝对标准误相同而点估计不同时，权重不同，几何均值不再是 0.10 L/h。
- **两名受试者、标准误不同**：0.05 L/h 配 0.005 L/h，0.20 L/h 配 0.10 L/h。权重是 100 与 4。几何均值是 `exp((100·ln 0.05 + 4·ln 0.20) / 104)`。对数方差是 `(ln 0.20 − ln 0.05)^2 / 2`。相对误差不超过 **1e-12**。
- **三名受试者**：0.05、0.10、0.20 L/h 配标准误 0.05、0.10、0.10 L/h。权重是 1、1、4。几何均值是 `0.1·√2`，对数方差是 `(7/6)·(ln 2)^2`。相对误差不超过 **1e-12**。该结果与不加权汇总不同。
- **单位**：50 与 200 mL/h，标准误 10 与 40 mL/h，换到 L/h 后与 0.05、0.20 L/h 的等权重结果一致。
- **不退回**：标准误缺失、为 0 或为负时拒绝，且不调用不加权公式代替。其中一组失败时不返回其他组的均值。

## 0.2.8 周转抑制、效应室、三室和已验收记录

- **恒定抑制**：`tests/regression/inhibited_turnover.json` 取 kin = 2、kout = 1、R0 = 0、I = 0.5。稳态是 1，`R(ln 2) = 0.5`，相对误差不超过 **1e-12**。`I = 0` 时与 `simulate_linear_turnover` 的相对误差不超过 **1e-12**。`I = 1` 或 `estimate_parameters=True` 时拒绝。
- **效应室**：`tests/regression/effect_compartment.json` 取 Dose = 1 mg、V = 1 L、CL = 1 L/h、ke0 = 2 h⁻¹。`Ce(0) = 0`，`Ce(ln 2) = 0.5 mg/L`。`ke0 = kel = 1` 时 `Ce(1) = 1/e mg/L`。相对误差不超过 **1e-12**。血浆浓度来自一室静脉推注闭式解。
- **三室**：`tests/regression/three_compartment_iv_bolus.json` 取 Dose = 1 mg、CL = 0.5 L/h、V1 = 1 L、Q2 = 1.5 L/h、V2 = 4 L、Q3 = 0.25 L/h、V3 = 10 L。`C(0) = 1 mg/L`。系数积分 AUC∞ = 2 mg·h/L，相对误差不超过 **1e-12**。无噪声回收和 `fit_three_compartment_iv_bolus_profile` 的相对误差不超过 **1e-4**。拟合报告标准误和 95% 置信区间。点数不足时诊断为失败，估计值为空。
- **命令与结果包**：`accepted-record` 写出的 JSON 与直接调用的相对误差不超过 **1e-12**。结果包文件是 `analysis/accepted_record.json`。`kind` 为 `engine_response` 时拒绝。群体结果表请求不调用 `subprocess.run`。

## 0.2.10 浓度驱动周转、三室单次模拟和审计拟合

- **浓度驱动周转**：`tests/regression/concentration_turnover.json`。`kin = kout = kel = 1`、`R0 = 0` 时 `R(1) = 1/e`。`kin = kout = 2`、`kel = 1` 时，`R(ln 2)` 的驱动项是 0.5。相对误差不超过 **1e-12**。血浆浓度来自一室静脉推注闭式解。`estimate_parameters=True` 时拒绝。
- **三室单次模拟**：`tests/regression/three_compartment_exposure.json`。`C(0) = 1 mg/L`，`AUC∞ = 2 mg·h/L`，相对误差不超过 **1e-12**。终末速率名是 `gamma`。`n_doses = 2` 时请求无法建成。
- **审计拟合选项**：三室拟合导出的 `CL` 证据与 0.5 L/h 的相对误差不超过 **1e-4**。`audit-predictions --help` 列出 `three_compartment_iv_bolus`。

## 明确不做的符合性声明

字段设计参考 ICH M10 的定量限和 QC 概念，以及 ICH M15 的 question of interest 与 context of use。测试通过不等于生物分析方法学验证，也不等于 MIDD 证据已经被监管机构接受。CDISC SEND 导出、21 CFR Part 11 审计追踪和人类受试者访问控制不在本版本。
