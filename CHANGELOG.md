# Changelog

## 0.2.6 — 2026-11-10

- **两阶段对数参数汇总（`summarize_log_parameters`）**：对调用方已经给出的个体 `CL` 与 `V`，或 `CL/F` 与 `V/F`，计算几何均值和对数尺度样本方差（分母 n−1）。标准误只保存，不作为权重。混合绝对参数与表观参数、少于两名受试者、同一参数内受试者重复、非正估计或量纲不符时，整次调用不给出均值。这不是 NLME 后验，也不是给药建议。
- **显式生理表（`validate_physiology_table`）**：体重、心输出量、组织体积、血流量、分配系数和每项来源说明都由调用方给出。流量换算到 L/h 后必须与心输出量一致。通过校验的记录没有浓度、AUC 或清除率。不查物种表，未列出的组织不补零。
- **线性周转闭式解（`simulate_linear_turnover`）**：`dR/dt = kin − kout·R`。证据等级是 `predicted`。`estimate_parameters=True`、非正速率或时间不严格递增时拒绝。曲线不连到临床终点。
- **肝加一个组织的稳态代数解（`steady_state_one_tissue`）**：`fu = 1`。`CL = Q_liver · CLint / (Q_liver + CLint)`。请求时间过程、组织不是正好两个、或没有名为 `liver` 的组织时拒绝。清除率名称保持调用方声明的 `CL` 或 `CL/F`。
- **调用方群体结果表（`load_population_result_table`）**：核对受试者、参数名、估计、单位、估计方法和证据等级。证据等级只保留 `estimated` 或 `predicted`。`source_parameter_name` 或 `source_evidence_grade` 与声明字段不一致时拒绝整张表。读入不创建子进程。
- **验收补齐**：混合 `V` 与 `V/F`、一组失败则整次无均值、标准误不改变几何均值、`mL/min` 与 `mL` 换到约定单位、流量相对容差 1e-6 与绝对容差 1e-9、缺少名为 `liver` 的组织、证据等级被改标、出现 `constitutes_dosing_recommendation`，以及负的 `R0`，都有失败用例或手算参考。
- **仍然拒绝**：`population_pk`、`pbpk`、`qsp` 的 `run_in_process` 继续拒绝进程内 FOCE、SAEM、MCMC、全身动态 PBPK 和 QSP 网络。这三个模块的执行状态改为 `partial`。nlmixr2 与 PK-Sim 仍不导入、不按名字启动。暴露模拟与引擎交换契约保持 `0.2.2`。新增五份结果契约，`schemas/` 共 19 份，研究契约仍是 `0.1.0`。

## 0.2.5 — 2026-11-10

- **调用方进程外引擎槽（`run_caller_engine`）**：可执行文件、版本字符串和许可证文本三者都在时，才启动子进程。结果记录这三项、可执行文件 SHA-256、请求 SHA-256 和响应 SHA-256。缺任何一项则不创建子进程。子进程输出按原样保存，不解释成 PharmacoWeave 的估计。命令是 `pharmacoweave engine-slot`。
- **内置后端不变**：`supported_backends()` 仍然只有 `pharmacoweave-analytical`。`run_backend("nlmixr2")` 与 `run_backend("pk-sim")` 在启动子进程之前拒绝。这两个引擎不进入依赖，包内也不重写 NONMEM 或 PK-Sim。

## 0.2.4 — 2026-11-10

- **尿药排泄 NCA（`run_urine_excretion`）**：单独配置 `analysis: urine_excretion`。区间排泄量是尿药浓度乘以记录的尿液体积。BLQ、超定量上限和缺失体积保持为空，不写成 0。只有全部区间都有排泄量且区间连续时才汇总总量。肾清除是该总量除以调用方给出的血浆 AUC，单位为 L/h，证据类型是 `renal_clearance`，不是血浆 `clearance`。
- **与血浆 NCA 分开**：尿药浓度不进入血浆浓度–时间轮廓，也不使用血浆梯形积分。结果写入 `analysis/urine_excretion.parquet`。采集区间和体积在可选表 `urine_collections`；没有尿药采集时不写出该表。命令是 `pharmacoweave urine-nca`。

## 0.2.3 — 2026-11-10

- **一室恒速输注拟合（`fit_one_compartment_iv_infusion`）**：在已有闭式解 `iv_infusion_concentration` 上，用对数参数的 Levenberg–Marquardt 估计 `CL` 与 `V`，并给出渐近标准误和 95% 置信区间。权重仍是 `uniform`、`1/y`、`1/y^2`。不新写 ODE 求解器。参数语义是静脉 `CL` 与 `V`，不是 `CL/F`。
- **失败时不输出清除率**：输注时长缺失或非正、时长不短于最后定量时间、输注后浓度不下降、或曲率不足以辨识两个参数时，诊断为失败，估计值为空。BLQ 不进入拟合，也不被写成 0。
- **命令**：`pharmacoweave fit --model one_compartment_iv_infusion` 使用加权非线性最小二乘。对数线性方法只保留给静脉推注。

## 0.2.2 — 2026-11-05

- **等剂量、等间隔 NCA**：`repeated_equal` 计算最后一个由下一次给药闭合的间隔；两端都有定量点才报告 `AUC_tau`，缺端点不外推。该设计不把 Dose/AUC_tau 写成清除率。`steady_state` 由调用方声明该间隔已达稳态，此时静脉清除率是 Dose/AUC_tau，血管外是 Dose/(CL/F)。间隔内不报告 AUC∞。剂量不等、间隔不等、或输注时长不短于间隔时拒绝。配置 `multiple_dose` 仍然无效。
- **对数参数相关矩阵**：模拟请求可以带 `parameter_correlation`（如 `clearance|volume`）。相关作用在对数参数尺度上，用 Cholesky 生成联合扰动。矩阵不正定或相关系数不在 (-1, 1) 时拒绝。未提供相关时，各参数仍独立抽样。这不是群体后验，也不是 NLME。

## 0.2.1 — 2026-11-02

- **重复给药与稳态解析暴露（`simulation/exposure.py`）**：在 0.2.0 的单次给药解析解上，用线性叠加实现等剂量、等间隔的重复给药，并用闭式解实现稳态区间。覆盖一室静脉推注、恒速输注、血管外吸收和二室静脉推注。有限次数的 AUC∞ 等于次数乘以单次 AUC；稳态不报告 AUC∞，而报告一个给药间隔的 `auc_tau`（静脉为 Dose/CL，血管外为 Dose/(CL/F)）。稳态时间必须落在一个给药间隔内。输注时长大于或等于给药间隔时拒绝。输出仍是 `predicted`，不构成给药方案。个体间协方差、PBPK 和 QSP 继续拒绝。
- **契约**：暴露模拟与引擎交换 Schema 升到 `0.2.1`。既有研究与结果包契约仍是 `0.1.0`。

## 0.2.0 — 2026-11-01

- **探索性单次给药暴露区间（`simulation/exposure.py`）**：在已经验收的解析解上模拟一室静脉推注、一室恒速输注、一室血管外吸收和二室静脉推注。输出保持 `evidence_type="predicted"`，`constitutes_dosing_recommendation` 固定为 false。血管外给药只报告 `CL/F` 与 `V/F`。`ka < kel` 时终末半衰期记在 `ka` 上；`ka = kel` 时使用 Bateman 极限形式。可选的独立对数正态参数扰动给出 5/50/95 分位，并写明这不是置信区间、预测区间或治疗窗。重复给药、稳态、群体、PBPK 和 QSP 继续拒绝。`dosing_recommendation` 与 `human_starting_dose_from_simulation` 直接拒绝。
- **进程外参考引擎（`integrations/external_backends.py`）**：唯一可执行后端是 `pharmacoweave-analytical`（MIT）。它通过子进程读写 JSON，记录引擎名、版本、许可证和请求 SHA-256。`nlmixr2`（GPL-2.0-or-later）和 `pk-sim`（GPL-2.0-only）只进入目录，不导入、不启动。`population_pk`、`pbpk` 和 `qsp` 的进程内计算仍然拒绝。
- **契约与命令**：新增 `pharmacoweave simulate` 与 `pharmacoweave engine-run`。结果包可写入 `analysis/exposure_simulation.json`。导出 `schemas/exposure_simulation.schema.json` 与 `schemas/engine_exchange.schema.json`（契约版本 `0.2.0`）。既有研究与结果包契约保持 `0.1.0`。

## 0.1.0 — 2026-10-31

- **研究目录观测向药效与毒理终点的显式桥接（`core/observation_bridge.py`）**：打通研究表 `observations`（`observation_type='effect'` 与 `'toxicity'`）向强类型 `EffectObservation` 与 `ToxEndpoint` 的显式桥接转换（`PharmacologyStudy.bridge_effects` 与 `bridge_toxicology`）；自适应解析关联 `Biosample`、`DoseEvent`、`Assay` 及给药时间差与生物基质，智能推导量纲（如 `%` / `percent` → `PERCENT`、`ratio` / `fraction` → `FRACTION`、浓度与速率量纲等），自动附加强制性 `RESEARCH_USE_LIMITATION` 与 `TOX_RESEARCH_LIMITATION`，并支持调用方细粒度参数重载（`endpoint_overrides`）。
- **跨尺度时序对齐与端到端可追溯性（`integrations/alignment.py`）**：交付跨项目完整可追溯示范流（`dose_event → in vivo exposure → EphysWeave / OmicsWeave PD endpoint → source artifact`）；对齐流程强制执行：① 跨项目受试者标识在 `identity_map` 中必须存在且通过审核（`review_status="accepted"`，未审核或拒绝时直接拒联）；② 声明的 `timebase_id` 必须在研究元数据时间基准表中已注册；③ 显式校准时钟同步抖动（`sync_error_ms`），超限（`sync_error_ms > max_allowed_sync_error_ms`）时拒绝关联；④ 严守负向科学红线：**代谢通量绝不是浓度**（Omics 反应通量试图声明为浓度或附带浓度量纲时直接触发 `ContractViolation` 拒算）。
- **声明式原始异构表映射工具（`ingestion/mapping.py`）**：支持通过 `ingest_mapping.yaml`（`IngestMappingConfig`）声明式映射异构湿实验导出表格（CSV/TSV/Parquet），实现列重命名、常量注入与枚举映射转换（如 `<LLOQ` / `BLQ` 自动映射为规范 `ObservationStatus.BLQ` 并抹除非实测数值），在零信息损失前提下将原始文件 SHA-256 摘要严格持久化至 `study.source_files` 确保证据链不可篡改，并提供 `PharmacologyStudy.from_mapping` 一键构建。
- **离线结果包核验与差异对比工具链（`core/bundle_verify.py` & `core/bundle_diff.py`）**：CLI 新增 `verify-bundle` 与 `diff` 命令；`verify-bundle` 离线校验 `manifest.json` 与 `schemas/result_bundle.schema.json` 契约，重新逐一哈希全量落盘文件比对 SHA-256，排查文件篡改、文件丢失与未列盘多余文件；`diff` 命令对比两次分析运行结果包的元数据差异、文件增删改状态、Parquet 分析参数表数值变化（自动计算最大绝对/相对偏差并标注显著差异参数）与 JSON 诊断指标变动，支持 `--json` 格式化机器可读输出。
- **契约与 Schema 升级（10 份落盘）**：版本号正式推进至 `0.1.0`；在 `schemas/` 中导出并冻结第 10 份数据契约 `ingest_mapping.schema.json`；提供端到端示范脚本 `examples/multiscale_ephys_omics_alignment.py` 与 `examples/ingest_mapping_demo.py`。

## 0.0.7 — 2026-10-31

- 扩展体外 ADME 记录与单位体系（`dmpk`）：`AdmeAssayKind` 新增血浆游离分数 `fu_plasma`、微粒体结合游离分数 `fu_mic` 与全血/血浆分配比 `blood_to_plasma_ratio`；单位系统引入 `fraction` 与 `ratio` 规范化映射，严格校验数值边界（`fu` 强制在 `(0, 1]` 范围且支持 `percent`/`pct` 转换，`blood_to_plasma_ratio > 0`），维持拒绝血浆清除率单位与未实测行绝不填充数值等安全防线。
- 显式生理假设的研究级探索性 IVIVE（`dmpk/ivive.py`）：仅在调用方显式提供特定物种生理参数（`PhysiologicalScalingParams`：体重、肝重、肝血流 $Q_H$、微粒体蛋白含量 MPPGL 或肝细胞密度 HPGL，核心库不包含隐式物种默认表）以及结合参数（`fu_plasma`、`fu_inc`、`blood_to_plasma_ratio`，标量输入显式记录为假设）时，按 Yang et al. (2007) Well-Stirred 肝清除率模型计算探索性预测（`IviveClearancePrediction`）。预测记录强制锁定为 `evidence_grade="predicted"`，`scale_to_in_vivo_clearance` 保持拒绝单条体外记录换算，`human_starting_dose_from_ivive` 坚决拒绝推导首次人体推荐剂量。
- 先验预测与后验体内估算折叠误差审计（`diagnostics/prediction_audit.py`）：对挂载的预测性先验证据（`TheraWeave` CADD、探索性 IVIVE 等 `evidence_type="predicted"`）与体内估算参数（`NcaParameterSet` 或房室拟合生成的 `estimated`）进行严格双向折叠误差校准审计（$P/O$、$\max(P/O, O/P)$、$\log_{10}$ 误差、2 倍/3 倍阈值标记与按 `measurement_type` 分组的 AFE、AAFE 诊断，防止跨参数正负偏差相消），保持输入记录严格只读不可变；严格语义隔离绝对清除率与表观清除率（$CL$ vs $CL/F$）、体积与表观体积（$V$ vs $V/F$）、通用参数与两室参数（排除 $Q$、$V2$）以及母药与代谢物，支持基于受试者体重的剂量/清除率单位换算与显式分子量摩尔转换。
- 命令行与结果包全链路落地：新增 `pharmacoweave ivive` 与 `pharmacoweave audit-predictions` CLI 子命令；结果包持久化 `analysis/ivive_predictions.parquet` 与 `analysis/prediction_audit.json`，HTML 验证报告新增探索性 IVIVE 与预测审计展示；导出 `schemas/prediction_audit.schema.json`。

## 0.0.6 — 2026-10-31

- 开放单次血管外给药（`oral`、`subcutaneous`、`intramuscular`）与单次静脉恒速输注（`iv_infusion`）NCA（`nca_extravascular`、`nca_iv_infusion`、`nca_profile`）：血管外给药从 `(0, 0)` 起点积分且不执行静脉 `C0` 反向外推，终末相候选窗口限制在严格达峰后点（`t > Tmax`），生物利用度 `F` 未知时强制报告表观清除率 `CL/F` 与表观终末分布容积 `Vz/F`；静脉恒速输注强制校验正 `infusion_duration`，从 `(0, 0)` 积分并将终末相限制在 `t >= max(Tmax, infusion_duration)`。同步冻结 `tests/regression/oral_quiet_nca.json` 与 `tests/regression/iv_infusion_quiet_nca.json` 两份 50 位十进制独立手算参考。
- 在 `pk_models` 中新增一室血管外一级吸收加权非线性拟合（`fit_one_compartment_extravascular`，估计 `ka`、`CL/F`、`V/F`）与二室静脉推注双指数加权非线性拟合（`fit_two_compartment_iv_bolus`，估计 `CL`、`V1`、`Q`、`V2`），提供渐近标准误与对数尺度 95% 置信区间，并开放 `weight="uniform"`；当吸收相观测缺失、`ka ≈ ke` 发生 flip-flop 不可辨识，或二室快慢相退化为单指数（`alpha / beta < 1.25`、快/慢相贡献过低或海森曲率比退化）时诊断为失败，不输出参数点估计。`pharmacoweave fit` 新增 `--model` 与 `--weight uniform` 选项。
- 开放显式分子量支撑的质量–摩尔量纲互转（`convert_with_molecular_weight`、`Compound.convert_quantity`），词表新增摩尔浓度原子（`fM`、`pM`、`nM`、`uM`、`mM`、`M`）、`fmol` 及分子量单位（`Da`、`kDa`、`g/mol`）；未提供正分子量或调用纯同量纲 `convert()` 时继续拒绝质量与摩尔互转。
- 扩展 `PharmacologyEvidence` 与跨项目证据守卫以承载表观药代参数（`apparent_clearance`、`apparent_volume`、`ka`），确保血管外 `CL_over_F` / `Vz_over_F` 绝不与静脉真实 `CL` / `V` 混同；同步更新 `schemas/` 导出的 JSON Schema 与 HTML/SVG 报告语义展示。

## 0.0.5 — 2026-10-30

- 直接效应与 Sigmoid Emax 默认改用原始效应尺度的阻尼 Levenberg–Marquardt 非线性最小二乘（`direct_emax_nls_v1` / `sigmoid_emax_nls_v1`），消除低浓度噪声放大以及带噪低浓度观测略低于固定基线 `E0` 时的拒算脆弱性；新增可选 Sigmoid Hill 系数估计（`fit_sigmoid_emax` / `estimate_hill`），同时保留 `method="linearized"` 供显式对比。高浓度平台覆盖不足、过渡区无观测、海森曲率比退化或效应方向矛盾时诊断为失败，不输出 `EC50` 或 `hill`。
- 基于收敛点 $J^T W J$ 协方差矩阵、Student-$t$ 分位数与正参数对数尺度 Delta 变换，为一室静脉推注加权非线性拟合（`CL`、`V`）与原始尺度 Emax/Sigmoid 拟合（`Emax`、`EC50`、`hill`）填充渐近标准误 `standard_error` 与 95% 置信区间（`ci_lower`、`ci_upper`、`ci_level=0.95`）。
- 新增纯标准库半对数浓度–时间图（`figures/concentration_time_semilog.svg`，含 LLOQ 阈值线、拟合曲线叠加与终末相选点标记）、加权相对残差图（`figures/weighted_residuals.svg`）及浓度–效应拟合曲线图（`figures/concentration_effect.svg`）；解除 `pharmacoweave nca --bundle` 对 `fit_iv_bolus` 的隐式耦合（默认纯 NCA，可通过 `--with-fit` 显式附加对数线性拟合）。
- 打通 Phase 2 全链路落盘与 Schema 导出：在结果包中持久化 `analysis/derived_observations.parquet`、`analysis/emax_parameters.parquet`、`analysis/emax_predictions.parquet`、`experimental/in_vitro_adme.parquet`、`experimental/tox_endpoints.parquet`，自动汇集并经证据守卫核验写入 `evidence/pharmacology_evidence.json`；导出 `pd_records.schema.json`、`dmpk_records.schema.json`、`toxicology_records.schema.json`；新增 `collapse_replicates` 可审计平行样聚合（算术/几何均值与 CV%）及 `pharmacoweave emax` CLI 子命令。

## 0.0.4 — 2026-10-30

- 效应观测记录时间尺度和量纲，并保持实测与预测分开。直接效应 Emax（基线固定、Hill 系数为 1）在无噪声合成曲线上回收 Emax 和 EC50。点数不足、效应平坦、效应等于基线，或方向与 Emax 符号相反时诊断为失败，不输出 EC50。AUC 和 Cmax 不能当作直接效应的浓度。BLQ 不替换成 0。滞后和间接反应不估计。
- 增加体外通透性、微粒体清除率和肝细胞清除率记录。血浆清除率单位会被拒绝。本版本不把这些结果换算成体内 CL。
- 增加研究用毒理终点。记录可以链接剂量或暴露，调用给药建议会失败。
- PharmacologyEvidence 用 `derived_from` 引用暴露或分析来源。预测与实测在 JSON 往返后等级不变。不导入 TheraWeave、OmicsWeave 或 EphysWeave。

## 0.0.3 — 2026-10-30

- 一室静脉推注拟合写出观测值–预测值和残差–时间 SVG。BLQ 不画成 0。`pharmacoweave nca --bundle` 和 `pharmacoweave fit --bundle` 把这两张图写入结果包的 `figures/`。
- 增加权重为 `1/y` 和 `1/y^2` 的一室静脉推注非线性最小二乘。带噪声合成例子在 10% 相对误差内回收 CL 和 V。浓度不下降或不可辨识时诊断为失败，不输出清除率。对数线性拟合仍然可用。
- 在 `tests/regression/iv_bolus_quiet_nca.json` 冻结无噪声静脉推注 NCA 的十进制参考值。相对误差上限写在 `docs/validation.md`。

## 0.0.2 — 2026-10-30

- 计算单次静脉推注、连续采血的 NCA：linear-up/log-down、对数线性回推 C0、adjusted R² 终末相。
- BLQ 和高于 ULOQ 的点不参与梯形积分，也不替换成数字。
- 外推比例超过配置上限时仍保存 AUC∞，并由 `auc_inf_is_reportable` 标成不宜报告。
- 用对数线性最小二乘拟合一室静脉推注，回收 CL 和 V。浓度不下降时不给出清除率。
- 结果包可以写入 `analysis/nca_parameters.parquet` 和浓度–时间图。`pharmacoweave nca` 运行该计算。

## 0.0.1 — 2026-10-29

研究用数据契约的第一版骨架。

- 增加研究、受试者、给药、样本、观测、检测和身份链接模型。
- 增加封闭单位词表、关系校验、可审计的 BLQ 策略和结果包。
- 增加无噪声一室静脉推注合成研究，以及一室静脉推注解析解。
- 增加 NCA 配置和结果 Schema。计算器尚未实现。
- 增加 OmicsWeave、TheraWeave、EphysWeave 的文件级证据守卫。
- 后期分析包可以导入，并标记为 reserved。
