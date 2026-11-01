# Changelog

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
