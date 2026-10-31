# 路线

阶段按科学验收切片，不按日历。

## 已完成：Phase 0 骨架（0.0.1）

- 可安装的 `pharmacoweave` 包、CI、MIT 许可和引用文件。
- 研究实体、封闭单位、关系校验、BLQ 策略。
- 合成静脉推注例子和严格校验。
- 结果包、HTML 校验清单、JSON Schema。
- 跨项目证据守卫和 NCA 配置拒绝超范围设计。
- 一室静脉推注解析解。

## 已完成：第一种 NCA 和一室静脉拟合（0.0.2）

- 单次静脉推注、连续采血、linear-up/log-down。
- 终末相规则 `pharmacoweave.lambda_z.adjusted_r2.v1`。点数不足或 adjusted R² 不够时不报告 λz、t½、AUC∞。
- 外推比例超过配置上限时保留数值，并由 `auc_inf_is_reportable` 返回 false。
- 对数线性拟合回收 CL 与 V。斜率不下降时诊断为失败，不输出清除率。
- 结果包在传入 NCA 结果后写出参数表和浓度–时间图，`analysis_executed` 为 true。

## 已完成：残差图、加权拟合和冻结 NCA 参考（0.0.3）

- 一室静脉推注拟合写出观测值–预测值和残差–时间 SVG。BLQ 不进入图，也不画成 0。
- 加权非线性最小二乘（权重 `1/y` 和 `1/y^2`）在带噪声合成曲线上回收 CL 与 V。浓度不下降或参数不可辨识时诊断为失败，不输出清除率。对数线性拟合仍然可用。
- `tests/regression/iv_bolus_quiet_nca.json` 是无噪声静脉推注 NCA 的第二份冻结参考。相对误差上限写在 `docs/validation.md`。这份数字用十进制独立算出，不把 PKNCA 或 Phoenix 引进依赖。

## 已完成：直接效应、体外 ADME、毒理记录和证据往返（0.0.4）

- 效应观测带时间尺度（`direct`、`interval`、`time_to_event`）和量纲（`fraction`、`percent`、`concentration`、`rate`、`time`）。它引用一条暴露，不把预测写成实测。
- 直接效应 Emax：基线固定，Hill 系数固定为 1。浓度–效应对不足、效应不随浓度变化、或声明方向与 Emax 符号相反时，诊断为失败，不输出 EC50。AUC 和 Cmax 不能当作这个模型的浓度。滞后和间接反应不估计。
- 体外 ADME 记录：通透性（`cm/s`、`um/s`、`nm/s`）、微粒体清除率（`uL/min/mg_protein`）和肝细胞清除率（`uL/min/million_cells`）。`L/h` 会被拒绝。不做体内外推。
- 毒理终点链接剂量事件或暴露 id。不能生成临床给药建议。
- `PharmacologyEvidence` 用 `derived_from` 指向剂量、暴露或拟合。预测 ADMET 在往返后仍是 `predicted`。不导入 TheraWeave。

## 已完成：统计加固与 Phase 2 全链路收口（0.0.5）

- **直接效应与 Sigmoid Emax 原始尺度非线性拟合**：将默认估计器升级为原始效应尺度的阻尼 Levenberg–Marquardt 非线性最小二乘（`direct_emax_nls_v1` / `sigmoid_emax_nls_v1`），消除双倒数线性化的低浓度噪声放大以及带噪低浓度观测略低于基线 `E0` 时的拒算脆弱性；开放可选 Sigmoid Hill 系数估计（`fit_sigmoid_emax` / `estimate_hill`），同时保留 `method="linearized"` 供对比。高浓度平台未覆盖、过渡区缺失观测、海森曲率比退化或效应方向矛盾时诊断为失败，不输出 `EC50` 与 `hill`。
- **参数渐近标准误与 95% 置信区间**：基于收敛点 $J^T W J$ 协方差矩阵、Student-$t$ 分位数与正参数对数尺度 Delta 变换，为一室加权非线性拟合（`CL`、`V`）与直接/Sigmoid Emax（`Emax`、`EC50`、`hill`）填充 `standard_error`、`ci_lower`、`ci_upper` 和 `ci_level=0.95`；自由度不足或海森矩阵奇异时保持空值并记录警告。
- **半对数图、加权残差图、浓度–效应图与 CLI 解耦**：新增半对数浓度–时间 SVG（标示 LLOQ 虚线、拟合叠加曲线与终末相选点框）、加权相对残差–时间 SVG 及浓度–效应拟合 SVG；解除 `pharmacoweave nca --bundle` 对 `fit_iv_bolus` 的默认隐式耦合（可通过 `--with-fit` 显式开启）。
- **Phase 2 结果包落盘、Schema 导出与派生/平行样收口**：结果包新增 `analysis/derived_observations.parquet`、`analysis/emax_parameters.parquet`、`analysis/emax_predictions.parquet`、`experimental/in_vitro_adme.parquet`、`experimental/tox_endpoints.parquet`，并自动汇总核验写入 `evidence/pharmacology_evidence.json`；同步生成 `pd_records.schema.json`、`dmpk_records.schema.json`、`toxicology_records.schema.json`；在 `core/profiles.py` 提供 `collapse_replicates` 可审计算术/几何均值聚合与 CV% 统计，并补齐 `pharmacoweave emax` CLI 子命令。

## 已完成：血管外给药、恒速输注与二室模型（0.0.6，Phase 1 余量收口）

- **单次血管外给药（`oral` / `subcutaneous` / `intramuscular`）与静脉恒速输注（`iv_infusion`）NCA**：血管外给药以 `(0, 0)` 作为起点积分，不执行静脉 `C0` 反向外推，终末相候选窗口限制在 `t > Tmax`，生物利用度 `F` 未知时强制报告表观参数 `CL/F`（`CL_over_F`）与 `Vz/F`（`Vz_over_F`）；静脉恒速输注强制核验正 `infusion_duration`，终末相限制在 `t >= max(Tmax, infusion_duration)`。两类设计均在 `tests/regression/oral_quiet_nca.json` 与 `tests/regression/iv_infusion_quiet_nca.json` 建立 50 位十进制独立手算冻结参考。
- **一室血管外吸收拟合与二室静脉推注拟合**：在 `pk_models` 中新增一室一级吸收加权非线性拟合（`ka`、`CL/F`、`V/F`）及二室静脉推注双指数加权非线性拟合（`CL`、`V1`、`Q`、`V2`），填充渐近标准误与对数尺度 95% 置信区间；吸收相缺失、`ka ≈ ke` 发生 flip-flop 不可辨识，或二室快慢相坍缩为单指数时强制诊断为 `fail` 并撤回参数估计。
- **显式分子量浓度与剂量换算**：词表新增摩尔浓度（`fM`–`M`）、`fmol` 与分子量（`Da`、`kDa`）原子，仅在显式传入正 `Compound.molecular_weight` 时通过 `convert_with_molecular_weight` / `Compound.convert_quantity` 开放质量与摩尔量纲互转；缺失分子量或调用常规 `convert()` 时继续拒绝自动互转。

## 已完成：蛋白结合、显式假设 IVIVE 与预测偏差审计（0.0.7，Phase 2 余量收口）

- **蛋白结合与全血/血浆分配比记录**：在 `dmpk` 中新增血浆游离分数（`fu_plasma`）、微粒体非特异性结合分数（`fu_mic`）与全血/血浆浓度比（`blood_to_plasma_ratio`）记录；单位词表补充 `fraction` 与 `ratio` 映射，校验无量纲数值边界（`fu` 规范化为 `(0, 1]`，`blood_to_plasma_ratio > 0`），维持拒绝血浆清除率单位与未实测行无数值等安全规则。
- **显式生理假设的研究级探索性 IVIVE**：在 `dmpk/ivive.py` 中实现 Well-Stirred 肝清除率预测（`predict_hepatic_clearance_well_stirred`）；仅在调用方显式传入物种生理缩放参数（`PhysiologicalScalingParams`：体重、肝重、肝血流 $Q_H$、微粒体蛋白含量 MPPGL 或肝细胞密度 HPGL，核心库绝不包含隐式物种默认表）以及结合参数（`fu_plasma`、`fu_inc`、`blood_to_plasma_ratio`，标量输入显式记录为假设）时完成推导。输出强制标记 `evidence_grade="predicted"` 并附带显式假设与局限性，不混入实测/估算体内 `CL`，`scale_to_in_vivo_clearance` 维持拒绝单记录缩放，`human_starting_dose_from_ivive` 坚决拒绝给出首次人体或安全给药推荐。
- **先验预测与后验实测校准审计**：在 `diagnostics/prediction_audit.py` 中实现先验预测与体内估算参数的双向折叠误差校准审计（`audit_study_predictions` / `audit_evidence_predictions`）；输出比值 $P/O$、对称折叠误差 $\max(P/O, O/P)$、$\log_{10}$ 误差、2 倍/3 倍阈值标记以及按 `measurement_type` 分组的 AFE、AAFE 诊断，防止跨参数正负偏差相互抵消；严格语义隔离绝对清除率与表观清除率（$CL$ vs $CL/F$）、体积与表观体积（$V$ vs $V/F$）、通用参数与两室参数（排除 $Q$、$V2$）以及母药与代谢物，输入证据保持严格只读不可变。
- **命令行与结果包全链路贯通**：新增 `pharmacoweave ivive` 与 `pharmacoweave audit-predictions` CLI 子命令；结果包持久化 `analysis/ivive_predictions.parquet` 与 `analysis/prediction_audit.json`，HTML 验证报告呈现探索性 IVIVE 与预测审计诊断表格；导出并冻结 `schemas/prediction_audit.schema.json`。

## 已完成：Phase 3 组学与电生理终点对齐、声明式映射与结果包核验（0.1.0）

- **研究目录观测向药效与毒理终点的显式桥接**：建立研究表 `observations`（`effect`、`toxicity`）到 `EffectObservation` 与 `ToxEndpoint` 的显式桥接转换（`PharmacologyStudy.bridge_effects` / `bridge_toxicology`）；自适应解析关联 `Biosample`、`DoseEvent`、`Assay` 及给药时间差与生物基质，智能推导量纲，自动附加强制性 `RESEARCH_USE_LIMITATION` 与 `TOX_RESEARCH_LIMITATION`，并支持细粒度参数重载。
- **跨尺度时序对齐与端到端可追溯示范流**：在 `integrations/alignment.py` 中交付跨项目完整可追溯链条：`dose_event → exposure → EphysWeave / OmicsWeave PD endpoint → source artifact`。对齐流程强制执行：① 跨项目受试者标识在 `identity_map` 中必须存在且通过审核（`review_status="accepted"`，未审核或拒绝时直接拒绝关联）；② 声明的 `timebase_id` 必须在研究元数据时间基准表中已注册；③ 显式校准时钟同步抖动（`sync_error_ms`），超限（`sync_error_ms > max_allowed_sync_error_ms`）时拒绝关联；④ 严守负向科学红线：**代谢通量绝不是浓度**（Omics 反应通量试图声明为浓度或附带浓度量纲时直接触发 `ContractViolation` 拒算）。
- **声明式异构实验室表格映射（`ingest_mapping.yaml`）**：在 `ingestion/mapping.py` 中实现通过 YAML 配置映射异构湿实验导出表格（CSV/TSV/Parquet），实现列重命名、常量注入与枚举映射转换（如 `<LLOQ` / `BLQ` 自动映射为规范 `ObservationStatus.BLQ` 并抹除非实测数值），在零信息损失前提下将原始文件 SHA-256 摘要严格持久化至 `study.source_files` 确保证据链不可篡改，并提供 `PharmacologyStudy.from_mapping`。
- **离线结果包核验与差异对比工具链（`verify-bundle` & `diff`）**：CLI 新增 `verify-bundle` 与 `diff` 命令；`verify-bundle` 离线校验 `manifest.json` 与 `schemas/result_bundle.schema.json` 契约，重新逐一哈希全量落盘文件比对 SHA-256，排查文件篡改、文件丢失与未列盘多余文件；`diff` 命令对比两次分析运行结果包的元数据差异、文件增删改状态、Parquet 分析参数表数值变化（自动计算最大绝对/相对偏差并标注显著差异参数）与 JSON 诊断指标变动，支持 `--json` 格式化机器可读输出。
- **10 份落盘 JSON Schema**：在 `schemas/` 中导出并冻结第 10 份数据契约 `ingest_mapping.schema.json`。

## 当前基线（0.1.0）模块状态与边界

- `available`：`core`、`ingestion`、`bioanalysis`、`pk_nca`、`integrations`。
- `partial`：`pk_models`（一室静脉推注、一室血管外吸收、二室静脉推注）、`diagnostics`（SVG 图形与先验预测校准审计）、`reporting`、`pd_models`（原始尺度直接与 Sigmoid Emax、研究观测桥接）、`pkpd`、`dmpk`（体外记录、蛋白结合与显式假设 IVIVE）、`toxicology`（毒理终点记录与研究观测桥接）。
- `reserved`：`population_pk`、`pbpk`、`qsp`、`simulation`。
- 当前仍显式拒绝：稀疏/破坏性采样、多次给药叠加、稳态给药、尿药排泄、间接反应、效应滞后、隐式物种生理查找表、非显式假设体内外推、SEND 导出和临床给药建议。

## 计划：Phase 4 外部 PopPK、PBPK、QSP 与研究设计模拟

- **进程外引擎适配器**：通过独立子进程与文件交换协议接入经过选择的外部引擎（如 `nlmixr2`、`Open Systems Pharmacology / PK-Sim & MoBi`），不把 copyleft 引擎引进 `pharmacoweave` 的运行时 `import` 图。
- **科学验收与模拟边界**：验收包含适配器契约测试、数值交叉验证、外部引擎名称/版本/许可证溯源记录，以及强制写明的 `context_of_use` 和 `limitations`。本项目不实现第二套 PK-Sim 或 NONMEM；研究设计模拟（`simulation`）仅输出显式模型假设下的探索性暴露区间，不构成临床给药方案。

## 持续约束

- 新算法先有基准和失败用例，再出现在公开函数里。
- 负结果、删失、排除原因和局限性留在输出里。
- 人类临床数据保持研究用标识，例子里不放可识别的受试者信息。
- 跨项目协作保持文件级交换，零运行时耦合；预测与实测、静脉参数与表观参数保持类型隔离。
