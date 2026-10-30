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
- 直接效应 Emax：基线固定，Hill 系数固定为 1。`E = E0 + Emax * C / (EC50 + C)`。浓度–效应对不足、效应不随浓度变化、效应等于基线，或声明方向与 Emax 符号相反时，诊断为失败，不输出 EC50。AUC 和 Cmax 不能当作这个模型的浓度。滞后和间接反应不估计。
- 体外 ADME 记录：通透性（`cm/s`、`um/s`、`nm/s`）、微粒体清除率（`uL/min/mg_protein`）和肝细胞清除率（`uL/min/million_cells`）。`L/h` 会被拒绝。不做体内外推。
- 毒理终点链接剂量事件或暴露 id。不能生成临床给药建议。
- `PharmacologyEvidence` 用 `derived_from` 指向剂量、暴露或拟合。预测 ADMET 在往返后仍是 `predicted`。不导入 TheraWeave。

## 当前基线（0.0.4）模块状态与边界

- `available`：`core`、`ingestion`、`bioanalysis`、`pk_nca`（仅单次静脉推注连续采血）。
- `partial`：`pk_models`、`diagnostics`、`reporting`、`integrations`、`pd_models`、`pkpd`、`dmpk`、`toxicology`。
- `reserved`：`population_pk`、`pbpk`、`qsp`、`simulation`。
- 当前仍显式拒绝：稀疏/破坏性采样、稳态、尿药、静脉输注、血管外给药、多次给药、二室模型、Sigmoid Hill 估计、间接反应、效应滞后、血浆蛋白结合、体内外推（IVIVE）、SEND 导出和临床给药建议。

## 计划：统计加固与 Phase 2 全链路收口（0.0.5）

- **直接效应 Emax 原始尺度非线性拟合与可选 Hill 估计**：将 `pd_models` 中的双倒数线性化替换为原始效应空间的阻尼非线性最小二乘，消除低浓度噪声放大以及带噪观测跨越基线 `E0` 时的拒算脆弱性；新增可选 Sigmoid Emax（估计 `hill` 系数）。高浓度平台覆盖不足、海森曲率比退化或效应方向矛盾时诊断为失败，不输出 `EC50` 与 `hill`。先冻结带噪声合成基准与退化失败用例再开放接口。
- **参数渐近标准误与置信区间**：利用收敛点 $J^T W J$ 协方差矩阵，为一室加权非线性拟合与直接效应 Emax 填充 `ParameterEstimate` 的 `standard_error`、`ci_lower` 和 `ci_upper`。信息矩阵奇异或自由度不足时保持空值并记录诊断，不伪造有限区间。
- **半对数图、加权残差图与 CLI 解耦**：新增半对数浓度–时间 SVG（标示 LLOQ 阈值线、拟合曲线与终末相选点窗口）、加权/相对残差–时间 SVG 及浓度–效应拟合 SVG；解除 `pharmacoweave nca --bundle` 对 `fit_iv_bolus` 的隐式耦合。
- **Phase 2 结果包落盘、Schema 导出与派生层收口**：将 `EffectObservation`、`DirectEmaxFit`、`InVitroAdmeRecord`、`ToxEndpoint` 及新生成的 `PharmacologyEvidence` 写入 `PharmacologyResultBundle`、同步至 `schemas/` 并补齐 CLI 入口；在 `core/profiles.py` 串联 `bioanalysis` 的可审计 BLQ 派生策略与 `replicate_id` 平行样聚合，落盘派生分析表。

## 计划：血管外给药、恒速输注与二室模型（0.0.6，Phase 1 余量）

- **单次血管外给药（`oral` / `sc`）与静脉恒速输注（`iv_infusion`）NCA**：血管外给药以 `(0, 0)` 作为起点积分，不执行静脉 `C0` 反向外推，生物利用度 `F` 未知时强制报告 `CL/F` 与 `Vz/F`；静脉恒速输注核验显式 `infusion_duration`。两类设计均先建立 50 位十进制独立手算冻结参考 JSON。
- **一室血管外吸收拟合与二室静脉推注拟合**：在 `pk_models` 中增加一室一级吸收拟合（`ka`、`CL/F`、`V/F`，当吸收相点数不足或 `ka ≈ ke` 发生 flip-flop 不可辨识时诊断为失败）及二室静脉推注双指数拟合（`CL`、`V1`、`Q`、`V2`）。
- **显式分子量浓度换算**：仅在 `Compound.molecular_weight` 已提供并校验通过时，开放质量浓度与摩尔浓度之间的显式换算；缺失分子量时继续拒绝自动互转。
- **仍拒绝范围**：稀疏/破坏性采样复合曲线、多次给药叠加、稳态给药与尿药排泄仍显式拒绝。

## 计划：蛋白结合、显式假设 IVIVE 与预测偏差审计（0.0.7，Phase 2 余量）

- **蛋白结合与全血分配记录**：在 `dmpk` 中新增血浆游离分数（`fu_plasma`）、微粒体非特异性结合分数（`fu_mic`）与全血/血浆浓度比（`blood_to_plasma_ratio`）记录，校验数值区间与 QC 状态。
- **显式生理参数的研究级 IVIVE**：仅在调用方显式传入物种生理缩放参数（如微粒体蛋白含量、肝细胞密度、肝血流量）与游离分数时，计算探索性预测清除率；输出强制标记 `evidence_type="predicted"` 并附带模型假设与 `limitations`，不混入实测 `CL`，也不给出人体起始剂量。
- **先验预测与后验实测校准审计**：对比挂载的外部先验预测（来自 `TheraWeave` 或外部 CADD 导出的 `evidence_type="predicted"` 记录）与体内估计参数（`estimated`），输出折叠误差（Fold Error）、平均绝对折叠误差（AAFE）与平均折叠误差（AFE）诊断表，原始记录保持只读。

## 计划：Phase 3 组学与电生理终点对齐（0.1.0）

- **研究目录观测与跨尺度可追溯示例**：建立研究表 `observations`（`effect`、`toxicity`）到 `EffectObservation` 与 `ToxEndpoint` 的显式桥接；在 `identity_map` 审核通过（`review_status="accepted"`）且同步误差存在的前提下，交付端到端可追溯例子：`dose_event → exposure → EphysWeave / OmicsWeave PD endpoint → source artifact`。对齐必须显式校验 `timebase_id` 与 `sync_error_ms`；通量仍然不是浓度，同步误差超限时拒绝关联。
- **声明式原始表映射与结果包核验工具**：支持可选 `ingest_mapping.yaml`，在保留原始实验文件 SHA-256 摘要的前提下映射异构实验室导出表；CLI 增加离线结果包哈希与 Schema 复核（`verify-bundle`）及两次运行结果包差异对比（`diff`）。

## 计划：Phase 4 外部 PopPK、PBPK、QSP 与研究设计模拟

- **进程外引擎适配器**：通过独立子进程与文件交换协议接入经过选择的外部引擎（如 `nlmixr2`、`Open Systems Pharmacology / PK-Sim & MoBi`），不把 copyleft 引擎引进 `pharmacoweave` 的运行时 `import` 图。
- **科学验收与模拟边界**：验收包含适配器契约测试、数值交叉验证、外部引擎名称/版本/许可证溯源记录，以及强制写明的 `context_of_use` 和 `limitations`。本项目不实现第二套 PK-Sim 或 NONMEM；研究设计模拟（`simulation`）仅输出显式模型假设下的探索性暴露区间，不构成临床给药方案。

## 持续约束

- 新算法先有基准和失败用例，再出现在公开函数里。
- 负结果、删失、排除原因和局限性留在输出里。
- 人类临床数据保持研究用标识，例子里不放可识别的受试者信息。
- 跨项目协作保持文件级交换，零运行时耦合；预测与实测、静脉参数与表观参数保持类型隔离。
