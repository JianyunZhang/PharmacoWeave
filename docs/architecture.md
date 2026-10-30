# 架构

PharmacoWeave 0.0.5 是一个研究用的定量药理学数据、统计估计与证据工具包。计算模块依赖 `core`；`core` 在模块顶层不依赖任何分析业务包。集成模块导出证据并做文件级守卫，不导入其他 Weave 项目或外部 GPL/专有求解器。

```mermaid
flowchart TB
    subgraph available [可执行]
        ingestion[ingestion]
        bioanalysis[bioanalysis BLQ 与派生表]
        core[core 契约/单位/平行样/结果包]
        nca[pk_nca 单次静脉推注]
    end
    subgraph partial [契约与局部估计器]
        models[pk_models OLS / WNLS + CI]
        pd[pd_models 原始尺度直接/Sigmoid Emax]
        pkpd[pkpd 显式浓度连接]
        dmpk[dmpk 体外 ADME 记录]
        tox[toxicology 研究用终点]
        diagnostics[diagnostics]
        reporting[reporting HTML + 5 类原生 SVG]
        integrations[integrations 证据汇集与守卫]
    end
    reserved[population_pk / pbpk / qsp / simulation]

    ingestion --> core
    bioanalysis --> core
    nca --> core
    models --> core
    pd --> core
    pkpd --> pd
    dmpk --> core
    tox --> core
    diagnostics --> core
    reporting --> core
    integrations --> core
    reserved -.-> core
```

## 运行路径

1. `PharmacologyStudy.from_files` 读取 `study_metadata.json` 以及 CSV、TSV 或 Parquet 表。
2. 单条记录的规则在 Pydantic 模型里执行：时区、剂量维度、BLQ 不得带数值、输注必须有时长、正参数置信区间必须为正。
3. `study.check` 执行跨记录规则：外键、时钟、LLOQ、基质、身份链接和外部证据。
4. `collapse_replicates` 可选对同一给药后时间点的定量浓度平行样执行算术或几何均值聚合，记录 CV% 与溯源 ID，绝不把 BLQ 纳入均值或改为 0。
5. `run_nca`、`fit_iv_bolus` / `fit_iv_bolus_wnls`、`fit_linked_direct_emax` / `fit_sigmoid_emax` 在各自契约边界内执行计算。
6. `export_bundle` 在校验通过后写出新的不可变结果目录，自动物化 `analysis/derived_observations.parquet`，汇集并核验 `evidence/pharmacology_evidence.json`。目录已存在时拒绝覆盖。

## 记录规则

- 研究、受试者、化合物、分析物、检测、给药、样本、观测是分开的类型。
- 受试者是统计单位。技术重复用 `replicate_id`，不是新的受试者。
- 剂量的主字段是实际给药量（药物质量或物质的量）。`mg/kg` 写在归一化字段里，并与体重核对。
- 样本同时保存绝对时间和相对给药时间，并声明 `timebase_id`。负的相对时间必须标成 predose。
- 浓度观测必须指向分析物和浓度型检测。检测物基质与样本基质不一致时校验失败。
- 模型说明必须写下 question of interest 和 context of use。诊断或拟合要标成 pass，必须同时写下局限性，且 `converged=true`。

## 结果包

写出的目录包含研究表、实验表、派生观测表、分析参数与残差表、SVG 诊断图、来源、软件版本、运行时、证据、身份链接、校验 JSON 和 HTML 报告。每个文件都有 SHA-256，`manifest.json` 自身的摘要写在 `manifest.sha256`。

未运行分析时 `analysis_executed` 为 `false`，未执行的表与图列在 `omitted` 并说明缺席原因。`pharmacoweave nca --bundle` 默认只执行纯 NCA（写入 `nca_parameters.parquet`、`derived_observations.parquet`、线性与半对数浓度–时间图及 NCA 证据），不隐式运行房室拟合；需要同时拟合并输出残差图时可显式传入 `--with-fit` 或使用 `pharmacoweave fit --bundle`。`pharmacoweave emax --study ... --bundle ...` 写出 Emax 参数表、预测表和浓度–效应拟合图。

## 证据边界

外部或本运行生成的 JSON 均通过 `EvidenceExchange` 与 `screen_exchange` 守卫：

- 实测浓度只存在于观测表。证据交换里的 `observed` + `concentration` 会被拒绝。
- OmicsWeave 的通量、丰度或通路不能标成浓度、清除率、AUC、半衰期或分布容积。
- TheraWeave 可以附带 `predicted` ADMET，不能提供 `observed` 或 `estimated` 的 PK。
- 预测、文献和假设不能使用 `qc_status=accepted`。
- EphysWeave 的神经或行为终点必须带 `subject_id`、`session_id`、`timebase_id` 和 `sync_error_ms`。读取函数不移动时钟，也不做插值。
- `promote_to_observation` 一律拒绝。PharmacoWeave 写出的证据必须带非空 `derived_from`。模型预测效应始终保持 `predicted`。

## 依赖方向

允许：`ingestion`、`bioanalysis`、`pk_nca`、`pk_models`、`pd_models`、`dmpk`、`toxicology`、`diagnostics`、`reporting` 依赖 `core`。`pkpd` 依赖 `pd_models`。`integrations` 依赖 `core` 和这些记录的公开类型，用来写出证据。

不允许：`core` 在模块顶层依赖上述业务包；任何模块导入 OmicsWeave、TheraWeave、EphysWeave、PKNCA、nlmixr2、PK-Sim 或 NONMEM。
