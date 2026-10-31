# 数据契约

Schema 版本是 `0.1.0`。Pydantic 模型是唯一真源，`schemas/` 是导出的 10 份 JSON Schema。重新生成：

```bash
python -m pharmacoweave.schema_export
```

导出的 Schema 清单：

- `schemas/result_bundle.schema.json`（`ResultBundleManifest`）
- `schemas/evidence_exchange.schema.json`（`EvidenceExchange`）
- `schemas/pharmacology_study.schema.json`（`PharmacologyStudy`）
- `schemas/analysis_records.schema.json`（`AnalysisRecords`，含 PK 参数渐近标准误与置信区间字段）
- `schemas/nca_config.schema.json`（`NcaConfig`，支持 `iv_bolus`、`iv_infusion`、`oral`、`subcutaneous`、`intramuscular` 单次连续采血设计）
- `schemas/pd_records.schema.json`（`PdAnalysisRecords`，含 `EffectObservation`、`DirectEffectPair`、`EffectParameter`、`PredictedEffectRecord`）
- `schemas/dmpk_records.schema.json`（`InVitroAdmeSet`，含 `PhysiologicalScalingParams` 与 `IviveClearancePrediction`）
- `schemas/toxicology_records.schema.json`（`ToxEndpointSet`）
- `schemas/prediction_audit.schema.json`（`PredictionAuditReport`）
- `schemas/ingest_mapping.schema.json`（`IngestMappingConfig`，声明式异构实验室表格映射契约）

## 研究目录

| 文件 | 要求 |
|---|---|
| `study_metadata.json` | 必需 |
| `subjects.csv` / `.tsv` / `.parquet` | 必需，三者只能有一个 |
| `compounds` | 必需 |
| `analytes` | 必需 |
| `assays` | 必需 |
| `dosing_events` | 必需 |
| `biosamples` | 必需 |
| `observations` | 必需 |
| `identity_map` | 可省略 |

编码为 UTF-8，不带 BOM。空单元格是缺失，不是 0。`qc_flags` 用 `|` 连接，条目本身不能含 `|`。绝对时间是带时区偏移的 ISO-8601 字符串；Parquet 里同样用字符串，避免各语言的 timestamp 时区不一致。

列顺序以 `pharmacoweave.core.tables` 为准。多一列或少一列都会在读入时失败。

### 主要列

受试者：`subject_id`、`study_id`、`group_id`、`species`、`sex`、`body_weight_value`、`body_weight_unit`。

给药：`dose_event_id`、`subject_id`、`compound_id`、`absolute_time`、`dose_value`、`dose_unit`、`route`、`formulation`、`infusion_duration`、`infusion_duration_unit`、`body_weight_normalized_value`、`body_weight_normalized_unit`。

样本：`sample_id`、`subject_id`、`collection_time_absolute`、`time_since_dose`、`time_since_dose_unit`、`dose_event_id`、`matrix`、`collection_site`、`timebase_id`、`is_predose`。

观测：`observation_id`、`sample_id`、`analyte_id`、`observation_type`、`status`、`raw_value`、`reported_value`、`numeric_value`、`unit`、`assay_id`、`qc_flags`、`exclusion_reason`、`replicate_id`、`assay_batch_id`、`analyte_form`。

`observation_type` 取 `concentration`、`effect`、`toxicity`。`status` 取 `quantified`、`BLQ`、`above_ULOQ`、`missing`、`invalid`。只有 `quantified` 同时拥有 `numeric_value` 和 `unit`。

### Phase 2 领域记录

- **PD 效应与 Emax**：`EffectObservation` 强制携带 `time_scale`、`dimension`、`exposure_ref` 和 `evidence_grade="observed"`；`EffectParameter` 记录 `value`、`unit`、`fixed`、`estimation_method` 及可选 `standard_error`、`ci_lower`、`ci_upper`、`ci_level`；模型预测点固定为 `evidence_grade="predicted"`。
- **体外 ADME、蛋白结合与探索性 IVIVE**：`InVitroAdmeRecord` 支持通透性（`cm/s` 等）、微粒体/肝细胞清除率（`uL/min/mg_protein`、`uL/min/million_cells`）、血浆/微粒体游离分数（`fu_plasma`、`fu_mic`，无量纲，规范化在 `(0, 1]` 范围）及全血/血浆比（`blood_to_plasma_ratio > 0`）；`InVitroAdmeSet` 扩展容纳调用方显式声明的物种生理缩放参数（`PhysiologicalScalingParams`）与 Well-Stirred 探索性预测（`IviveClearancePrediction`，锁定 `predicted` 等级）。拒绝体内血浆清除率单位 `L/h`。
- **毒理终点**：`ToxEndpoint` 必须引用 `dose_event_id` 或 `exposure_ref`，并携带非临床给药建议声明。
- **预测校准审计**：`PredictionAuditReport` 结构化呈现对比条目 `entries`（比值 $P/O$、折叠误差、$\log_{10}$ 误差、2/3 倍标记）、未匹配条目 `skipped`（含具体原因）及按 `measurement_type` 分组的 `summaries`（AFE、AAFE、2/3 倍比例与异常值数）。

### Phase 3 领域记录与跨尺度契约

- **声明式异构表映射（`IngestMappingConfig`）**：`tables` 字典定义各研究表的 `source_file`、`column_mapping`、`constants`、`value_transforms` 及 `delimiter`；摄取时自动计算每个原始物理文件的 SHA-256 摘要并注入 `study.source_files`（`role="raw_<table_name>"`）。
- **跨尺度对齐终点（`AlignedEndpointRecord` & `MultiScaleAlignmentReport`）**：记录经时钟校准与受试者映射对齐后的电生理/组学终点，显式携带 `canonical_subject_id`、`dose_event_id`、`exposure_ref`、`timebase_id`、`sync_error_ms`，以及完整端到端可追溯路径 `traceability_path`（`dose -> exposure -> foreign PD endpoint -> source artifact`）。
- **结果包离线核验（`BundleVerificationResult`）与差异对比（`BundleDiffResult`）**：结构化输出结果包的哈希完整性诊断、Schema 依从性，以及跨运行两次结果包的元数据差异、文件树增删改状态与 Parquet 表数值相对误差。

## 标识符与单位

标识符匹配 `^[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}$`。封闭单位词表的基准是小时、毫克药物、千克体重、升和摩尔，支持 `ug/mL`、`mg*h/L`、`L/h`、`mg/kg`、`1/h`、摩尔浓度（`fM`–`M`）、分子量（`g/mol`、`Da`、`kDa`）、`percent`、`fold`、`fraction`、`ratio` 等显式组合，未知原子（如 `IU`）被拒绝。常规 `convert()` 仅在同一量纲内换算；质量与摩尔量纲之间的互转必须通过 `convert_with_molecular_weight()` 或带正分子量的 `Compound.convert_quantity()` 显式执行。

## 结果包产物

`manifest.json` 的 `schema_name` 是 `pharmacology-result-bundle`，`quality_status` 固定为 `research_only`。

- 实验观测表：`experimental/concentration_observations.parquet`、`experimental/pd_observations.parquet`、`experimental/toxicology_observations.parquet`（研究原始表保持只读），以及可选传入的 `experimental/in_vitro_adme.parquet`、`experimental/tox_endpoints.parquet`。
- 派生与分析表（仅在对应分析运行后写出）：`analysis/derived_observations.parquet`、`analysis/nca_parameters.parquet`、`analysis/model_parameters.parquet`、`analysis/predictions.parquet`、`analysis/residuals.parquet`、`analysis/emax_parameters.parquet`、`analysis/emax_predictions.parquet`、`analysis/ivive_predictions.parquet`、`analysis/diagnostics.json`、`analysis/prediction_audit.json`。
- 诊断图（原生 SVG）：`figures/concentration_time.svg`、`figures/concentration_time_semilog.svg`、`figures/goodness_of_fit.svg`、`figures/residual_diagnostics.svg`、`figures/weighted_residuals.svg`、`figures/concentration_effect.svg`。
- 证据交换：`evidence/pharmacology_evidence.json`（自动汇集并经 `screen_exchange` 核验所有运行产出的 `estimated`、`predicted` 与体外/毒理 `observed` 证据）。
