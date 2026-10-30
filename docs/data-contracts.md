# 数据契约

Schema 版本是 `0.1.0`。Pydantic 模型是唯一真源，`schemas/` 是导出的 8 份 JSON Schema。重新生成：

```bash
python -m pharmacoweave.schema_export
```

导出的 Schema 清单：

- `schemas/result_bundle.schema.json`（`ResultBundleManifest`）
- `schemas/evidence_exchange.schema.json`（`EvidenceExchange`）
- `schemas/pharmacology_study.schema.json`（`PharmacologyStudy`）
- `schemas/analysis_records.schema.json`（`AnalysisRecords`，含 PK 参数渐近标准误与置信区间字段）
- `schemas/nca_config.schema.json`（`NcaConfig`）
- `schemas/pd_records.schema.json`（`PdAnalysisRecords`，含 `EffectObservation`、`DirectEffectPair`、`EffectParameter`、`PredictedEffectRecord`）
- `schemas/dmpk_records.schema.json`（`InVitroAdmeSet`）
- `schemas/toxicology_records.schema.json`（`ToxEndpointSet`）

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
- **体外 ADME**：`InVitroAdmeRecord` 允许 `cm/s`、`um/s`、`nm/s`、`uL/min/mg_protein`、`mL/min/mg_protein`、`uL/min/million_cells` 和 `mL/min/million_cells`，拒绝体内血浆清除率单位 `L/h`。
- **毒理终点**：`ToxEndpoint` 必须引用 `dose_event_id` 或 `exposure_ref`，并携带非临床给药建议声明。

## 标识符与单位

标识符匹配 `^[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}$`。封闭单位词表的基准是小时、毫克药物、千克体重、升和摩尔，支持 `ug/mL`、`mg*h/L`、`L/h`、`mg/kg`、`1/h`、`percent`、`fold` 等显式组合，未知原子（如 `IU`）被拒绝。

## 结果包产物

`manifest.json` 的 `schema_name` 是 `pharmacology-result-bundle`，`quality_status` 固定为 `research_only`。

- 实验观测表：`experimental/concentration_observations.parquet`、`experimental/pd_observations.parquet`、`experimental/toxicology_observations.parquet`（研究原始表保持只读），以及可选传入的 `experimental/in_vitro_adme.parquet`、`experimental/tox_endpoints.parquet`。
- 派生与分析表（仅在对应分析运行后写出）：`analysis/derived_observations.parquet`、`analysis/nca_parameters.parquet`、`analysis/model_parameters.parquet`、`analysis/predictions.parquet`、`analysis/residuals.parquet`、`analysis/emax_parameters.parquet`、`analysis/emax_predictions.parquet`、`analysis/diagnostics.json`。
- 诊断图（原生 SVG）：`figures/concentration_time.svg`、`figures/concentration_time_semilog.svg`、`figures/goodness_of_fit.svg`、`figures/residual_diagnostics.svg`、`figures/weighted_residuals.svg`、`figures/concentration_effect.svg`。
- 证据交换：`evidence/pharmacology_evidence.json`（自动汇集并经 `screen_exchange` 核验所有运行产出的 `estimated`、`predicted` 与体外/毒理 `observed` 证据）。
