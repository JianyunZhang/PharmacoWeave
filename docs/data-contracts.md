# 数据契约

Schema 版本是 `0.1.0`。Pydantic 模型是来源，`schemas/` 是导出的 JSON Schema。重新生成：

```bash
python -m pharmacoweave.schema_export
```

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

## 标识符

标识符匹配 `^[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}$`。研究、受试者、检测、给药、样本和观测在各自表内唯一。

## 单位

封闭词表的基准是小时、毫克药物、千克体重、升和摩尔。可组合，例如 `ug/mL`、`mg*h/L`、`L/h`、`mg/kg`、`1/h`。`kg` 是体重。未知原子，例如 `foot` 或 `IU`，会被拒绝。

## 结果包

`manifest.json` 的 `schema_name` 是 `pharmacology-result-bundle`，`quality_status` 固定为 `research_only`。`data` 列出实际写出的相对路径。`omitted` 列出故意没有写出的分析产物。`files` 带 SHA-256 和媒体类型。

实验观测按类型分开：

- `experimental/concentration_observations.parquet`
- `experimental/pd_observations.parquet`
- `experimental/toxicology_observations.parquet`

三张表使用同一观测列契约。空表仍然带列名。

证据文件使用 `pharmacology-evidence-exchange`。它引用外部结论，不把外部结论复制进观测表。

## 分析记录

`schemas/analysis_records.schema.json` 描述模型说明、拟合、参数、NCA 参数集和诊断。分析包把它们写到 `analysis/analysis_records.json`，并把 NCA 行写到 `analysis/nca_parameters.parquet`。只做校验的包不写这些表，并在 `omitted` 中说明原因。报告规则包括：终末相未合格时不能带 λz、t½、AUC∞；静脉与 CL/F 不能同时出现；血管外给药与 CL、V 不能同时出现。Tmax 可以为 0。Schema 版本仍是 `0.1.0`。

`schemas/nca_config.schema.json` 只接受单次静脉推注、连续采血和 `linear_up_log_down`。
