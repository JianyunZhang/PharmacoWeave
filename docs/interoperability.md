# 跨项目互操作

0.0.4 用文件交换证据。它不安装、不导入 OmicsWeave、TheraWeave 或 EphysWeave。PharmacoWeave 自己的证据必须带 `derived_from`。

夹具在 `examples/cross_project_demo/`。

## 可以附带的内容

| 来源 | 示例 | 结果 |
|---|---|---|
| TheraWeave | `theraweave_predicted_admet.json` | 预测清除率进入 `external_evidence`，观测表和化合物表不变 |
| OmicsWeave | `omics_flux.json` | 通量估计作为证据引用 |
| EphysWeave | `ephys_endpoint.json` | 读出终点和 2.5 ms 同步误差，不改剂量时钟 |

附带使用：

```python
study.attach_candidates("examples/cross_project_demo/theraweave_predicted_admet.json")
```

## PharmacologyEvidence 往返

`pharmacoweave.integrations.evidence_from_nca` 把可报告的 NCA 参数写成 `evidence_type=estimated`。外推比例超过报告上限时不写出 AUC∞。`evidence_from_emax` 把 Emax 和 EC50 写成 estimated，把曲线上的效应写成 predicted。`evidence_from_adme` 使用 `permeability`、`microsomal_clearance` 或 `hepatocyte_clearance`，不用 `clearance`。未测量的体外行不产生数值证据。`evidence_from_tox` 的适用范围写明这不是临床给药建议。

`round_trip_exchange` 把交换文件转成 JSON 再读回。预测清除率不会变成 estimated 或 observed。`promote_to_observation` 仍然一律拒绝。

## 会被拒绝的内容

- 把 OmicsWeave 通量标成 `measurement_type=concentration`（`omics_mislabeled_concentration.json`）。
- 证据交换中的实测浓度。实测浓度来自研究观测表。
- TheraWeave 提供的 observed 或 estimated PK。
- 预测、文献或假设使用 `qc_status=accepted`。
- 缺少 `sync_error_ms` 或 `timebase_id` 的电生理终点（`ephys_missing_sync.json`）。程序不会在时钟基准缺失时插值。
- 把任意证据记录提升为 `Observation`。
- PharmacoWeave 证据缺少 `derived_from`。
- 把预测效应写成 `EffectObservation`。

## 身份

`identity_map` 记录 `source_project`、`source_study_id`、`source_subject_id`、`canonical_subject_id`、`match_method` 和 `review_status`。

合成例子里，`EW-RAT-A` 经 `synthetic_protocol` 审核后对应 `SUBJ-01`。这条链接存在，是因为方案写明了对应关系。文本相同不会自动生成链接。未审核链接在 research 模式下是警告，在 strict 模式下是错误。

结果包里的 `evidence/cross_project_links.json` 重复这条政策：本包不执行连接。

## 其他引擎

`supported_backends()` 返回 `pharmacoweave-analytical`。该后端是本包解析模拟器的 MIT 子进程，交换 JSON 并记录引擎版本、许可证和请求 SHA-256。`describe_engine("nlmixr2")` 与 `describe_engine("pk-sim")` 记录 copyleft 许可证。`run_backend("nlmixr2")` 和 `run_backend("pk-sim")` 拒绝执行，并且不会导入或启动这两个引擎。
