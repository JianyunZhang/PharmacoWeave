# 跨项目夹具

这些 JSON 文件演示证据交换的接受和拒绝规则。它们不连接真实的 OmicsWeave、TheraWeave 或 EphysWeave 安装。

| 文件 | 期望 |
|---|---|
| `theraweave_predicted_admet.json` | 可作为 predicted 证据附带，不进入浓度表 |
| `omics_flux.json` | 可作为通量证据附带 |
| `omics_mislabeled_concentration.json` | 拒绝。通量不能标成药物浓度 |
| `ephys_endpoint.json` | 可读出，含同步误差，不移动剂量时钟 |
| `ephys_missing_sync.json` | 拒绝。缺少同步信息时不做对齐 |
