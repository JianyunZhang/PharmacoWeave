# 例子

`data/iv_bolus_single_dose/` 是可校验的合成 PK 研究。`cross_project_demo/` 是证据交换夹具，包含应当接受和应当拒绝的文件。

校验和分析都使用命令行。笔记本仍未加入。

```bash
pharmacoweave validate examples/data/iv_bolus_single_dose --strict
pharmacoweave analyze examples/data/iv_bolus_single_dose \
  --config configs/nca_single_dose_iv.yaml \
  --bundle results/iv-bolus-analysis \
  --strict
```
