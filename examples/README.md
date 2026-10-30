# 例子

`data/iv_bolus_single_dose/` 是可校验的无噪声合成 PK 研究。`data/iv_bolus_noisy/` 是同一结构上的带噪声曲线，用来核对加权拟合。`cross_project_demo/` 是证据交换夹具，包含应当接受和应当拒绝的文件。

笔记本会在 NCA 或拟合可以端到端运行之后再添加。当前用命令行校验：

```bash
pharmacoweave validate examples/data/iv_bolus_single_dose --strict
```
