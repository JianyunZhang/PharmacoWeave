# JSON Schema

这些文件由 Pydantic 模型导出。修改模型后运行：

```bash
python -m pharmacoweave.schema_export
```

测试会比较导出结果和这里的文件。直接编辑 JSON 会在 CI 中失败。

| 文件 | 模型 |
|---|---|
| `result_bundle.schema.json` | 结果包清单 |
| `evidence_exchange.schema.json` | 跨项目证据交换 |
| `pharmacology_study.schema.json` | 内存中的研究 |
| `analysis_records.schema.json` | 模型、拟合、NCA 和诊断记录 |
| `nca_config.schema.json` | 第一版 NCA 配置 |
