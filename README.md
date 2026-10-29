# PharmacoWeave

定量药理学研究用 Python 工具包。当前版本 **0.0.1** 锁定实验数据模型、单位、质量控制、结果包和跨项目证据边界。它可以从已知参数计算一室静脉推注浓度，并生成数据校验报告。它不计算 NCA，不拟合 PK 模型，也不给出给药建议。

**Research use only.** 本软件不是医疗器械，不构成监管申报、GxP 或临床决策工具。参考 ICH M15 的记录结构，并不表示本软件满足该指南。

## 安装

```bash
pip install -e ".[dev]"
pharmacoweave validate examples/data/iv_bolus_single_dose --strict
```

需要 Python 3.11 或更新版本。

## 这一版能做什么

```python
from pharmacoweave import PharmacologyStudy

study = PharmacologyStudy.from_files(
    "examples/data/iv_bolus_single_dose",
    mode="strict",
)
study.export_bundle("results/iv-bolus-validation", run_id="iv-bolus-validation")
```

`results/` 会得到带校验和的 Pharmacology Result Bundle：Parquet 表、来源、QC 和 HTML 校验报告。清单里的 `analysis_executed` 为 false，NCA 与拟合产物列在 `omitted` 中。

```bash
pharmacoweave capabilities
```

`core`、`ingestion`、`bioanalysis` 可以执行。`pk_nca`、`pk_models`、`diagnostics`、`reporting`、`integrations` 只提供已测试的契约或参考方程。`pd_models` 及之后的包可以导入，导入不会运行分析。

## 文档

| 文档 | 内容 |
|---|---|
| [可行性评估](docs/feasibility.md) | 对立项规划的评估和收窄决定 |
| [架构](docs/architecture.md) | 模块边界、数据流、技术选择 |
| [数据契约](docs/data-contracts.md) | 实体、表列、结果包 |
| [验证](docs/validation.md) | 已有测试和第一版 NCA 的验收靶值 |
| [互操作](docs/interoperability.md) | 与 OmicsWeave、TheraWeave、EphysWeave 的文件级边界 |
| [路线](docs/roadmap.md) | 分阶段交付 |
| [立项规划原文](docs/development-plan.md) | 评估所依据的 v0.1 规划 |

## 合成示例

`examples/data/iv_bolus_single_dose/` 是无噪声的一室静脉推注数据：两只大鼠，CL = 0.20 L/h/kg，V = 1.0 L/kg。24 h 样品低于 0.01 µg/mL，以 BLQ 保存，数值列为空。这不是动物实验。

第一版准备接受的 NCA 设计写在 `configs/nca_single_dose_iv.yaml`：单次静脉推注、连续采血、linear-up/log-down。加载该文件不会计算 NCA。
