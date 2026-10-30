# PharmacoWeave

定量药理学研究用 Python 工具包。当前版本 **0.0.3** 锁定实验数据模型、单位、质量控制、结果包和跨项目证据边界，并计算一种 NCA：单次静脉推注、连续采血。它可以用对数线性回归或加权非线性最小二乘估计该模型的 CL 和 V，并写出观测–预测图和残差图。它不计算其他试验设计，也不给出给药建议。

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

`core`、`ingestion`、`bioanalysis` 和 `pk_nca` 可以执行。`pk_nca` 只接受 `configs/nca_single_dose_iv.yaml` 里的设计。`pk_models` 提供一室静脉推注的对数线性拟合，以及权重为 `1/y` 或 `1/y^2` 的加权非线性最小二乘。`pd_models` 及之后的包可以导入，导入不会运行分析。

```bash
pharmacoweave nca examples/data/iv_bolus_single_dose \
  --config configs/nca_single_dose_iv.yaml \
  --bundle results/iv-bolus-nca \
  --run-id iv-bolus-nca
```

带 `--bundle` 时，这一命令同时写入对数线性拟合、`figures/goodness_of_fit.svg` 和 `figures/residual_diagnostics.svg`。BLQ 不画成 0。只跑加权拟合时：

```bash
pharmacoweave fit examples/data/iv_bolus_noisy \
  --method weighted_nls \
  --weight 1/y^2 \
  --bundle results/iv-bolus-fit \
  --run-id iv-bolus-fit
```

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

`examples/data/iv_bolus_noisy/` 是带乘性噪声的一室静脉推注。加权拟合应在 10% 相对误差内回收 CL = 0.40 L/h 和 V = 2.0 L。24 h 仍然是 BLQ。
