# PharmacoWeave

定量药理学研究用 Python 工具包。当前版本 **0.1.0** 在 Phase 0 的数据契约之上，计算锁定设计的 NCA，并拟合一室静脉推注的 CL 与 V。它不给出给药建议。

**Research use only.** 本软件不是医疗器械，不构成监管申报、GxP 或临床决策工具。参考 ICH M15 的记录结构，并不表示本软件满足该指南。

## 安装

```bash
pip install -e ".[dev]"
pharmacoweave validate examples/data/iv_bolus_single_dose --strict
```

需要 Python 3.11 或更新版本。

## 这一版能做什么

校验研究并写出只含数据清单的结果包。`analysis_executed` 保持 false：

```python
from pharmacoweave import PharmacologyStudy

study = PharmacologyStudy.from_files(
    "examples/data/iv_bolus_single_dose",
    mode="strict",
)
study.export_bundle("results/iv-bolus-validation", run_id="iv-bolus-validation")
```

对同一份合成研究运行已锁定的分析。NCA 只接受单次静脉推注、连续采血和 linear-up/log-down。拟合只估计一室静脉推注的 CL 与 V。成功写出分析包时 `analysis_executed` 为 true：

```bash
pharmacoweave analyze examples/data/iv_bolus_single_dose \
  --config configs/nca_single_dose_iv.yaml \
  --bundle results/iv-bolus-analysis \
  --strict
```

分析目录包含 NCA 表、参数估计、预测与残差、浓度–时间图，以及研究用 HTML 报告。终末相不合格、外推过大、拟合不收敛或不辨识时，相应产物列在 `omitted` 里，不会被标成通过。

```bash
pharmacoweave capabilities
```

`core`、`ingestion`、`bioanalysis`、`pk_nca`、`pk_models`、`diagnostics`、`reporting` 可以执行，范围限于上面的锁定设计。`integrations` 只有文件级证据守卫。`pd_models` 及之后的包可以导入，导入不会运行分析。

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

`configs/nca_single_dose_iv.yaml` 锁定第一种 NCA：单次静脉推注、连续采血、linear-up/log-down。口服、输注、多次给药和稀疏采样会被拒绝。数值容差写在 `docs/validation.md`。
