# PharmacoWeave

定量药理学研究用 Python 工具包。当前版本 **0.0.5** 提供单次静脉推注 NCA、带渐近标准误与 95% 置信区间的一室 CL/V 加权非线性拟合、原始效应尺度的直接与 Sigmoid Emax 暴露–效应拟合、可审计 BLQ 派生与平行样聚合、体外 ADME 与研究用毒理记录落盘，以及自动经证据守卫核验的 `PharmacologyEvidence` 与原生诊断 SVG。预测与实测严格分层。它不计算超出已声明契约的试验设计，不做体内外推，也不给出临床给药建议。

**Research use only.** 本软件不是医疗器械，不构成监管申报、GxP 或临床决策工具。参考 ICH M10 / M15 的记录结构，并不表示本软件满足该指南。

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

`results/` 会得到带 SHA-256 校验和的 Pharmacology Result Bundle：Parquet 表、来源、QC 和 HTML 校验报告。未运行分析时 `analysis_executed` 为 `false`，NCA 与拟合产物列在 `omitted` 中。

```bash
pharmacoweave capabilities
```

`core`、`ingestion`、`bioanalysis` 和 `pk_nca` 可以执行。`pk_nca` 只接受 `configs/nca_single_dose_iv.yaml` 里的设计。`pk_models` 提供一室静脉推注的对数线性拟合，以及带渐近标准误和 95% 置信区间的 `1/y`、`1/y^2` 加权非线性最小二乘。`pd_models` 在原始效应尺度估计固定基线的直接效应 Emax 和可选 Sigmoid Hill 系数。`pkpd` 只把直接效应连到显式浓度。`dmpk` 校验并持久化通透性、微粒体和肝细胞清除率记录。`toxicology` 记录研究用终点。`population_pk`、`pbpk`、`qsp` 和 `simulation` 仍是 `reserved`，导入不会运行分析。

纯 NCA 结果包（含线性与半对数浓度–时间图、派生观测表和自动生成的 NCA 证据；如需同时运行对数线性一室拟合可加 `--with-fit`）：

```bash
pharmacoweave nca examples/data/iv_bolus_single_dose \
  --config configs/nca_single_dose_iv.yaml \
  --bundle results/iv-bolus-nca \
  --run-id iv-bolus-nca
```

加权非线性一室静脉拟合（写入参数标准误与 95% CI、观测–预测图、残差图与加权相对残差图）：

```bash
pharmacoweave fit examples/data/iv_bolus_noisy \
  --method weighted_nls \
  --weight 1/y^2 \
  --bundle results/iv-bolus-fit \
  --run-id iv-bolus-fit
```

原始尺度直接效应或 Sigmoid Emax 拟合及结果包导出：

```bash
pharmacoweave emax examples/data/direct_effect_emax.json \
  --study examples/data/iv_bolus_single_dose \
  --bundle results/direct-emax \
  --run-id direct-emax
```

## 文档

| 文档 | 内容 |
|---|---|
| [可行性评估](docs/feasibility.md) | 对立项规划的评估和收窄决定 |
| [架构](docs/architecture.md) | 模块边界、数据流、技术选择 |
| [数据契约](docs/data-contracts.md) | 实体、表列、Schema 与结果包 |
| [验证](docs/validation.md) | 数值基准、冻结参考和验收靶值 |
| [互操作](docs/interoperability.md) | 与 OmicsWeave、TheraWeave、EphysWeave 的文件级边界 |
| [路线](docs/roadmap.md) | 分阶段交付与里程碑规划 |
| [项目评估报告](docs/project-evaluation.md) | 设计合理性、实现审计与阶段演进 |
| [工作流与生态策略](docs/workflow-and-ecosystem-strategy.md) | 研发执行流程与 SDGR / Weave 生态整合方案 |
| [立项规划原文](docs/development-plan.md) | 评估所依据的 v0.1 规划 |

## 合成示例

`examples/data/iv_bolus_single_dose/` 是无噪声的一室静脉推注数据：两只大鼠，CL = 0.20 L/h/kg，V = 1.0 L/kg。24 h 样品低于 0.01 µg/mL，以 BLQ 保存，数值列为空。这不是动物实验。

`examples/data/iv_bolus_noisy/` 是带乘性噪声的一室静脉推注。加权拟合在 10% 相对误差内回收 CL = 0.40 L/h 和 V = 2.0 L，并给出对数尺度渐近 95% 置信区间。24 h 仍然是 BLQ。

`examples/direct_effect_and_evidence.py` 与 `examples/noisy_emax_case.py` 分别是无噪声和带噪直接效应/Sigmoid Emax 基准，附带体外 ADME 和研究用毒理记录。这不是动物实验，也不是给药方案。
