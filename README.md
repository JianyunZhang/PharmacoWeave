# PharmacoWeave

定量药理学研究用 Python 工具包。当前版本 **0.2.10**。它保留既有的 NCA、房室拟合、Emax、IVIVE、预测审计和结果包核验，并可对等剂量、等间隔的重复给药与已声明的稳态间隔做 NCA。尿药排泄按记录的尿液体积计算排泄量和肾清除，不把尿药浓度当作血浆浓度。一室恒速输注可以估计 `CL` 与 `V`。解析暴露模拟仍可带调用方给出的对数参数相关。调用方自备的可执行文件只有在同时给出版本和许可证文本时才会启动，并记录文件哈希。它不把未声明的稳态写成已经证明，也不给出临床给药建议。nlmixr2 和 PK-Sim 不被导入，也不作为内置后端启动。

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

`core`、`ingestion`、`bioanalysis` 和 `pk_nca` 可以执行。`pk_nca` 接受单次静脉推注（`configs/nca_single_dose_iv.yaml`）、单次静脉恒速输注（`configs/nca_single_dose_iv_infusion.yaml`）以及单次血管外给药（`configs/nca_single_dose_oral.yaml`，含 `oral`、`subcutaneous`、`intramuscular`）设计。`pk_models` 提供一室静脉推注的对数线性拟合，以及带渐近标准误和 95% 置信区间的一室静脉推注、一室恒速输注（`CL`、`V`）、一室血管外吸收（`ka`、`CL/F`、`V/F`）和二室静脉推注（`CL`、`V1`、`Q`、`V2`）加权非线性最小二乘。`pd_models` 在原始效应尺度估计固定基线的直接效应 Emax 和可选 Sigmoid Hill 系数。`pkpd` 只把直接效应连到显式浓度。`dmpk` 校验并持久化通透性、微粒体和肝细胞清除率记录。`toxicology` 记录研究用终点。`simulation` 可以在显式 `context_of_use` 下计算单次给药解析暴露区间，输出保持 `predicted`。`population_pk` 对已有个体估计做不加权和按标准误加权的几何均值，并读取调用方群体结果表。`pbpk` 校验调用方生理表，并计算肝加一个组织的代数稳态。`qsp` 计算一条线性周转方程的闭式解。这三项是 `partial`：进程内群体后验、全身动态 PBPK 和 QSP 网络仍由 `run_in_process` 拒绝。可执行的外部后端只有 `pharmacoweave-analytical`。

纯 NCA 结果包（含线性与半对数浓度–时间图、派生观测表和自动生成的 NCA 证据；静脉推注或血管外设计如需同时附加一室拟合可加 `--with-fit`）：

```bash
pharmacoweave nca examples/data/oral_single_dose \
  --config configs/nca_single_dose_oral.yaml \
  --bundle results/oral-nca \
  --run-id oral-nca
```

尿药排泄（单独配置；肾清除还需要调用方提供血浆 AUC。尿液体积缺失或 BLQ 时该区间排泄量为空）：

```bash
pharmacoweave urine-nca STUDY \
  --config configs/urine_excretion.yaml \
  --plasma-auc plasma_auc.yaml
```

调用方自备的进程外引擎（缺可执行文件、版本或许可证文本时不启动子进程；内置后端仍只有 `pharmacoweave-analytical`）：

```bash
pharmacoweave engine-slot REQUEST \
  --executable PATH \
  --engine-version 1.0 \
  --license-file LICENSE.txt
```

加权非线性房室拟合（支持 `one_compartment_iv_bolus`、`one_compartment_iv_infusion`、`one_compartment_extravascular`、`two_compartment_iv_bolus`，写入参数标准误与 95% CI、观测–预测图、残差图与加权相对残差图）：

```bash
pharmacoweave fit examples/data/two_compartment_iv_bolus \
  --model two_compartment_iv_bolus \
  --method weighted_nls \
  --weight 1/y^2 \
  --bundle results/two-comp-fit \
  --run-id two-comp-fit
```

原始尺度直接效应或 Sigmoid Emax 拟合及结果包导出：

```bash
pharmacoweave emax examples/data/direct_effect_emax.json \
  --study examples/data/iv_bolus_single_dose \
  --bundle results/direct-emax \
  --run-id direct-emax
```

显式生理假设探索性 Well-Stirred IVIVE（输出锁定 `predicted` 等级，无隐式物种默认表）：

```bash
pharmacoweave ivive examples/data/ivive_prediction_audit_demo.json \
  --study examples/data/iv_bolus_single_dose \
  --bundle results/ivive-bundle \
  --run-id ivive-bundle
```

先验预测与体内估算双向折叠误差校准审计（对比 TheraWeave CADD 先验与体内估算）：

```bash
pharmacoweave audit-predictions examples/data/iv_bolus_single_dose \
  --priors examples/cross_project_demo/theraweave_quantitative_priors.json \
  --nca-config configs/nca_single_dose_iv.yaml \
  --bundle results/audit-bundle \
  --run-id audit-bundle
```

声明式原始表格映射（支持 CSV/TSV/Parquet，自动计算并保留原始文件 SHA-256 摘要）：

```bash
pharmacoweave ingest \
  --mapping examples/data/ingest_mapping_demo/ingest_mapping.yaml \
  --output-dir results/mapped_study \
  --bundle results/mapped_bundle \
  --run-id mapped-study-01
```

结果包离线哈希与 Schema 验证（`verify-bundle`）：

```bash
pharmacoweave verify-bundle results/iv-bolus-validation
```

两次运行结果包差异对比（`diff`，对比元数据、文件树增删改及 Parquet 参数数值相对误差）：

```bash
pharmacoweave diff results/bundle_run1 results/bundle_run2
```

探索性单次给药暴露区间（预测等级，不是给药方案）：

```bash
pharmacoweave simulate examples/data/exposure_simulation_iv_bolus.json
pharmacoweave engine-run examples/data/exposure_simulation_iv_bolus.json \
  --engine pharmacoweave-analytical
```

## 文档

| 文档 | 内容 |
|---|---|
| [可行性评估](docs/feasibility.md) | 对立项规划的评估和收窄决定 |
| [架构](docs/architecture.md) | 模块边界、数据流、技术选择 |
| [数据契约](docs/data-contracts.md) | 实体、表列、Schema 与结果包 |
| [验证](docs/validation.md) | 数值基准、冻结参考和验收靶值 |
| [互操作](docs/interoperability.md) | 与 OmicsWeave、TheraWeave、EphysWeave 的文件级边界 |
| [路线](docs/roadmap.md) | 已完成范围，以及 0.2.6 之后的计划入口 |
| [群体、PBPK 与 QSP 计划](docs/phase4-population-pbpk-qsp.md) | 五项有限函数已在 0.2.6 落地。`run_in_process` 仍拒绝重型算法 |
| [0.2.6 之后的计划](docs/next-steps-after-0.2.6.md) | 加权两阶段、周转抑制、效应室、接线与三室已落地。重型算法仍拒绝 |
| [项目评估报告](docs/project-evaluation.md) | 立项阶段与 0.2.6 代码的对照，以及仍拒绝的功能 |
| [工作流与生态策略](docs/workflow-and-ecosystem-strategy.md) | 研发执行流程与 SDGR / Weave 生态整合方案 |
| [立项规划原文](docs/development-plan.md) | v0.1 规划；文首和第 11 节注明 0.2.6 完成状态与仍拒绝的功能 |

## 合成示例

`examples/data/iv_bolus_single_dose/` 是无噪声的一室静脉推注数据：两只大鼠，CL = 0.20 L/h/kg，V = 1.0 L/kg。24 h 样品低于 0.01 µg/mL，以 BLQ 保存，数值列为空。这不是动物实验。

`examples/data/iv_bolus_noisy/` 是带乘性噪声的一室静脉推注。加权拟合在 10% 相对误差内回收 CL = 0.40 L/h 和 V = 2.0 L，并给出对数尺度渐近 95% 置信区间。24 h 仍然是 BLQ。

`examples/data/oral_single_dose/`、`examples/data/iv_infusion_single_dose/` 与 `examples/data/two_compartment_iv_bolus/` 分别是无噪声单次口服给药（输出表观 `CL/F` 与 `Vz/F`，含显式分子量）、0.5 h 恒速静脉输注与二室静脉推注双指数合成基准。

`examples/direct_effect_and_evidence.py` 与 `examples/noisy_emax_case.py` 分别是无噪声和带噪直接效应/Sigmoid Emax 基准，附带体外 ADME 和研究用毒理记录。这不是动物实验，也不是给药方案。

`examples/ivive_and_prediction_audit.py` 演示完整的体外结合测定、显式生理假设探索性 Well-Stirred IVIVE 计算，以及先验预测（TheraWeave CADD 与 IVIVE）对后验体内实测 NCA 参数的双向折叠误差（Fold Error、AFE、AAFE）校准审计全流程。不是给药方案，也不推导人体起始剂量。
