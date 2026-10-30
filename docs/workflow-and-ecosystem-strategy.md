# PharmacoWeave 流程优化、SDGR 对比与生态整合战略方案

> **文档版本**：v1.0  
> **编写日期**：2026-10-30  
> **适用版本**：PharmacoWeave 0.0.1 (Pre-Alpha, commit `02240d9`) 及后续演进规划  
> **关联文档**：`docs/project-evaluation.md`、`docs/architecture.md`、`docs/development-plan.md`、`docs/roadmap.md`  
> **文档定位**：针对运行时管线与研发步骤的修订方案、对标 SDGR（Schrödinger）软件套件的异同剖析、核心差异化创新路径及三层生态整合蓝图  
> **与当前 main 的关系**：文中的“现状”按 **0.0.1** 书写。**0.0.2**（commit `7fd9028`）已经落地单次静脉推注 NCA 和对数线性 CL/V 拟合，绘图使用标准库 SVG，没有把 SciPy 或 Matplotlib 放进必装依赖。派生分析数据集、方法学准入门禁、加权非线性拟合、残差图和证据包导出仍是后续建议。当前进度以 `docs/roadmap.md` 为准。

---

## 一、 执行流程与研发步骤的修改及改进方案

结合当前代码实现（`src/pharmacoweave/`）与既有路线规划（`docs/roadmap.md`、`docs/project-evaluation.md`），需从 **运行时数据与分析执行链路（Runtime Pipeline）** 和 **工程研发推进顺序（Development Steps）** 两个层面对现有流程进行修正。

### 1.1 运行时执行链路（Runtime Pipeline）的四项结构性修正

#### 修正 1：解耦只读输入 `PharmacologyStudy` 与单次运行产物 `StudyRunResult`
- **现状痛点**：
  - `Record` 基类配置为 `frozen=True`（不可变模型），`study.attach_candidates()` 返回新实例而非原地修改，这与 `docs/development-plan.md` 第 7.5 节预想的原地调用链存在语义冲突。
  - 更关键的是，进入 Phase 1 后，NCA 计算出的 `NcaParameterSet` 与模型拟合产出的 `ModelFit`、`ParameterEstimate`、`DiagnosticRecord`（汇总为 `AnalysisRecords`）若直接挂载回 `PharmacologyStudy` 实例，会污染原始研究输入对象，违背“同一份干净实验数据可使用不同配置/模型重复运行并输出独立 `run_id` 结果包”的设计原则。
- **改进方案**：
  引入独立的单次运行结果聚合对象 `StudyRunResult`（或 `AnalysisRun`），将运行流重构为纯函数式三段管道：

```mermaid
flowchart LR
    A["PharmacologyStudy<br/>(只读实验输入)"] --> C["Execution Engine<br/>(validate / run_nca / fit_model)"]
    B["AnalysisConfig / ModelSpec<br/>(分析配置与模型声明)"] --> C
    C --> D["StudyRunResult<br/>(AnalysisRecords + 派生表 + 图表 + 审计日志)"]
    D --> E["export_bundle()<br/>(写出带 SHA-256 的不可变目录)"]
```

- **推荐目标 API 形态**：
  ```python
  from pharmacoweave import PharmacologyStudy
  from pharmacoweave.pk_nca import load_nca_config, execute_nca

  study = PharmacologyStudy.from_files("data/experiment_001/", mode="research")
  study = study.attach_candidates("theraweave_candidates/admet.json")

  nca_cfg = load_nca_config("configs/nca_single_dose_iv.yaml")
  run_result = execute_nca(study, config=nca_cfg, run_id="study_001_nca_run_01")
  manifest = run_result.export_bundle("results/study_001_nca_run_01/")
  ```

#### 修正 2：补齐悬空的“分析数据集派生层（Derived Analysis Dataset / ADPC 层）”
- **现状痛点**：
  - `bioanalysis/qc.py` 中的 `apply_blq_policy()` 生成了不可变的 `AnalysisObservation`，但目前在整个运行链路中处于“悬空”状态：既未接入 `pk_nca/runner.py`，在 `core/bundle.py` 的 `_PATHS` 清单中也未定义落盘路径。
  - 缺少数值计算前的**单位归一化**（如同一研究中混有 `ng/mL` 与 `ug/mL` 观测记录）以及**技术平行样（`replicate_id`）聚合策略**（算术均值、几何均值或保留独立测量）。
- **改进方案**：
  在跨表校验通过之后、进入 NCA 或非线性拟合之前，显式插入 `derive_analysis_dataset()` 标准步骤，并将派生结果写入 `analysis/derived_observations.parquet`：
  1. **单位标准化**：按目标分析量纲统一转换时间（如 `h`）、浓度（如 `ug/mL` 或 `nmol/L`）与剂量单位；
  2. **BLQ 与超限规则派生**：调用 `apply_blq_policy()`，保留原始 `Observation` 不变，生成带 `rule` 与 `caution` 审计字段的派生列；
  3. **平行样处理**：显式按配置聚合同一 `sample_id` 下的 `replicate_id` 记录并记录聚合前离散度（CV%）；
  4. **异常与排除过滤**：对标记了 `exclusion_reason` 的点设置 `included=False`，确保不进入积分/拟合但完整保留在派生表与图表中（用特殊符号标示被剔除点）。

#### 修正 3：将重复校验重构为“三级递进校验门禁（3-Stage Validation Gate）”
- **现状痛点**：
  - 当前 `from_files()`、`run_nca()`、`export_bundle()` 会三次重复触发相同的 `validate_study()`，但缺少**实验设计与所选分析方法匹配度**的校验（例如将口服或多剂量数据集传入仅支持 `iv_bolus + single_dose` 的 `NcaConfig` 时，现有 `runner.py` 未检查 `study.dose_events` 的真实给药途径和给药频次）。
- **改进方案**：
  拆分为职责正交的三级门禁：
  - **L1 记录模式门禁（Schema & Record Gate）**：单行字段类型、时区、封闭单位词表、BLQ 禁填数值校验（文件读取时执行一次）；
  - **L2 跨表关系与证据守卫门禁（Relational & Evidence Gate）**：外键完整性、绝对/相对时钟一致性、Assay 基质与 LLOQ 核对、跨项目证据防篡改守卫（产出可复用的 `ValidationReport` 摘要）；
  - **L3 方法学适用性门禁（Method Eligibility Gate）**：在运行具体 NCA 或房室模型前，逐个 `Subject` 核查给药途径（IV vs Extravascular）、给药事件数（单次 vs 多次）、基质类型（血浆 vs 尿液/组织）与有效定量点数是否满足算法前置假设；不满足的受试者不抛未捕获异常中断整批任务，而是写入 `NcaParameterSet.excluded_reason` 并附明确诊断码。

#### 修正 4：前置增加“声明式原始实验表映射器（Declarative Ingest Mapper）”
- **现状痛点**：
  - `ingestion/tabular.py` 要求输入 CSV/Parquet 列名严格对齐内部标准字段，而实验室 LC-MS/MS 工作站（Sciex Analyst、MassLynx、Thermo Xcalibur）或 CRO 交付的原始宽表/异构长表列名差异极大，人工改表既易出错又会破坏原始文件的哈希溯源链。
- **改进方案**：
  支持可选的 `ingest_mapping.yaml` 映射配置：先对未经手工篡改的实验室原始导出文件计算 `source_files.json` SHA-256 摘要，再通过声明式列映射、单位绑定与宽表转长表（melt）规则在内存中构建 `PharmacologyStudy`。

---

### 1.2 研发推进步骤（Development Roadmap）的四项顺序调优

对比 `docs/project-evaluation.md` 原定的 Phase 1 顺序，建议做如下工程顺序调整：

```mermaid
flowchart TD
    S0["Step 1.0: 冻结外部黄金基准与病态失败集<br/>(Benchmark-First TDD)"] --> S1["Step 1.1: 可选依赖分层配置<br/>(保持 core 轻量，拆分 pk/plotting extras)"]
    S1 --> S2["Step 1.2: 分析数据集派生层与 L3 准入门禁<br/>(单位归一 / BLQ 派生落盘 / 平行样聚合)"]
    S2 --> S3["Step 1.3: NCA 核心计算引擎实现<br/>(Linear-up Log-down + λz 自动选点)"]
    S3 --> S4["Step 1.4: 一室静脉加权非线性拟合器<br/>(ParameterEstimate + 收敛与局限性诊断)"]
    S4 --> S5["Step 1.5: 诊断图表生成与结果包解封<br/>(线性/半对数图 + 残差图 + analysis_executed=True)"]
    S5 --> S6["Step 1.6: 提前打通 PharmacologyEvidence 导出<br/>(形成 PharmacoWeave → TheraWeave 最小闭环)"]
```

1. **Benchmark-First 顺序前置（将原 Step 1.6 提至 Step 1.0）**：
   严格落实 `CONTRIBUTING.md` 规范，在写任何一行 NCA 或拟合实现代码前，先在 `tests/numerical_benchmarks/` 和 `tests/regression/` 中固化由外部独立参考工具（如 R `PKNCA`）生成的标准真值 JSON，同时纳入**病态与拒算测试集**（浓度反弹、双峰曲线、仅 2 个定量点、终末外推 > 20%、口服数据误用 IV 配置）。
2. **核心包与重型科学计算包分层装配（修正依赖策略）**：
   避免把 `scipy`、`matplotlib` 强塞进基础必装依赖。基础安装 `pip install pharmacoweave` 保持极简，供 CI 与兄弟项目快速执行 Schema 及证据守卫校验；NCA 梯形法与对数线性回归仅依赖标准库 `math` 或轻量 `numpy`；非线性优化与静态绘图通过 `pyproject.toml` 的 `[project.optional-dependencies]`（如 `pharmacoweave[pk]`、`pharmacoweave[plotting]`、`pharmacoweave[standard]`）按需安装。
3. **提前在 Phase 1 末打通 `PharmacologyEvidence` 导出闭环**：
   不等到 Phase 2，只要单室 IV NCA/拟合跑通，即提供 `export_evidence_exchange()` 将合格的体内暴露参数（`CL`、`V`、`t½`、`AUC∞`，标记 `evidence_type="estimated"` 及 `derived_from` Bundle 哈希）导出为标准证据包，尽早验证与 `TheraWeave` 的往返联调。
4. **增强 CLI 全生命周期体验**：
   在 `cli.py` 现有 `validate`、`capabilities`、`export-schemas` 基础上，新增 `pharmacoweave init-study`（生成带注释的标准目录脚手架）与 `pharmacoweave verify-bundle`（离线核对已有结果包内全部文件的 SHA-256 完整性及 Schema 合规性）。

---

## 二、 与现有 SDGR（Schrödinger 薛定谔）系列软件套件的异同点剖析

> **对标基准说明**：以下以计算药物发现行业广泛使用的 **SDGR（Schrödinger, Inc.）软件套件**（涵盖 Maestro、Glide、Desmond、FEP+、QikProp、AutoQSAR/DeepAutoQSAR、Predictive Toxicology 及企业协同平台 LiveDesign）作为参照系展开对比。

### 2.1 相同点（业务交集与共同目标）

1. **共同服务于新药研发 DMTA（Design-Make-Test-Analyze）闭环与转化降险**：两者均致力于在临床前阶段尽早识别药代动力学缺陷与潜在毒性风险，提升候选化合物（PCC）成药率。
2. **高度关注 ADMET / DMPK 核心性质表征**：SDGR 通过 QikProp、机器学习模型及脱靶物理模拟评估溶解度、膜通透性、血脑屏障穿透、血浆蛋白结合与 CYP/hERG 风险；PharmacoWeave（协同 TheraWeave）同样围绕清除率、分布容积、暴露量、体外 ADME 与毒理信号建立定量描述。
3. **强调干湿实验协同与跨学科数据汇聚**：SDGR 的 **LiveDesign** 致力于打通计算预测（*In Silico*）与实验室生物测定（*In Vitro / In Vivo* Assay）的数据孤岛；PharmacoWeave 亦以连接上游计算候选预测与下游实测生物分析数据为核心场景。
4. **面向多治疗模态拓展**：SDGR 已从小分子拓展至多肽与生物大分子（BioLuminate）；PharmacoWeave 在架构规划（`docs/development-plan.md` 第 5.4 节）中也前瞻性地为小分子、多肽、mRNA/siRNA/ASO 及基因编辑预留了差异化动力学建模边界。

### 2.2 核心不同点（多维对比矩阵）

| 对比维度 | SDGR 系列软件套件（Maestro / FEP+ / QikProp / LiveDesign） | PharmacoWeave（及所属 Weave 开源软件家族） |
| :--- | :--- | :--- |
| **底层科学范式** | **自下而上（Bottom-Up）微观物理第一性原理与结构模拟**：以量子化学（QM/MM）、分子对接、分子动力学和自由能微扰（FEP+）为核心，从三维分子与靶点结构出发做正向先验预测。 | **自上而下与中间层（Top-Down / Middle-Out）实验定量药理学**：以活体/体外实验测定数据为核心，围绕生物分析质控、NCA、PK/PD 动态反演与跨尺度暴露–效应关系做定量推断。 |
| **DMTA 主战场** | 强于 **Design（分子设计与虚拟筛选）**——解决“合成哪个分子亲和力更高、预测类药性更好”。 | 强于 **Test → Analyze（实验解析与转化证据确证）**——解决“化合物进入生物体后真实暴露如何、是否覆盖靶点有效浓度、药效与毒性时间窗如何定量”。 |
| **对“预测 vs 实测”的认识论策略** | **宽表融合与效率优先**：在 LiveDesign 中倾向于将计算预测列与实验测定列并置，用于多参数优化打分（MPO），但在导出至下游分析脚本时容易模糊证据等级与适用前提。 | **认识论强隔离与类型守卫（Epistemic Guardrails）**：在 Schema 与代码层硬编码区分 `observed / estimated / predicted / literature / hypothesis`，严禁预测属性冒充实测暴露，严禁混淆 `CL` 与 `CL/F`。 |
| **时空尺度与生物学纵深** | 聚焦 **飞秒–微秒级原子热运动** 与微观结合口袋构象；不直接处理在体多组学代谢重编程或毫秒级神经电生理功能表型。 | 聚焦 **毫秒级（电生理）至小时/天/周级（PK/PD、毒理与组学响应）** 的活体系统动力学，原生支持与多组学（`OmicsWeave`）及神经功能终点（`EphysWeave`）跨尺度对齐。 |
| **软件形态、成本与可复现工程哲学** | 商业闭源重型套件、高昂商业授权壁垒、以专有 GUI / 企业云数据库为中心，外部审稿人或合作方难以零成本重现完整计算流水线。 | MIT 开源、Python 原生、Code-First & File-First，产出带 SHA-256 校验链的不可变 `PharmacologyResultBundle`，天然适配 Git 版本控制、CI/CD 与开放科学同行复现。 |

---

## 三、 PharmacoWeave 应当体现的六大差异性与创新性方向

PharmacoWeave **不应在三维分子模拟/QSAR 赛道模仿 SDGR，也不应在重型闭源 GxP 图形界面赛道追赶传统商业 PK 软件（如 Phoenix WinNonlin）**，而应聚焦以下六个高辨识度创新支点：

### 1. “认识论感知（Epistemically-Aware）”的防幻觉药理类型系统
在自动化工作流与 AI4Science 快速普及的背景下，药理分析最隐蔽的风险是 **“证据洗白（Evidence Laundering）”**——例如把 AI/QikProp 预测的清除率当作动物实测清除率混入训练表、把 BLQ 静默替换为 0 扭曲终末半衰期、或者在拟合未收敛/未声明局限性时直接给出合格结论。
PharmacoWeave 将 **封闭量纲词表、给药途径语义绑定（`CL` vs `CL/F`）、可审计非破坏性 BLQ 派生、终末相准入判据、跨项目证据守卫以及强制填写的 `context_of_use` / `limitations`** 直接编译进 Pydantic 类型系统与 JSON Schema，打造 **AI4Pharmacology 时代的确定性科学安全护栏**。

### 2. 跨时间尺度“体内药物暴露 × 神经电生理 × 动态代谢通量”联合建模锚点
这是商业 CADD 套件空白、传统 PK 软件也极少原生支持的交叉学科高地：
- **神经精神药理学（PK–Neural-PD）**：依托内置的 `timebase_id`、`sync_error_ms` 时钟同步契约和受试者 `identity_map`，严谨连接血浆/脑脊液/脑组织药物暴露曲线与 `EphysWeave` 提取的神经元发放率（Spike Rate）、局部场电位（LFP）频段功率及行为学终点，精准刻画血脑屏障穿透滞后与效应室动力学。
- **暴露驱动的代谢流药理学（PK–MFA）**：连接个体真实药物暴露与 `OmicsWeave` 的同位素代谢通量分析（MFA）及通路重编程证据，严格区分药物浓度与内源性代谢通量。

### 3. “先验计算预测 ↔ 后验体内实测”定量校准与偏差审计闭环（In Silico-to-In Vivo Audit）
传统流程中，前端计算预测（如 SDGR QikProp 或 TheraWeave 预测 ADMET）在分子送去合成并测出动物 PK 后往往被束之高阁，缺乏系统性的误差复盘。
PharmacoWeave 可基于已挂载的 `evidence_type="predicted"` 记录与实测计算的 `NcaParameterSet` / DMPK 数据，原生提供 **预测偏差定量审计器**：自动计算预测折叠误差（Fold Error）、平均绝对折叠误差（AAFE）、系统偏差（AFE）及结构系列相关性，自动标记“预测高清除但实测低清除”等异常化合物，反向指导上游计算模型重校准。

### 4. 面向新兴治疗模态的“多基质组织暴露–持久效应”原语
针对 siRNA/ASO 核酸药物、mRNA-LNP、多肽偶联物与基因编辑疗法，血浆游离浓度往往在数小时内降至 BLQ，而靶器官（肝、脑、肾等）内的活性成分与下游药效可持续数周至数月。
充分发挥 PharmacoWeave 原生区分 15+ 种生物基质（`plasma`、`csf`、`brain`、`liver` 等）与多重分析物角色（原药、递送载体、活性代谢物、靶标标志物）的优势，沉淀一套轻量开源的 **“血浆–靶组织–持久 PD 解耦建模标准模板”**。

### 5. 践行 ICH M15（MIDD）理念的“Bundle-as-Evidence”可计算证据包与差异比对
不同于传统软件输出静态 PDF 或散装 Excel，`PharmacologyResultBundle` 不仅完整封装输入文件哈希、运行环境、结构化 Parquet 表，还显式记录 **`omitted`（未执行的分析项及原因）** 与 **`limitations`（模型假设与适用边界）**。在此基础上开发 **`pharmacoweave diff <bundle_a> <bundle_b>`** 工具，可一键对比两次运行（如更换 BLQ 策略或更换误差模型前后）在入模样本、参数估计与结论风险上的精确差异。

### 6. 原生 Agent-Ready / MCP-Ready 定量药理工具内核
凭借完备的 JSON Schema 契约、结构化的 `ValidationIssue` 诊断码和拒绝静默猜值的报错机制，PharmacoWeave 天然适合封装为科研智能体（AI Agent / Model Context Protocol Server）的标准工具底座——让大模型负责自然语言交互与方案起草，由 PharmacoWeave 负责 100% 确定性的单位换算、实验数据质控、NCA/PKPD 求解与证据包签发。

---

## 四、 整合已有生态发挥软件套件优势的实施方案

建议采取 **“内环闭环筑壁垒、中环开源借算力、外环商业做互补”** 的三层同心圆生态整合策略：

```mermaid
flowchart TB
    subgraph Ring3 ["外环：商业 CADD 与实验室仪器数据生态（错位互补）"]
        SDGR["SDGR 套件导出<br/>(Maestro / QikProp / FEP+ / LiveDesign)"]
        LabRaw["实验室原始数据<br/>(LC-MS/MS 导出表 / CRO Excel / LIMS)"]
    end

    subgraph Ring1 ["内环：Weave 自有软件家族四剑客（核心证据飞轮）"]
        OW["OmicsWeave<br/>多组学差异 / 通路 / MFA 通量"]
        TW["TheraWeave<br/>靶点评估 / 候选分子 / 预测 ADMET"]
        PW["PharmacoWeave<br/>生物分析 QC / NCA / PKPD / 证据守卫"]
        EW["EphysWeave<br/>Plexon Spike / LFP / 行为表型"]

        OW -->|"疾病靶点与生物标志物证据"| TW
        TW -->|"候选化合物 + predicted ADMET"| PW
        PW -->|"estimated 实测暴露裕度 + 预测误差校准"| TW
        PW <-->|"连续暴露协变量 ↔ 分子/通量响应"| OW
        PW <-->|"给药时钟上下文 ↔ 神经电生理 PD 终点"| EW
    end

    subgraph Ring2 ["中环：开源数据科学、生信标准与成熟药理引擎（算力与底座复用）"]
        Lake["PyArrow + Parquet + DuckDB<br/>零部署跨百批次 SAR-PK 查询"]
        ExtEngines["PKNCA / nlmixr2 / OSP (PK-Sim & MoBi)<br/>进程外容器化重型引擎适配"]
        OpenStd["RDKit / NWB / AnnData / CDISC SEND 子集"]
    end

    SDGR -->|"CADD 预测属性桥接器 (标为 predicted)"| TW
    SDGR -->|"先验属性对照包"| PW
    LabRaw -->|"声明式 YAML 列映射器"| PW
    PW <-->|"零拷贝列式扫描"| Lake
    PW <-->|"标准文件交换子进程"| ExtEngines
    OpenStd -.->|"格式与元数据互操作"| Ring1
```

### 4.1 内环整合：构建 Weave 四剑客家族的统一交换协议与三大示范工作流

1. **制定轻量化《Weave 跨项目交换元规范（Weave Common Exchange Spec）》**：
   保持四个仓库代码零运行时耦合，但对齐以下四项元数据标准：
   - 统一跨项目实体 URI 命名空间：`weave://<project>/<study_id>/<entity_type>/<id>`；
   - 统一 `identity_map` 个体映射表结构与 `review_status`（`pending` / `accepted` / `rejected`）审核语义；
   - 统一 `manifest.json` 顶层溯源字段（`schema_name`、`schema_version`、`producer`、`files` SHA-256、`derived_from`）；
   - 统一 5 级 `evidence_type` 词表（`observed`、`estimated`、`predicted`、`literature`、`hypothesis`）。
2. **打造三条跨仓库联合端到端示范管线（Flagship Cross-Weave Pipelines）**：
   - **示范管线 1（分子设计—湿实验 PK 校准飞轮）**：`TheraWeave` 导出候选分子与先验 ADMET 预测包 $\rightarrow$ `PharmacoWeave` 导入动物体内 PK 数据完成 NCA 与一室拟合，计算有效暴露窗口并输出预测偏差审计包 $\rightarrow$ 回流 `TheraWeave` 更新候选分子优先级。
   - **示范管线 2（中枢神经药物 PK–Ephys–Behavior 联合分析）**：`PharmacoWeave` 管理给药时钟与血/脑暴露轨迹 + `EphysWeave` 提供带 `sync_error_ms` 的 LFP/Spike 特征序列 $\rightarrow$ 联合建立药物暴露–神经环路电活动–行为表型时间对齐分析包。
   - **示范管线 3（实测暴露驱动的多组学剂量响应解析）**：将 `PharmacoWeave` 产出的个体级实测暴露指标（`AUC_last`、`Cmax`）导出为协变量包传入 `OmicsWeave`，替代粗糙的名义剂量分组，精准剥离个体吸收/清除差异对转录组与代谢通量的混杂影响。

### 4.2 中环整合：无缝复用开源数据湖与专业建模引擎

1. **基于 Parquet + DuckDB 构建零运维“跨批次药理证据库”**：
   由于 `PharmacologyResultBundle` 已统一采用强 Schema 约束的 Parquet 表存储，无需引入重型关系型数据库服务，仅需提供一个轻量查询模块，即可利用 DuckDB 直接跨本地或对象存储上的数百个历史 `ResultBundle` 目录执行联邦 SQL 查询（例如一键汇总某靶点项目所有批次化合物的清除率、生物利用度与 Assay BLQ 比例）。
2. **与开源标准格式及专业建模引擎分层对接**：
   - **化学信息学可选联动（RDKit）**：在可选扩展中利用 RDKit 校验 `Compound` 的 SMILES/InChIKey 一致性与精确分子量，为质量浓度与摩尔浓度安全换算提供经过核验的分子量来源；
   - **领域标准互操作**：对接 **NWB（Neurodata Without Borders）** 的时间序列同步元数据（衔接 `EphysWeave`）、**AnnData / MuData** 的样本元数据表（衔接 `OmicsWeave`），并预留面向 **CDISC SEND（PC/EX/DM 域）** 的单向导出转换器；
   - **重型引擎进程外编排**：在 Phase 4 通过标准化子进程与容器协议封装 **nlmixr2**（群体 PK/PD NLME）与 **Open Systems Pharmacology（PK-Sim / MoBi）**（全身生理药代 PBPK），由 PharmacoWeave 负责前端数据质控守卫与后端不可变结果包归档。

### 4.3 外环整合：提供 SDGR 等主流商业套件的“错位互补桥接器”

针对已部署 SDGR（Schrödinger）套件的科研课题组或生物技术团队，不采取排他替代策略，而是在 `integrations/` 下规划通用的 **外部 CADD 属性导入适配器（如 `cadd_importers/schrodinger.py`）**：
- **自动解析与证据降级保护**：读取 Maestro / QikProp / LiveDesign 导出的 SDF 或 CSV 表格，将 `QPlogPo/w`、`QPlogS`、`QPPCaco`、`QPlogBB`、`QPlogKhsa` 及 FEP+ 预测结合自由能/亲和力自动映射为标准 `EvidenceRecord`，强制绑定 `evidence_type="predicted"`、`qc_status="not_applicable"` 及工具版本溯源，严防直接混入实测表。
- **形成“上半场 SDGR 物理预测 + 下半场 Weave 实验确证”的黄金组合**：让计算化学团队继续使用熟悉的 SDGR 完成虚拟筛选，随后无缝将候选分子移交至 **PharmacoWeave + OmicsWeave + EphysWeave** 开展体内外 DMPK 质控、PK/PD 建模、组学机制解析与预测准确度闭环复盘，最大化整套软件生态的科研落地价值。

---

## 五、 近期优先落地行动项汇总

| 优先级 | 行动项 | 涉及模块 / 路径 | 预期产出与验收标准 |
| :---: | :--- | :--- | :--- |
| **P0** | 冻结外部 NCA 与一室拟合黄金基准集 + 病态失败集 | `tests/numerical_benchmarks/`、`tests/regression/` | 包含标准真值与至少 5 类拒算/警告边界用例 |
| **P0** | 重构 `StudyRunResult` 并补齐 `derive_analysis_dataset()` | `core/study.py`、`core/bundle.py`、`bioanalysis/qc.py` | 原始 `Study` 保持只读，BLQ 与单位归一派生表落盘为 Parquet |
| **P0** | 拆分可选依赖并实现 NCA + 一室 IV 非线性拟合与诊断图 | `pyproject.toml`、`pk_nca/`、`pk_models/`、`reporting/` | 核心包保持轻量；`analysis_executed=True` 端到端产出完整报告包 |
| **P1** | 补充声明式实验室原始表列映射器与 `PharmacologyEvidence` 导出 | `ingestion/`、`integrations/` | 支持原始 LC-MS/MS 导出表直接入库溯源，跑通与 `TheraWeave` 证据往返 |
| **P1** | 构建“先验预测–后验实测”偏差审计模块及外部 CADD 导入桥接 | `diagnostics/`、`integrations/` | 支持挂载 `TheraWeave` 或 SDGR 预测包，自动输出 Fold Error / AAFE 校准报告 |
