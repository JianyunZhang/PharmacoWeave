# 路线

阶段按科学验收切片，不按日历。

## 已完成：Phase 0 骨架（0.0.1）

- 可安装的 `pharmacoweave` 包、CI、MIT 许可和引用文件。
- 研究实体、封闭单位、关系校验、BLQ 策略。
- 合成静脉推注例子和严格校验。
- 结果包、HTML 校验清单、JSON Schema。
- 跨项目证据守卫和 NCA 配置拒绝超范围设计。
- 一室静脉推注解析解。

## 已完成：Phase 1 第一种 NCA 和一室静脉拟合（0.1.0）

`run_nca` 只计算 `configs/nca_single_dose_iv.yaml` 锁定的设计：单次静脉推注、连续采血、linear-up/log-down。口服、输注、多次给药和稀疏采样仍然拒绝。

- 梯形规则、C0 回推和终末相选点有独立参考转录，并与冻结数字比较。容差见 `docs/validation.md`。
- 终末相少于 3 点或调整 R² 低于 0.90 时，省略 λz、t½ 和 AUC∞。外推比例高于 0.20 时省略 AUC∞、CL 和 Vz。
- BLQ 保持删失。合成研究的 24 h 点不是 0，也不是定量点。
- 一室静脉推注拟合回收 CL 与 V。不收敛、贴边界或不可辨识时诊断为失败，估计缺席。
- 分析路径写出浓度–时间图。分析包的 `analysis_executed` 为 true。只做校验的包仍为 false。
- `schema_version` 保持 0.1.0。JSON Schema 没有新增字段。

仍延后：PD、DMPK、毒理、群体 PK、PBPK、QSP，以及这些阶段的笔记本。

## Phase 2：PD、DMPK、毒理和 TheraWeave 往返

在 Phase 1 的暴露结果之上，增加效应观测的时间尺度和量纲、典型体外 ADME 记录，以及带 `derived_from` 的 `PharmacologyEvidence`。预测与实测继续分开。`pd_models`、`pkpd`、`dmpk`、`toxicology` 从 reserved 改为有测试的实现之后才能对外声称可用。

## Phase 3：组学与电生理终点

在身份链接和同步误差已经存在的前提下，做一个可追溯例子：`dose_event → exposure → PD endpoint → source artifact`。对齐必须显式使用 `sync_error_ms` 和 `timebase_id`。通量仍然不是浓度。

## Phase 4：外部 PopPK、PBPK 和 QSP

通过进程外适配器接入经过选择的引擎。验收是适配器测试、数值交叉验证、版本记录和写明的 context of use。本项目不实现第二套 PK-Sim 或 NONMEM。

## 持续约束

- 新算法先有基准和失败用例，再出现在公开函数里。
- 负结果、删失、排除原因和局限性留在输出里。
- 人类临床数据保持研究用标识，例子里不放可识别的受试者信息。
