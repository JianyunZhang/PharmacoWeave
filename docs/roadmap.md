# 路线

阶段按科学验收切片，不按日历。

## 已完成：Phase 0 骨架（0.0.1）

- 可安装的 `pharmacoweave` 包、CI、MIT 许可和引用文件。
- 研究实体、封闭单位、关系校验、BLQ 策略。
- 合成静脉推注例子和严格校验。
- 结果包、HTML 校验清单、JSON Schema。
- 跨项目证据守卫和 NCA 配置拒绝超范围设计。
- 一室静脉推注解析解。

## 已完成：第一种 NCA 和一室静脉拟合（0.0.2）

- 单次静脉推注、连续采血、linear-up/log-down。
- 终末相规则 `pharmacoweave.lambda_z.adjusted_r2.v1`。点数不足或 adjusted R² 不够时不报告 λz、t½、AUC∞。
- 外推比例超过配置上限时保留数值，并由 `auc_inf_is_reportable` 返回 false。
- 对数线性拟合回收 CL 与 V。斜率不下降时诊断为失败，不输出清除率。
- 结果包在传入 NCA 结果后写出参数表和浓度–时间图，`analysis_executed` 为 true。

## 已完成：残差图、加权拟合和冻结 NCA 参考（0.0.3）

- 一室静脉推注拟合写出观测值–预测值和残差–时间 SVG。BLQ 不进入图，也不画成 0。
- 加权非线性最小二乘（权重 `1/y` 和 `1/y^2`）在带噪声合成曲线上回收 CL 与 V。浓度不下降或参数不可辨识时诊断为失败，不输出清除率。对数线性拟合仍然可用。
- `tests/regression/iv_bolus_quiet_nca.json` 是无噪声静脉推注 NCA 的第二份冻结参考。相对误差上限写在 `docs/validation.md`。这份数字用十进制独立算出，不把 PKNCA 或 Phoenix 引进依赖。

## 之后的 Phase 1 余量

稀疏采样、稳态、尿药、输注、血管外给药和多次给药仍拒绝计算。二室模型不在拟合库里。

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
