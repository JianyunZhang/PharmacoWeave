# Changelog

## 0.0.5 — 2026-10-30

- 直接效应与 Sigmoid Emax 默认改用原始效应尺度的阻尼 Levenberg–Marquardt 非线性最小二乘（`direct_emax_nls_v1` / `sigmoid_emax_nls_v1`），消除低浓度噪声放大以及带噪低浓度观测略低于固定基线 `E0` 时的拒算脆弱性；新增可选 Sigmoid Hill 系数估计（`fit_sigmoid_emax` / `estimate_hill`），同时保留 `method="linearized"` 供显式对比。高浓度平台覆盖不足、过渡区无观测、海森曲率比退化或效应方向矛盾时诊断为失败，不输出 `EC50` 或 `hill`。
- 基于收敛点 $J^T W J$ 协方差矩阵、Student-$t$ 分位数与正参数对数尺度 Delta 变换，为一室静脉推注加权非线性拟合（`CL`、`V`）与原始尺度 Emax/Sigmoid 拟合（`Emax`、`EC50`、`hill`）填充渐近标准误 `standard_error` 与 95% 置信区间（`ci_lower`、`ci_upper`、`ci_level=0.95`）。
- 新增纯标准库半对数浓度–时间图（`figures/concentration_time_semilog.svg`，含 LLOQ 阈值线、拟合曲线叠加与终末相选点标记）、加权相对残差图（`figures/weighted_residuals.svg`）及浓度–效应拟合曲线图（`figures/concentration_effect.svg`）；解除 `pharmacoweave nca --bundle` 对 `fit_iv_bolus` 的隐式耦合（默认纯 NCA，可通过 `--with-fit` 显式附加对数线性拟合）。
- 打通 Phase 2 全链路落盘与 Schema 导出：在结果包中持久化 `analysis/derived_observations.parquet`、`analysis/emax_parameters.parquet`、`analysis/emax_predictions.parquet`、`experimental/in_vitro_adme.parquet`、`experimental/tox_endpoints.parquet`，自动汇集并经证据守卫核验写入 `evidence/pharmacology_evidence.json`；导出 `pd_records.schema.json`、`dmpk_records.schema.json`、`toxicology_records.schema.json`；新增 `collapse_replicates` 可审计平行样聚合（算术/几何均值与 CV%）及 `pharmacoweave emax` CLI 子命令。

## 0.0.4 — 2026-10-30

- 效应观测记录时间尺度和量纲，并保持实测与预测分开。直接效应 Emax（基线固定、Hill 系数为 1）在无噪声合成曲线上回收 Emax 和 EC50。点数不足、效应平坦、效应等于基线，或方向与 Emax 符号相反时诊断为失败，不输出 EC50。AUC 和 Cmax 不能当作直接效应的浓度。BLQ 不替换成 0。滞后和间接反应不估计。
- 增加体外通透性、微粒体清除率和肝细胞清除率记录。血浆清除率单位会被拒绝。本版本不把这些结果换算成体内 CL。
- 增加研究用毒理终点。记录可以链接剂量或暴露，调用给药建议会失败。
- PharmacologyEvidence 用 `derived_from` 引用暴露或分析来源。预测与实测在 JSON 往返后等级不变。不导入 TheraWeave、OmicsWeave 或 EphysWeave。

## 0.0.3 — 2026-10-30

- 一室静脉推注拟合写出观测值–预测值和残差–时间 SVG。BLQ 不画成 0。`pharmacoweave nca --bundle` 和 `pharmacoweave fit --bundle` 把这两张图写入结果包的 `figures/`。
- 增加权重为 `1/y` 和 `1/y^2` 的一室静脉推注非线性最小二乘。带噪声合成例子在 10% 相对误差内回收 CL 和 V。浓度不下降或不可辨识时诊断为失败，不输出清除率。对数线性拟合仍然可用。
- 在 `tests/regression/iv_bolus_quiet_nca.json` 冻结无噪声静脉推注 NCA 的十进制参考值。相对误差上限写在 `docs/validation.md`。

## 0.0.2 — 2026-10-30

- 计算单次静脉推注、连续采血的 NCA：linear-up/log-down、对数线性回推 C0、adjusted R² 终末相。
- BLQ 和高于 ULOQ 的点不参与梯形积分，也不替换成数字。
- 外推比例超过配置上限时仍保存 AUC∞，并由 `auc_inf_is_reportable` 标成不宜报告。
- 用对数线性最小二乘拟合一室静脉推注，回收 CL 和 V。浓度不下降时不给出清除率。
- 结果包可以写入 `analysis/nca_parameters.parquet` 和浓度–时间图。`pharmacoweave nca` 运行该计算。

## 0.0.1 — 2026-10-29

研究用数据契约的第一版骨架。

- 增加研究、受试者、给药、样本、观测、检测和身份链接模型。
- 增加封闭单位词表、关系校验、可审计的 BLQ 策略和结果包。
- 增加无噪声一室静脉推注合成研究，以及一室静脉推注解析解。
- 增加 NCA 配置和结果 Schema。计算器尚未实现。
- 增加 OmicsWeave、TheraWeave、EphysWeave 的文件级证据守卫。
- 后期分析包可以导入，并标记为 reserved。
