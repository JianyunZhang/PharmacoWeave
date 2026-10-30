# Changelog

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
