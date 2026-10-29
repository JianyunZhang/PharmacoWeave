# 配置

`nca_single_dose_iv.yaml` 是当前唯一接受的 NCA 设计：单次静脉推注、连续采血、linear-up/log-down。终末相至少 3 点，调整 R² 不低于 0.90，AUC∞ 外推比例的报告上限是 0.20。BLQ 规则是 `keep_as_censored`。

`PharmacologyStudy.run_nca` 和 `pharmacoweave analyze` 会按该文件计算。其他途径、多次给药和稀疏采样在计算前被拒绝。
