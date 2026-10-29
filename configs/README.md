# 配置

`nca_single_dose_iv.yaml` 是 0.0.1 唯一接受的 NCA 设计：单次静脉推注、连续采血、linear-up/log-down。BLQ 规则是 `keep_as_censored`。

加载配置会检查设计。`PharmacologyStudy.run_nca` 在计算开始前停止，因为计算器尚未实现。
