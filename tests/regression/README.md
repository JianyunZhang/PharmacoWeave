# 回归测试

`iv_bolus_quiet_nca.json` 是无噪声单次静脉推注 NCA 的冻结参考。它用 50 位十进制运算按 linear-up/log-down 和 `pharmacoweave.lambda_z.adjusted_r2.v1` 手算，不导入 PharmacoWeave、PKNCA 或 Phoenix。相对误差上限见 `docs/validation.md`，并由 `test_iv_bolus_quiet_reference.py` 再写死一遍。

24 h 的 BLQ 记在 `omitted_observations` 里，不作为 0 参加积分。
