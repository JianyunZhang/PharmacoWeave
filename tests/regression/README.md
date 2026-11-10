# 回归测试

本目录保存三份单次给药 NCA 的 50 位十进制独立手算冻结参考，均按 linear-up/log-down 梯形法则与 `pharmacoweave.lambda_z.adjusted_r2.v1` 计算，不导入 PharmacoWeave、PKNCA 或 Phoenix。另有一份尿药排泄手算参考，只做浓度乘体积，不使用血浆梯形积分：

- `iv_bolus_quiet_nca.json`：无噪声单次静脉推注 NCA 冻结参考（含对数线性回推 `C0`），由 `test_iv_bolus_quiet_reference.py` 校验。
- `oral_quiet_nca.json`：无噪声单次口服给药 NCA 冻结参考（从 `(0, 0)` 起点积分，终末相候选窗口限制在 `t > Tmax`，输出表观参数 `CL/F` 与 `Vz/F`），由 `test_oral_quiet_reference.py` 校验。
- `iv_infusion_quiet_nca.json`：无噪声单次静脉恒速输注 NCA 冻结参考（显式校验输注时长，从 `(0, 0)` 积分，终末相候选窗口限制在 `t >= max(Tmax, infusion_duration)`），由 `test_iv_infusion_quiet_reference.py` 校验。
- `urine_excretion_quiet.json`：两段尿液采集的排泄量与肾清除手算参考（2.5 mg、2.0 mg、总量 4.5 mg、血浆 AUC 9 mg·h/L、肾清除 0.5 L/h），由 `test_urine_excretion_nca.py` 校验。

相对误差上限见 `docs/validation.md`。各例子中的 BLQ 样品均记录在 `omitted_observations` 中，不作为 `0` 参与积分。
