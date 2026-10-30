# 带噪声的合成单次静脉推注

一名受试者，剂量 1 mg，CL = 0.40 L/h，V = 2.0 L。定量浓度由 `C(t) * exp(0.10 * z)` 抽一次（`random.Random(20261030)`）并保留 8 位小数。24 h 低于 0.01 µg/mL，以 BLQ 保存，`numeric_value` 为空。CL 和 V 是绝对值，不是按体重归一的参数。

这不是动物实验，也不能用于给药。加权非线性最小二乘应在 10% 相对误差内回收 CL 和 V。重新生成：

```bash
python -m examples.iv_bolus_noisy
```
