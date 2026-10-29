# 验证

0.0.1 验证的是数据契约和参考方程。它不验证一个 NCA 实现，因为计算器还不存在。

## 已自动执行的检查

- 单位维度和换算。
- 单条记录规则：时区、BLQ、零浓度、输注时长、CL 与 CL/F。
- 关系规则：外键、时钟、LLOQ、基质、未审核身份链接、临床数据的研究用提示。
- BLQ 策略不修改原观测；有偏替换必须 `acknowledge_bias=true`。
- 合成研究在 strict 模式下零错误、零警告。
- CSV、TSV、Parquet 往返后数值和校验结果保持一致。
- 结果包的文件哈希、`research_only`、以及 NCA 产物缺席。
- 证据夹具符合证据交换 Schema；错误升级会被拒绝。
- 一室静脉推注解析解和未来 NCA 的靶值。

本地与 CI 运行：

```bash
pip install -e ".[dev]"
ruff check src tests examples
mypy src/pharmacoweave
pytest
```

CI 在 Python 3.11 和 3.12 上执行 lint、类型检查和测试。

## 合成基准

`examples/data/iv_bolus_single_dose/`：

| 量 | 值 |
|---|---|
| 给药 | 1 mg/kg 静脉推注 |
| CL | 0.20 L/h/kg |
| V | 1.0 L/kg |
| 受试者体重 | 0.25 kg 与 0.30 kg |
| C(0) | 1 µg/mL |
| t½ | ln(2) / 0.2 h |
| AUC∞ | 5 µg·h/mL |
| LLOQ | 0.01 µg/mL |
| 24 h | 低于 LLOQ，状态 BLQ，无数值 |

`tests/numerical_benchmarks/test_iv_bolus_truth.py` 锁定这些靶值。将来的 NCA 应在量化点上、并把 24 h 保持为删失之后，接近 Dose/CL 和 ln(2)·V/CL。允许误差要在实现该计算器时单独写明，不能用本文件的解析解代替梯形积分测试。

## 第一版 NCA 还要补上的测试

1. 选定并写死终末相规则，包括最少点数、调整 R² 和外推比例的处理方式。
2. 用解析曲线的量化点回收 AUC_last、AUC∞、Cmax、Tmax、λz 和 t½，并记录与解析真值的容差。
3. 把 24 h BLQ 放进输入，确认它没有被当成 0 或当成一个定量点。
4. 不规则时间、终末相点数不足、外推比例过高、单位错误和优化式 NCA 不适用的设计都要失败或省略相应参数。
5. 对照一份外部参考输出的数字，不把该参考实现作为库依赖引进来。

房室拟合的第一份测试应是：用上述解析曲线加已知误差模型，回收 CL 和 V；不可辨识或不收敛时不得标成 pass。

## 明确不做的符合性声明

字段设计参考 ICH M10 的定量限和 QC 概念，以及 ICH M15 的 question of interest 与 context of use。测试通过不等于生物分析方法学验证，也不等于 MIDD 证据已经被监管机构接受。CDISC SEND 导出、21 CFR Part 11 审计追踪和人类受试者访问控制不在本版本。
