# 放热可逆反应的连续搅拌釜工况分析

中试装置中有一台连续搅拌釜，液相反应物 A 在釜内可逆转化为 B。运行组需要在开车前核对不同温度下的反应状态，并找出在设计换热条件下可能存在的全部稳态，供后续控制和联锁方案使用。

现有 `reactor_report.py` 只是早期的筛选版本，遇到温度变化、换热器给水温度和多重稳态时会给出不可靠的结果。请修改它，使它按下面的输入输出约定完成计算。

## 模型和单位

- 反应为液相基元可逆反应 `A <-> B`，体积流量、体积和液体密度恒定，忽略相变和混合热。
- A 的消耗速率采用 `r = k_f C_A - k_r C_B`，单位为 `kmol/(m^3 s)`；正反应速率常数符合 Arrhenius 关系，逆反应速率常数由同一温度下的反应平衡决定。
- 平衡常数随温度变化，采用恒定的标准反应焓计算；输入中的活化能和标准反应焓单位均为 `kJ/kmol`。
- 反应器物料衡算采用稳态连续搅拌釜假设；反应器热量衡算同时考虑进料显热和换热器移除的热量。
- 局部稳定性应根据未施加控制的釜内 A 浓度和温度所满足的非稳态物料、能量方程进行判断。

输入 JSON 的字段和单位如下：

```json
{
  "case_id": "example",
  "reaction": {
    "forward_pre_exponential_s^-1": 2.2e13,
    "forward_activation_energy_kj_kmol": 100000.0,
    "equilibrium_constant_at_reference": 8.3,
    "reference_temperature_k": 350.0,
    "standard_enthalpy_kj_kmol": -90000.0
  },
  "reactor": {
    "volume_m3": 1.0,
    "feed_volumetric_flow_m3_s": 0.00909090909090909,
    "feed_a_kmol_m3": 2.35,
    "feed_b_kmol_m3": 0.0,
    "feed_temperature_k": 303.3,
    "utility_temperature_k": 299.0,
    "overall_heat_transfer_kw_k": 0.58,
    "volumetric_heat_capacity_kj_m3_k": 1300.0
  },
  "point_check_temperatures_k": [300, 325, 350, 375],
  "temperature_scan_k": {"minimum_k": 275, "maximum_k": 425, "points": 4001}
}
```

## 需要交付的内容

继续使用 `python3 reactor_report.py --case <case.json>` 运行程序。程序把结果 JSON 写到标准输出，不能依赖第三方 Python 包，也不能把输入案例的数值写死在代码里。

输出必须包含：

1. `point_checks`：对 `point_check_temperatures_k` 中的每个温度，给出该温度下定温稳态的
   - `temperature_k`
   - `forward_rate_constant_s^-1`
   - `reverse_rate_constant_s^-1`
   - `equilibrium_constant`
   - `a_concentration_kmol_m3`
   - `b_concentration_kmol_m3`
   - `conversion`
   - `net_rate_kmol_m3_s`
   - `heat_generation_kw`
   - `utility_heat_removal_kw`
   - `energy_residual_kw`

   `energy_residual_kw` 等于反应放热减去进料显热和换热介质移除的热量；它不是要求在这些检查温度处等于零。
2. `steady_states`：在 `temperature_scan_k` 范围内找出所有满足物料衡算和热量衡算的温度点。每个温度点给出
   - `temperature_k`
   - `a_concentration_kmol_m3`
   - `b_concentration_kmol_m3`
   - `conversion`
   - `net_rate_kmol_m3_s`
   - `heat_generation_kw`
   - `utility_heat_removal_kw`
   - `stable`
3. `multiple_steady_states`：当且仅当扫描范围内存在两个或更多稳态时为 `true`。

`steady_states` 按温度从低到高排列。数值解要足够精确，避免把同一个根重复列出；工况可能只有一个稳态，也可能有多个稳态。

`temperature_scan_k.points` 只是初始采样密度，不能保证把彼此非常接近的稳态分开；程序需要自行细化，避免漏掉任何稳态。

## 约束

- 所有物理量都使用上述单位。
- 温度必须使用热力学温度计算指数项。
- 组分的稳态浓度必须与进料组成、停留时间和可逆反应速率同时一致。
- 换热器热量使用反应器温度与换热介质温度之差计算，不能把进料温度当作换热介质温度。
- 稳定性判断应对应未控制的局部线性化行为，而不能只依据稳态转化率高低。
- 只修改程序实现；输入文件中的参数是工况数据，不应被改写。
