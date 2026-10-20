#!/usr/bin/env python3

import argparse
import json
import math


R = 8.31446261815324


def load_case(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def residence_time(case):
    reactor = case["reactor"]
    return reactor["volume_m3"] / reactor["feed_volumetric_flow_m3_s"]


def forward_rate_constant(case, temperature_k):
    reaction = case["reaction"]
    return reaction["forward_pre_exponential_s^-1"] * math.exp(
        -reaction["forward_activation_energy_kj_kmol"] / (R * temperature_k)
    )


def equilibrium_constant(case, temperature_k):
    reaction = case["reaction"]
    delta_h = reaction["standard_enthalpy_kj_kmol"]
    reference_temperature = reaction["reference_temperature_k"]
    return reaction["equilibrium_constant_at_reference"] * math.exp(
        -delta_h / R * (1.0 / temperature_k - 1.0 / reference_temperature)
    )


def fixed_temperature_state(case, temperature_k):
    reactor = case["reactor"]
    tau = residence_time(case)
    k_forward = forward_rate_constant(case, temperature_k)
    k_reverse = k_forward / equilibrium_constant(case, temperature_k)

    a_feed = reactor["feed_a_kmol_m3"]
    b_feed = reactor["feed_b_kmol_m3"]
    a_concentration = (a_feed / tau + k_reverse * (b_feed + a_feed)) / (
        1.0 / tau + k_forward + k_reverse
    )
    b_concentration = b_feed + a_feed - a_concentration
    conversion = (a_feed - a_concentration) / a_feed
    net_rate = (
        k_forward * a_concentration - k_reverse * b_concentration
    )
    heat_generation = (
        -case["reaction"]["standard_enthalpy_kj_kmol"]
        * net_rate
        * reactor["volume_m3"]
    )
    sensible_heat = (
        reactor["feed_volumetric_flow_m3_s"]
        * reactor["volumetric_heat_capacity_kj_m3_k"]
        * (temperature_k - reactor["feed_temperature_k"])
    )
    utility_heat = reactor["overall_heat_transfer_kw_k"] * (
        temperature_k - reactor["utility_temperature_k"]
    )

    return {
        "temperature_k": temperature_k,
        "forward_rate_constant_s^-1": k_forward,
        "reverse_rate_constant_s^-1": k_reverse,
        "equilibrium_constant": equilibrium_constant(case, temperature_k),
        "a_concentration_kmol_m3": a_concentration,
        "b_concentration_kmol_m3": b_concentration,
        "conversion": conversion,
        "net_rate_kmol_m3_s": net_rate,
        "heat_generation_kw": heat_generation,
        "utility_heat_removal_kw": utility_heat,
        "energy_residual_kw": heat_generation - sensible_heat - utility_heat,
    }


def energy_residual_kw(case, temperature_k):
    state = fixed_temperature_state(case, temperature_k)
    return state["energy_residual_kw"]


def local_stability(case, temperature_k):
    reactor = case["reactor"]
    reaction = case["reaction"]
    tau = residence_time(case)
    volumetric_heat_capacity = reactor["volumetric_heat_capacity_kj_m3_k"]
    state = fixed_temperature_state(case, temperature_k)
    a_concentration = state["a_concentration_kmol_m3"]
    b_concentration = state["b_concentration_kmol_m3"]
    k_forward = state["forward_rate_constant_s^-1"]
    k_reverse = state["reverse_rate_constant_s^-1"]

    activation_energy = reaction["forward_activation_energy_kj_kmol"]
    delta_h = reaction["standard_enthalpy_kj_kmol"]
    d_k_forward = k_forward * activation_energy / (R * temperature_k**2)
    d_k_reverse = (
        k_reverse
        * (activation_energy - delta_h)
        / (R * temperature_k**2)
    )
    d_rate_d_temperature = (
        a_concentration * d_k_forward - b_concentration * d_k_reverse
    )
    heat_reaction = -delta_h

    j_ac = -1.0 / tau - (k_forward + k_reverse)
    j_at = -d_rate_d_temperature
    j_tc = heat_reaction * (-(k_forward + k_reverse)) / volumetric_heat_capacity
    j_tt = (
        -1.0 / tau
        + heat_reaction * d_rate_d_temperature / volumetric_heat_capacity
        - reactor["overall_heat_transfer_kw_k"]
        / (reactor["volume_m3"] * volumetric_heat_capacity)
    )
    trace = j_ac + j_tt
    determinant = j_ac * j_tt - j_at * j_tc
    return trace < 0.0 and determinant > 0.0


def find_steady_state_temperatures(case):
    scan = case["temperature_scan_k"]
    low = float(scan["minimum_k"])
    high = float(scan["maximum_k"])
    points = max(200001, int(scan.get("points", 4001)))
    span = high - low
    residual_values = []
    temperatures = []
    for index in range(points):
        temperature = low + span * index / (points - 1)
        temperatures.append(temperature)
        residual_values.append(energy_residual_kw(case, temperature))

    roots = []
    for index, temperature in enumerate(temperatures):
        if abs(residual_values[index]) <= 1.0e-10:
            roots.append(temperature)

    for index in range(points - 1):
        left_value = residual_values[index]
        right_value = residual_values[index + 1]
        if left_value * right_value >= 0.0:
            continue
        left = temperatures[index]
        right = temperatures[index + 1]
        for _ in range(90):
            middle = 0.5 * (left + right)
            middle_value = energy_residual_kw(case, middle)
            if left_value * middle_value <= 0.0:
                right = middle
            else:
                left = middle
                left_value = middle_value
        roots.append(0.5 * (left + right))

    deduplicated = []
    for root in sorted(roots):
        if not deduplicated or abs(root - deduplicated[-1]) > max(
            1.0e-5, span / points * 1.0e-3
        ):
            deduplicated.append(root)
    return deduplicated


def steady_state_report(case, temperature_k):
    state = fixed_temperature_state(case, temperature_k)
    reactor = case["reactor"]
    reaction = case["reaction"]
    heat_generation = (
        -reaction["standard_enthalpy_kj_kmol"]
        * state["net_rate_kmol_m3_s"]
        * reactor["volume_m3"]
    )
    utility_heat = reactor["overall_heat_transfer_kw_k"] * (
        temperature_k - reactor["utility_temperature_k"]
    )
    return {
        "temperature_k": temperature_k,
        "a_concentration_kmol_m3": state["a_concentration_kmol_m3"],
        "b_concentration_kmol_m3": state["b_concentration_kmol_m3"],
        "conversion": state["conversion"],
        "net_rate_kmol_m3_s": state["net_rate_kmol_m3_s"],
        "heat_generation_kw": heat_generation,
        "utility_heat_removal_kw": utility_heat,
        "stable": local_stability(case, temperature_k),
    }


def build_report(case):
    point_checks = [
        fixed_temperature_state(case, temperature_k)
        for temperature_k in case.get("point_check_temperatures_k", [])
    ]
    temperatures = find_steady_state_temperatures(case)
    steady_states = [steady_state_report(case, temperature_k) for temperature_k in temperatures]
    return {
        "case_id": case.get("case_id", ""),
        "point_checks": point_checks,
        "steady_states": steady_states,
        "multiple_steady_states": len(steady_states) > 1,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True)
    args = parser.parse_args()
    report = build_report(load_case(args.case))
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
