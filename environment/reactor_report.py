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
    activation_energy = reaction["forward_activation_energy_kj_kmol"]
    pre_exponential = reaction["forward_pre_exponential_s^-1"]
    return pre_exponential * math.exp(-activation_energy / (R * temperature_k))


def equilibrium_constant(case, temperature_k):
    return case["reaction"]["equilibrium_constant_at_reference"]


def fixed_temperature_state(case, temperature_k):
    reaction = case["reaction"]
    reactor = case["reactor"]
    tau = residence_time(case)
    k_forward = forward_rate_constant(case, temperature_k)
    k_reverse = k_forward / equilibrium_constant(case, temperature_k)

    a_feed = reactor["feed_a_kmol_m3"]
    b_feed = reactor["feed_b_kmol_m3"]
    denominator = 1.0 / tau + k_forward + k_reverse
    a_concentration = (a_feed / tau + k_reverse * b_feed) / denominator
    b_concentration = b_feed + a_feed - a_concentration
    conversion = (a_feed - a_concentration) / a_feed
    net_rate = k_forward * a_concentration - k_reverse * b_concentration
    heat_generation = (
        -reaction["standard_enthalpy_kj_kmol"]
        * net_rate
        * reactor["volume_m3"]
    )
    sensible_heat = (
        reactor["feed_volumetric_flow_m3_s"]
        * reactor["volumetric_heat_capacity_kj_m3_k"]
        * (temperature_k - reactor["feed_temperature_k"])
    )
    utility_heat = reactor["overall_heat_transfer_kw_k"] * (
        temperature_k - reactor["feed_temperature_k"]
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
    return None


def find_steady_state_temperatures(case):
    scan = case["temperature_scan_k"]
    low = scan["minimum_k"]
    high = scan["maximum_k"]
    guess = 0.5 * (low + high)

    for _ in range(80):
        residual = energy_residual_kw(case, guess)
        step = max(1.0e-3, abs(guess) * 1.0e-6)
        slope = (
            energy_residual_kw(case, guess + step) - energy_residual_kw(case, guess - step)
        ) / (2.0 * step)
        if abs(slope) < 1.0e-14:
            break
        next_guess = guess - residual / slope
        if not low <= next_guess <= high:
            break
        if abs(next_guess - guess) < 1.0e-7:
            guess = next_guess
            break
        guess = next_guess

    if low <= guess <= high and abs(energy_residual_kw(case, guess)) < 1.0e-5:
        return [guess]
    return []


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
        temperature_k - reactor["feed_temperature_k"]
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
