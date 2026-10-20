#!/usr/bin/env python3

import copy
import json
import math
import os
import subprocess
import tempfile
import traceback


R = 8.31446261815324
APP_PATH = os.environ.get("APP_PATH", "/app/reactor_report.py")
LOG_DIR = os.environ.get("LOG_DIR", "/logs/verifier")
REWARD_FILE = os.path.join(LOG_DIR, "reward.txt")
REPORT_FILE = os.path.join(LOG_DIR, "report.json")


BASE_CASE = {
    "case_id": "base-isomerization",
    "reaction": {
        "forward_pre_exponential_s^-1": 2.2e13,
        "forward_activation_energy_kj_kmol": 100000.0,
        "equilibrium_constant_at_reference": 8.3,
        "reference_temperature_k": 350.0,
        "standard_enthalpy_kj_kmol": -90000.0,
    },
    "reactor": {
        "volume_m3": 1.0,
        "feed_volumetric_flow_m3_s": 0.00909090909090909,
        "feed_a_kmol_m3": 2.35,
        "feed_b_kmol_m3": 0.0,
        "feed_temperature_k": 303.3,
        "utility_temperature_k": 299.0,
        "overall_heat_transfer_kw_k": 0.58,
        "volumetric_heat_capacity_kj_m3_k": 1300.0,
    },
    "point_check_temperatures_k": [300.0, 325.0, 350.0, 375.0],
    "temperature_scan_k": {
        "minimum_k": 275.0,
        "maximum_k": 425.0,
        "points": 4001,
    },
}


def hidden_point_case():
    case = copy.deepcopy(BASE_CASE)
    case["case_id"] = "hidden-point-case"
    case["reactor"].update(
        {
            "volume_m3": 0.8,
            "feed_volumetric_flow_m3_s": 0.0125,
            "feed_b_kmol_m3": 0.35,
        }
    )
    case["point_check_temperatures_k"] = [302.0, 337.0, 368.0]
    return case


def hidden_thermal_case():
    case = copy.deepcopy(BASE_CASE)
    case["case_id"] = "hidden-thermal-case"
    case["reactor"].update(
        {
            "volume_m3": 1.2,
            "feed_volumetric_flow_m3_s": 0.010,
            "feed_a_kmol_m3": 2.0,
            "feed_b_kmol_m3": 0.25,
            "feed_temperature_k": 301.5,
            "utility_temperature_k": 295.0,
            "overall_heat_transfer_kw_k": 16.0,
            "volumetric_heat_capacity_kj_m3_k": 1550.0,
        }
    )
    case["temperature_scan_k"] = {
        "minimum_k": 280.0,
        "maximum_k": 430.0,
        "points": 4001,
    }
    return case


def hidden_multiplicity_case():
    case = copy.deepcopy(BASE_CASE)
    case["case_id"] = "hidden-near-bifurcation-case"
    case["reaction"].update(
        {
            "forward_pre_exponential_s^-1": 180144209626732.0,
            "forward_activation_energy_kj_kmol": 109974.82159025133,
            "equilibrium_constant_at_reference": 1.6817307350817081,
            "standard_enthalpy_kj_kmol": -98773.91979538265,
        }
    )
    case["reactor"].update(
        {
            "volume_m3": 2.633597382858819,
            "feed_volumetric_flow_m3_s": 0.007055473449685984,
            "feed_a_kmol_m3": 2.00002406478424,
            "feed_b_kmol_m3": 0.012818651734272657,
            "feed_temperature_k": 301.6585775611144,
            "utility_temperature_k": 285.9074345504981,
            "overall_heat_transfer_kw_k": 1.406822,
            "volumetric_heat_capacity_kj_m3_k": 1973.4725060300707,
        }
    )
    case["temperature_scan_k"] = {
        "minimum_k": 270.0,
        "maximum_k": 430.0,
        "points": 17,
    }
    return case


def equilibrium_constant(case, temperature_k):
    reaction = case["reaction"]
    return reaction["equilibrium_constant_at_reference"] * math.exp(
        -reaction["standard_enthalpy_kj_kmol"]
        / R
        * (
            1.0 / temperature_k
            - 1.0 / reaction["reference_temperature_k"]
        )
    )


def expected_fixed_state(case, temperature_k):
    reaction = case["reaction"]
    reactor = case["reactor"]
    tau = reactor["volume_m3"] / reactor["feed_volumetric_flow_m3_s"]
    k_forward = reaction["forward_pre_exponential_s^-1"] * math.exp(
        -reaction["forward_activation_energy_kj_kmol"] / (R * temperature_k)
    )
    k_reverse = k_forward / equilibrium_constant(case, temperature_k)
    a_feed = reactor["feed_a_kmol_m3"]
    b_feed = reactor["feed_b_kmol_m3"]
    a_concentration = (a_feed / tau + k_reverse * (b_feed + a_feed)) / (
        1.0 / tau + k_forward + k_reverse
    )
    b_concentration = b_feed + a_feed - a_concentration
    conversion = (a_feed - a_concentration) / a_feed
    net_rate = k_forward * a_concentration - k_reverse * b_concentration
    heat_generation = (
        -case["reaction"]["standard_enthalpy_kj_kmol"]
        * net_rate
        * case["reactor"]["volume_m3"]
    )
    sensible_heat = (
        case["reactor"]["feed_volumetric_flow_m3_s"]
        * case["reactor"]["volumetric_heat_capacity_kj_m3_k"]
        * (temperature_k - case["reactor"]["feed_temperature_k"])
    )
    utility_heat = case["reactor"]["overall_heat_transfer_kw_k"] * (
        temperature_k - case["reactor"]["utility_temperature_k"]
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


def expected_residual(case, temperature_k):
    state = expected_fixed_state(case, temperature_k)
    return state["energy_residual_kw"]


def expected_roots(case):
    scan = case["temperature_scan_k"]
    low = scan["minimum_k"]
    high = scan["maximum_k"]
    points = max(200001, int(scan["points"]))
    temperatures = [
        low + (high - low) * index / (points - 1)
        for index in range(points)
    ]
    residuals = [expected_residual(case, value) for value in temperatures]
    roots = []
    for index in range(points - 1):
        left_value = residuals[index]
        right_value = residuals[index + 1]
        if left_value * right_value >= 0.0:
            continue
        left = temperatures[index]
        right = temperatures[index + 1]
        for _ in range(90):
            middle = 0.5 * (left + right)
            middle_value = expected_residual(case, middle)
            if left_value * middle_value <= 0.0:
                right = middle
            else:
                left = middle
                left_value = middle_value
        roots.append(0.5 * (left + right))
    return roots


def expected_stability(case, temperature_k):
    reactor = case["reactor"]
    reaction = case["reaction"]
    state = expected_fixed_state(case, temperature_k)
    tau = reactor["volume_m3"] / reactor["feed_volumetric_flow_m3_s"]
    capacity = reactor["volumetric_heat_capacity_kj_m3_k"]
    k_forward = state["forward_rate_constant_s^-1"]
    k_reverse = state["reverse_rate_constant_s^-1"]
    a_concentration = state["a_concentration_kmol_m3"]
    b_concentration = state["b_concentration_kmol_m3"]
    activation_energy = reaction["forward_activation_energy_kj_kmol"]
    delta_h = reaction["standard_enthalpy_kj_kmol"]
    rate_temperature_derivative = (
        a_concentration
        * k_forward
        * activation_energy
        / (R * temperature_k**2)
        - b_concentration
        * k_reverse
        * (activation_energy - delta_h)
        / (R * temperature_k**2)
    )
    heat_reaction = -delta_h
    j_ac = -1.0 / tau - (k_forward + k_reverse)
    j_at = -rate_temperature_derivative
    j_tc = (
        heat_reaction
        * (-(k_forward + k_reverse))
        / capacity
    )
    j_tt = (
        -1.0 / tau
        + heat_reaction * rate_temperature_derivative / capacity
        - reactor["overall_heat_transfer_kw_k"]
        / (reactor["volume_m3"] * capacity)
    )
    trace = j_ac + j_tt
    determinant = j_ac * j_tt - j_at * j_tc
    return trace < 0.0 and determinant > 0.0


def invoke(case):
    with tempfile.TemporaryDirectory(prefix="cstr-case-") as directory:
        case_path = os.path.join(directory, "case.json")
        with open(case_path, "w", encoding="utf-8") as handle:
            json.dump(case, handle)
        completed = subprocess.run(
            ["python3", APP_PATH, "--case", case_path],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    if completed.returncode != 0:
        raise AssertionError(
            "process failed with code "
            f"{completed.returncode}: {completed.stderr[-2000:]}"
        )
    try:
        return json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise AssertionError(
            f"standard output is not one JSON object: {completed.stdout[-2000:]}"
        ) from error


def close(actual, expected, *, abs_tolerance, rel_tolerance=1.0e-8):
    if not math.isclose(
        float(actual),
        float(expected),
        abs_tol=abs_tolerance,
        rel_tol=rel_tolerance,
    ):
        raise AssertionError(f"actual={actual!r}, expected={expected!r}")


def point_checks_pass():
    case = hidden_point_case()
    report = invoke(case)
    actual_points = report["point_checks"]
    expected_points = [
        expected_fixed_state(case, value)
        for value in case["point_check_temperatures_k"]
    ]
    if len(actual_points) != len(expected_points):
        raise AssertionError(
            f"expected {len(expected_points)} point checks, got {len(actual_points)}"
        )
    for actual, expected in zip(actual_points, expected_points):
        close(actual["temperature_k"], expected["temperature_k"], abs_tolerance=1.0e-9)
        close(
            actual["forward_rate_constant_s^-1"],
            expected["forward_rate_constant_s^-1"],
            abs_tolerance=1.0e-12,
            rel_tolerance=2.0e-5,
        )
        close(
            actual["reverse_rate_constant_s^-1"],
            expected["reverse_rate_constant_s^-1"],
            abs_tolerance=1.0e-12,
            rel_tolerance=2.0e-5,
        )
        close(
            actual["equilibrium_constant"],
            expected["equilibrium_constant"],
            abs_tolerance=1.0e-9,
            rel_tolerance=2.0e-5,
        )
        close(
            actual["a_concentration_kmol_m3"],
            expected["a_concentration_kmol_m3"],
            abs_tolerance=2.0e-5,
        )
        close(
            actual["b_concentration_kmol_m3"],
            expected["b_concentration_kmol_m3"],
            abs_tolerance=2.0e-5,
        )
        close(
            actual["conversion"],
            expected["conversion"],
            abs_tolerance=1.0e-5,
        )


def thermal_balance_pass():
    case = hidden_thermal_case()
    report = invoke(case)
    actual_points = report["point_checks"]
    if len(actual_points) != 4:
        raise AssertionError(
            f"expected four heat-balance checks, got {len(actual_points)}"
        )
    reactor = case["reactor"]
    reaction = case["reaction"]
    for actual in actual_points:
        temperature = actual["temperature_k"]
        heat_generation = (
            -reaction["standard_enthalpy_kj_kmol"]
            * actual["net_rate_kmol_m3_s"]
            * reactor["volume_m3"]
        )
        sensible_heat = (
            reactor["feed_volumetric_flow_m3_s"]
            * reactor["volumetric_heat_capacity_kj_m3_k"]
            * (temperature - reactor["feed_temperature_k"])
        )
        utility_heat = reactor["overall_heat_transfer_kw_k"] * (
            temperature - reactor["utility_temperature_k"]
        )
        close(
            actual["heat_generation_kw"],
            heat_generation,
            abs_tolerance=5.0e-3,
            rel_tolerance=2.0e-4,
        )
        close(
            actual["utility_heat_removal_kw"],
            utility_heat,
            abs_tolerance=5.0e-3,
            rel_tolerance=2.0e-4,
        )
        close(
            actual["energy_residual_kw"],
            heat_generation - sensible_heat - utility_heat,
            abs_tolerance=5.0e-3,
            rel_tolerance=2.0e-4,
        )


def multiplicity_pass():
    case = hidden_multiplicity_case()
    report = invoke(case)
    expected_temperatures = expected_roots(case)
    if len(expected_temperatures) != 3:
        raise AssertionError("test case construction did not produce three roots")
    actual_states = report["steady_states"]
    if len(actual_states) != len(expected_temperatures):
        raise AssertionError(
            f"expected {len(expected_temperatures)} steady states, "
            f"got {len(actual_states)}"
        )
    if report["multiple_steady_states"] is not True:
        raise AssertionError("multiple_steady_states must be true")
    for actual, expected_temperature in zip(actual_states, expected_temperatures):
        expected = expected_fixed_state(case, expected_temperature)
        close(
            actual["temperature_k"],
            expected_temperature,
            abs_tolerance=5.0e-3,
        )
        close(
            actual["conversion"],
            expected["conversion"],
            abs_tolerance=2.0e-4,
        )


def stability_pass():
    case = hidden_multiplicity_case()
    report = invoke(case)
    expected_temperatures = expected_roots(case)
    actual_states = report["steady_states"]
    if len(actual_states) != len(expected_temperatures):
        raise AssertionError("steady-state count does not permit stability comparison")
    expected_labels = [
        expected_stability(case, temperature)
        for temperature in expected_temperatures
    ]
    actual_labels = [state.get("stable") for state in actual_states]
    if actual_labels != expected_labels:
        raise AssertionError(
            f"actual stability labels={actual_labels!r}, "
            f"expected={expected_labels!r}"
        )


def write_results(results, score):
    os.makedirs(LOG_DIR, exist_ok=True)
    with open(REWARD_FILE, "w", encoding="utf-8") as handle:
        handle.write(f"{score:.3f}\n")
    with open(REPORT_FILE, "w", encoding="utf-8") as handle:
        json.dump(
            {"score": score, "checks": results},
            handle,
            ensure_ascii=False,
            indent=2,
        )


def main():
    checks = [
        ("temperature_dependent_reversible_kinetics", 0.30, point_checks_pass),
        ("coupled_heat_balance", 0.30, thermal_balance_pass),
        ("all_steady_states", 0.20, multiplicity_pass),
        ("local_stability", 0.20, stability_pass),
    ]
    results = []
    score = 0.0
    for name, weight, check in checks:
        try:
            check()
            passed = True
            detail = "passed"
            score += weight
        except Exception as error:
            passed = False
            detail = f"{type(error).__name__}: {error}"
        results.append(
            {
                "name": name,
                "weight": weight,
                "passed": passed,
                "detail": detail,
            }
        )
    score = round(score, 3)
    write_results(results, score)
    print(json.dumps({"score": score, "checks": results}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        os.makedirs(LOG_DIR, exist_ok=True)
        with open(REWARD_FILE, "w", encoding="utf-8") as handle:
            handle.write("0\n")
        with open(REPORT_FILE, "w", encoding="utf-8") as handle:
            handle.write(traceback.format_exc())
        raise
