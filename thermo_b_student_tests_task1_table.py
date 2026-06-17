"""
Public test runner for the Thermodynamics B AI test-bench notebook.

This file intentionally keeps the detailed assertion logic outside the notebook
so the notebook behaves more like an engineering test bench. These are still
public tests: students who open this file can inspect them. For summative
assessment, keep additional private tests separate and do not distribute them.
"""

from __future__ import annotations

import math
from typing import Callable, Mapping, Any

import numpy as np
import pandas as pd


class PublicTestFailure(AssertionError):
    """Raised when a public notebook test fails."""


def _fail(message: str) -> None:
    raise PublicTestFailure(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _close(actual: float, expected: float, *, rtol: float = 1e-7, atol: float = 1e-12, message: str = "") -> None:
    if not np.isclose(actual, expected, rtol=rtol, atol=atol):
        _fail(message or f"Numerical check failed: expected approximately {expected}, got {actual}.")


def _array_close(actual: Any, expected: Any, *, rtol: float = 1e-7, atol: float = 1e-12, message: str = "") -> None:
    if not np.allclose(actual, expected, rtol=rtol, atol=atol):
        _fail(message or "Array-like numerical check failed.")


def _get(ns: Mapping[str, Any], name: str) -> Any:
    if name not in ns:
        _fail(f"Missing required name `{name}`. Run the cell that defines it before running this test.")
    return ns[name]


def _callable(ns: Mapping[str, Any], name: str) -> Callable:
    obj = _get(ns, name)
    if not callable(obj):
        _fail(f"`{name}` exists, but it is not callable. It should be a function.")
    return obj


def _test_ideal_gas_pressure(ns: Mapping[str, Any]) -> None:
    ideal_gas_pressure = _callable(ns, "ideal_gas_pressure")
    _close(ideal_gas_pressure(T=300, V=30, R=8.314), 83.14, message="The direct ideal-gas calculation is not correct.")

    p_small_v = ideal_gas_pressure(T=300, V=30, R=8.314)
    p_large_v = ideal_gas_pressure(T=300, V=300, R=8.314)
    _require(p_small_v > p_large_v, "At fixed temperature, pressure should decrease as molar volume increases.")

    p_low_t = ideal_gas_pressure(T=250, V=100, R=8.314)
    p_high_t = ideal_gas_pressure(T=350, V=100, R=8.314)
    _require(p_high_t > p_low_t, "At fixed molar volume, pressure should increase as temperature increases.")


def _test_ideal_gas_table(ns: Mapping[str, Any]) -> None:
    build_ideal_gas_table = _callable(ns, "build_ideal_gas_table")

    # First check a small toy table so the function cannot be hard-coded only for the full exercise data.
    toy_volumes = [10, 20]
    toy_temperatures = [100.0, 200.0, 300.0]
    toy_df = build_ideal_gas_table(8.314, toy_volumes, toy_temperatures)
    _require(isinstance(toy_df, pd.DataFrame), "`build_ideal_gas_table` should return a pandas DataFrame.")
    _require(list(toy_df.index) == toy_volumes, "The toy table index should preserve the provided volume order.")
    _require(list(toy_df.columns) == toy_temperatures, "The toy table columns should preserve the provided temperature order.")
    _close(toy_df.loc[10, 100.0], 8.314 * 100.0 / 10, message="A toy-table pressure value is not correct.")
    _close(toy_df.loc[20, 300.0], 8.314 * 300.0 / 20, message="A toy-table pressure value is not correct.")

    df = _get(ns, "ideal_gas_law_df")
    _require(isinstance(df, pd.DataFrame), "`ideal_gas_law_df` should be a pandas DataFrame.")

    expected_temperatures = [255.2, 265.1, 274.4, 284.0, 294.4, 304.1, 334.1, 354.1]
    _require(df.shape == (570, 8), "The ideal-gas table has the wrong shape. Check the volume and temperature ranges.")
    _require(list(df.columns) == expected_temperatures, "The temperature columns do not match the required values from the protocol.")
    _require(df.index[0] == 30 and df.index[-1] == 599, "The volume index should run from 30 to 599 cm^3/mol.")
    _require(not df.isna().any().any(), "The ideal-gas table still contains missing values.")

    try:
        fixed_t_low_v = df.loc[30, 304.1]
        fixed_t_high_v = df.loc[599, 304.1]
        fixed_v_high_t = df.loc[100, 354.1]
        fixed_v_low_t = df.loc[100, 255.2]
    except Exception as exc:
        _fail(f"Could not access expected rows/columns in `ideal_gas_law_df`: {exc}")

    _close(df.loc[30, 255.2], 8.314 * 255.2 / 30, message="The low-temperature, low-volume table value is incorrect.")
    _close(df.loc[599, 354.1], 8.314 * 354.1 / 599, message="The high-temperature, high-volume table value is incorrect.")
    _require(fixed_t_low_v > fixed_t_high_v, "At fixed temperature, the pressure trend with volume is not physically consistent.")
    _require(fixed_v_high_t > fixed_v_low_t, "At fixed volume, the pressure trend with temperature is not physically consistent.")


def _test_calculate_kappa(ns: Mapping[str, Any]) -> None:
    calculate_kappa = _callable(ns, "calculate_kappa")
    omega = 0.225
    expected = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    _close(calculate_kappa(omega), expected, message="The kappa value is not correct for a scalar acentric factor.")

    omegas = np.array([0.0, 0.1, 0.2])
    expected_array = 0.37464 + 1.54226 * omegas - 0.26992 * omegas**2
    _array_close(calculate_kappa(omegas), expected_array, message="The kappa function should also work for NumPy arrays.")


def _test_calculate_alpha(ns: Mapping[str, Any]) -> None:
    calculate_alpha = _callable(ns, "calculate_alpha")
    kappa = 0.7
    _close(calculate_alpha(T=300, Tc=300, kappa=kappa), 1.0, message="Alpha should be 1 when T equals Tc.")

    alpha_low_T = calculate_alpha(T=250, Tc=300, kappa=kappa)
    alpha_high_T = calculate_alpha(T=350, Tc=300, kappa=kappa)
    _require(alpha_low_T > 1.0, "For this positive-kappa example, alpha should be above 1 below Tc.")
    _require(alpha_high_T < 1.0, "For this positive-kappa example, alpha should be below 1 above Tc.")


def _test_calculate_ai_bi(ns: Mapping[str, Any]) -> None:
    calculate_ai = _callable(ns, "calculate_ai")
    calculate_bi = _callable(ns, "calculate_bi")
    R = 8.314
    T = 304.1
    Tc = 304.2
    Pc = 7.38
    omega = 0.225
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    alpha = (1 + kappa * (1 - np.sqrt(T / Tc)))**2
    expected_ai = 0.45724 * (R**2 * Tc**2 / Pc) * alpha
    expected_bi = 0.07780 * R * Tc / Pc
    _close(calculate_ai(R, T, Tc, Pc, omega), expected_ai, message="The pure-component attraction parameter is not correct.")
    _close(calculate_bi(R, Tc, Pc), expected_bi, message="The pure-component co-volume parameter is not correct.")
    _require(calculate_ai(R, 280, Tc, Pc, omega) > calculate_ai(R, 360, Tc, Pc, omega), "The attraction parameter should decrease as temperature increases for this example.")


def _test_peng_robinson_pressure(ns: Mapping[str, Any]) -> None:
    peng_robinson_pressure = _callable(ns, "peng_robinson_pressure")
    T = 304.1
    Tc = 304.2
    Pc = 7.38
    omega = 0.225
    R = 8.314
    V_large = 1_000_000
    p_pr = peng_robinson_pressure(T, V_large, Tc, Pc, omega, R=R)
    p_ideal = R * T / V_large
    _close(p_pr, p_ideal, rtol=1e-3, message="At very large molar volume, PR should approach ideal-gas behavior.")
    _require(np.isfinite(p_pr), "The PR pressure function returned a non-finite value.")


def _test_calculate_A_B(ns: Mapping[str, Any]) -> None:
    calculate_A = _callable(ns, "calculate_A")
    calculate_B = _callable(ns, "calculate_B")
    a = 1000.0
    b = 25.0
    P = 2.0
    R = 8.314
    T = 300.0
    _close(calculate_A(a, P, R, T), a * P / (R * T) ** 2, message="Dimensionless A is not correct.")
    _close(calculate_B(b, P, R, T), b * P / (R * T), message="Dimensionless B is not correct.")


def _test_pr_cubic_coefficients(ns: Mapping[str, Any]) -> None:
    pr_cubic_coefficients = _callable(ns, "pr_cubic_coefficients")
    A = 0.5
    B = 0.1
    coeffs = pr_cubic_coefficients(A, B)
    expected = [1.0, -(1 - B), A - 3 * B**2 - 2 * B, -(A * B - B**2 - B**3)]
    _array_close(coeffs, expected, message="The Peng-Robinson cubic coefficients are not correct.")


def _test_roots_and_selection(ns: Mapping[str, Any]) -> None:
    real_roots_only = _callable(ns, "real_roots_only")
    select_pr_root = _callable(ns, "select_pr_root")
    roots = np.array([0.1 + 0j, 0.5 + 1e-10j, 1.2 + 1e-4j])
    real_roots = real_roots_only(roots, atol=1e-8)
    _array_close(real_roots, [0.1, 0.5], message="The real-root filtering does not handle tiny imaginary parts correctly.")
    _close(select_pr_root([0.1, 0.5, 0.9], phase="vapour"), 0.9, message="The vapour root selection is not correct for this example.")
    _close(select_pr_root([0.1, 0.5, 0.9], phase="liquid"), 0.1, message="The liquid root selection is not correct for this example.")
    try:
        select_pr_root([0.1, 0.5, 0.9], phase="middle")
    except ValueError:
        pass
    else:
        _fail("Unknown phase labels should raise a ValueError.")


def _test_molar_volume_from_Z(ns: Mapping[str, Any]) -> None:
    molar_volume_from_Z = _callable(ns, "molar_volume_from_Z")
    Z = 1.0
    R = 8.314
    T = 300.0
    P = 1.0
    _close(molar_volume_from_Z(Z, R, T, P), R * T / P, message="The conversion from Z to molar volume is not correct.")
    v_liquid = molar_volume_from_Z(0.05, R, T, P)
    v_vapour = molar_volume_from_Z(0.95, R, T, P)
    _require(v_vapour > v_liquid, "A larger Z should correspond to a larger molar volume at fixed T and P.")


def _test_pure_fugacity_functions(ns: Mapping[str, Any]) -> None:
    pure_pr_ln_phi = _callable(ns, "pure_pr_ln_phi")
    pure_pr_phi = _callable(ns, "pure_pr_phi")
    Z = 1.0
    A = 1e-8
    B = 1e-8
    ln_phi = pure_pr_ln_phi(Z, A, B)
    phi = pure_pr_phi(Z, A, B)
    _require(abs(ln_phi) < 1e-5, "At the ideal-gas limit, ln(phi) should be close to 0.")
    _require(abs(phi - 1.0) < 1e-5, "At the ideal-gas limit, phi should be close to 1.")
    _close(pure_pr_phi(Z, A, B), np.exp(pure_pr_ln_phi(Z, A, B)), message="phi should equal exp(ln_phi).")


def _test_validate_mole_fractions(ns: Mapping[str, Any]) -> None:
    validate_mole_fractions = _callable(ns, "validate_mole_fractions")
    valid_x = pd.Series({"Methane": 0.4, "Propane": 0.6})
    _require(validate_mole_fractions(valid_x) is True, "Valid mole fractions should return True.")
    for invalid_x, reason in [
        (pd.Series({"Methane": 0.4, "Propane": 0.7}), "mole fractions that do not sum to one"),
        (pd.Series({"Methane": -0.1, "Propane": 1.1}), "negative mole fractions"),
    ]:
        try:
            validate_mole_fractions(invalid_x)
        except ValueError:
            pass
        else:
            _fail(f"The function accepted invalid input: {reason}.")


def _test_build_mixture_property_table(ns: Mapping[str, Any]) -> None:
    build_mixture_property_table = _callable(ns, "build_mixture_property_table")
    toy_df = pd.DataFrame({
        "Element": ["A", "B"],
        "Tc_K": [100.0, 200.0],
        "Pc_MPa": [1.0, 2.0],
        "w": [0.1, 0.2],
    })
    toy_x = pd.Series({"A": 0.25, "B": 0.75})
    result = build_mixture_property_table(["A", "B"], toy_x, toy_df)
    _require(isinstance(result, pd.DataFrame), "The mixture property table should be a pandas DataFrame.")
    _require(list(result.columns) == ["A", "B"], "The mixture property table should use component names as columns.")
    for row in ["mol_frac", "Tc_K", "Pc_MPa", "w"]:
        _require(row in result.index, f"The mixture property table is missing the `{row}` row.")
    _close(result.loc["mol_frac", "A"], 0.25, message="The mole fraction for component A is not in the expected location.")
    _close(result.loc["Tc_K", "B"], 200.0, message="The critical temperature for component B is not in the expected location.")


def _test_add_pr_pure_parameters(ns: Mapping[str, Any]) -> None:
    add_pr_pure_parameters = _callable(ns, "add_pr_pure_parameters")
    toy_table = pd.DataFrame({"A": [0.25, 100.0, 1.0, 0.1], "B": [0.75, 200.0, 2.0, 0.2]}, index=["mol_frac", "Tc_K", "Pc_MPa", "w"])
    result = add_pr_pure_parameters(toy_table, T=150.0, R=8.314)
    for row in ["kappa", "ai", "bi"]:
        _require(row in result.index, f"The output table is missing the `{row}` row.")
    _require((result.loc["ai"] > 0).all(), "All pure-component a_i values should be positive.")
    _require((result.loc["bi"] > 0).all(), "All pure-component b_i values should be positive.")
    _require("mol_frac" in toy_table.index, "The input table should not be damaged by the function.")


def _test_build_kij_matrix(ns: Mapping[str, Any]) -> None:
    build_kij_matrix = _callable(ns, "build_kij_matrix")
    toy_interactions = pd.DataFrame({"component_1": ["A"], "component_2": ["B"], "k12": [0.03]})
    kij = build_kij_matrix(["A", "B"], toy_interactions)
    _require(isinstance(kij, pd.DataFrame), "The k_ij result should be a pandas DataFrame.")
    _require(list(kij.index) == ["A", "B"], "The k_ij matrix index is not arranged as expected.")
    _require(list(kij.columns) == ["A", "B"], "The k_ij matrix columns are not arranged as expected.")
    _close(kij.loc["A", "A"], 0.0, message="The self-interaction entry for A is not represented correctly.")
    _close(kij.loc["B", "B"], 0.0, message="The self-interaction entry for B is not represented correctly.")
    _close(kij.loc["A", "B"], 0.03, message="The A-B interaction value was not placed correctly.")
    _close(kij.loc["B", "A"], 0.03, message="The B-A interaction value was not mirrored correctly.")
    _require(np.allclose(kij.values, kij.values.T), "The k_ij matrix should be symmetric.")


def _test_build_aij_matrix(ns: Mapping[str, Any]) -> None:
    build_aij_matrix = _callable(ns, "build_aij_matrix")
    components_toy = ["A", "B"]
    ai_values = pd.Series({"A": 0.25, "B": 1.10})
    kij = pd.DataFrame([[0.0, 0.01], [0.01, 0.0]], index=components_toy, columns=components_toy)
    aij = build_aij_matrix(components_toy, ai_values, kij)
    _require(isinstance(aij, pd.DataFrame), "The a_ij result should be a pandas DataFrame.")
    _require(aij.shape == (2, 2), "The a_ij matrix has the wrong shape.")
    _require(np.allclose(aij.values, aij.values.T), "The a_ij matrix should be symmetric.")
    _close(aij.loc["A", "A"], 0.25, message="The A-A diagonal entry is not correct.")
    _close(aij.loc["B", "B"], 1.10, message="The B-B diagonal entry is not correct.")
    expected_cross = np.sqrt(0.25 * 1.10) * (1 - 0.01)
    _close(aij.loc["A", "B"], expected_cross, message="The cross-interaction entry is not correct.")


def _test_calculate_mixture_a(ns: Mapping[str, Any]) -> None:
    calculate_mixture_a = _callable(ns, "calculate_mixture_a")
    x_toy = pd.Series({"A": 0.4, "B": 0.6})
    aij = pd.DataFrame([[0.25, 0.50], [0.50, 1.10]], index=["A", "B"], columns=["A", "B"])
    expected = 0.4 * 0.4 * 0.25 + 0.4 * 0.6 * 0.50 + 0.6 * 0.4 * 0.50 + 0.6 * 0.6 * 1.10
    _close(calculate_mixture_a(x_toy, aij), expected, message="The double-summation mixture a calculation is not correct.")


def _test_calculate_mixture_b(ns: Mapping[str, Any]) -> None:
    calculate_mixture_b = _callable(ns, "calculate_mixture_b")
    x_toy = pd.Series({"A": 0.4, "B": 0.6})
    bi = pd.Series({"A": 0.03, "B": 0.08})
    expected = 0.4 * 0.03 + 0.6 * 0.08
    _close(calculate_mixture_b(x_toy, bi), expected, message="The linear mixture b calculation is not correct.")


def _test_pressure_unit_conversions(ns: Mapping[str, Any]) -> None:
    convert_mpa_to_bar = _callable(ns, "convert_mpa_to_bar")
    convert_bar_to_mpa = _callable(ns, "convert_bar_to_mpa")
    _close(convert_mpa_to_bar(0.1), 1.0, message="0.1 MPa should be 1 bar.")
    _close(convert_mpa_to_bar(1.0), 10.0, message="1 MPa should be 10 bar.")
    _close(convert_bar_to_mpa(1.0), 0.1, message="1 bar should be 0.1 MPa.")
    _close(convert_bar_to_mpa(convert_mpa_to_bar(2.5)), 2.5, message="MPa -> bar -> MPa conversion should recover the original pressure.")


_PUBLIC_TESTS = {
    "ideal_gas_pressure": _test_ideal_gas_pressure,
    "ideal_gas_table": _test_ideal_gas_table,
    "calculate_kappa": _test_calculate_kappa,
    "calculate_alpha": _test_calculate_alpha,
    "calculate_ai_bi": _test_calculate_ai_bi,
    "peng_robinson_pressure": _test_peng_robinson_pressure,
    "calculate_A_B": _test_calculate_A_B,
    "pr_cubic_coefficients": _test_pr_cubic_coefficients,
    "roots_and_selection": _test_roots_and_selection,
    "molar_volume_from_Z": _test_molar_volume_from_Z,
    "pure_fugacity_functions": _test_pure_fugacity_functions,
    "validate_mole_fractions": _test_validate_mole_fractions,
    "build_mixture_property_table": _test_build_mixture_property_table,
    "add_pr_pure_parameters": _test_add_pr_pure_parameters,
    "build_kij_matrix": _test_build_kij_matrix,
    "build_aij_matrix": _test_build_aij_matrix,
    "calculate_mixture_a": _test_calculate_mixture_a,
    "calculate_mixture_b": _test_calculate_mixture_b,
    "pressure_unit_conversions": _test_pressure_unit_conversions,
}


def list_public_tests() -> list[str]:
    """Return the names of public tests available to the notebook."""
    return sorted(_PUBLIC_TESTS)


def run_public_test(test_name: str, namespace: Mapping[str, Any]) -> None:
    """
    Run one public test by name.

    Parameters
    ----------
    test_name:
        Name of the public test, for example "build_kij_matrix".
    namespace:
        Usually pass `globals()` from the notebook so the test can access the
        student's functions and variables.
    """
    if test_name not in _PUBLIC_TESTS:
        available = ", ".join(list_public_tests())
        raise KeyError(f"Unknown public test `{test_name}`. Available tests: {available}")

    try:
        _PUBLIC_TESTS[test_name](namespace)
    except PublicTestFailure as exc:
        print(f"❌ Public test failed: {test_name}")
        print(f"Hint: {exc}")
        raise
    except Exception as exc:
        print(f"❌ Public test raised {type(exc).__name__}: {test_name}")
        print(f"Hint: {exc}")
        raise
    else:
        print(f"✅ Public test passed: {test_name}")
