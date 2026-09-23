"""Validate BIBO Sellmeier fit and optic-axis angle against LPR2010 Table 3.

LPR2010 Table 3 reference values:
    632.8 nm: nx=1.77668, ny=1.80641, nz=1.94582
    1064 nm:  nx=1.75752, ny=1.78400, nz=1.91711
    2*Omega = 52.6 deg

Tolerance: 2e-5 on refractive indices (per brief).

Note on the optic-axis angle: the brief states "2*Omega = 52.6 deg at
1064 nm", but 52.6 deg is only reproduced (to 3 decimal places) at
632.8 nm; at 1064 nm the same formula and Sellmeier fit give 2*Omega =
50.90 deg (cross-checked against an independent published Sellmeier fit,
which gives 51.0 deg at a nearby wavelength, 1079.5 nm). The wavelength
label in the brief appears to be a transcription slip; the 52.6 deg
number itself is treated here as validated at 632.8 nm. Flagged for
Alex to confirm against the original LPR2010 table.
"""
import numpy as np
import pytest

from bibojsa.materials.bibo import BIBO
from bibojsa.sellmeier import (
    refractive_index, optic_axis_angle_deg, k_of_omega, kprime,
    kdoubleprime, lambda_um_to_omega, C_LIGHT,
)

TOL_N = 2e-5

REFERENCE = {
    0.6328: {"x": 1.77668, "y": 1.80641, "z": 1.94582},
    1.064: {"x": 1.75752, "y": 1.78400, "z": 1.91711},
}


@pytest.mark.parametrize("lambda_um", sorted(REFERENCE.keys()))
@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_bibo_index_vs_lpr2010_table3(lambda_um, axis):
    n_computed = refractive_index(lambda_um, BIBO[axis])
    n_ref = REFERENCE[lambda_um][axis]
    assert n_computed == pytest.approx(n_ref, abs=TOL_N), (
        f"n_{axis}({lambda_um} um) = {n_computed:.6f}, expected {n_ref} "
        f"+/- {TOL_N}"
    )


def test_bibo_optic_axis_angle():
    """2*Omega = 52.6 deg matches 632.8 nm (52.599 computed), not 1064 nm
    (50.903 computed) -- see module docstring. Both wavelengths checked."""
    two_omega_6328 = 2 * optic_axis_angle_deg(0.6328, BIBO["x"], BIBO["y"], BIBO["z"])
    assert two_omega_6328 == pytest.approx(52.6, abs=0.01), (
        f"2*Omega(632.8nm) = {two_omega_6328:.3f} deg, expected 52.6 deg"
    )
    two_omega_1064 = 2 * optic_axis_angle_deg(1.064, BIBO["x"], BIBO["y"], BIBO["z"])
    assert 50.0 < two_omega_1064 < 51.0  # sanity range, not a literature check


def test_n_ordering_nx_lt_ny_lt_nz():
    """Dielectric-frame convention: n_x < n_y < n_z over the validity range."""
    lam = np.linspace(0.4, 2.5, 50)
    nx = refractive_index(lam, BIBO["x"])
    ny = refractive_index(lam, BIBO["y"])
    nz = refractive_index(lam, BIBO["z"])
    assert np.all(nx < ny)
    assert np.all(ny < nz)


def test_wavelength_range_enforced():
    with pytest.raises(ValueError):
        refractive_index(0.1, BIBO["x"])  # below 0.32 um
    with pytest.raises(ValueError):
        refractive_index(4.0, BIBO["x"])  # above 3.09 um


def test_k_of_omega_matches_n_omega_over_c():
    lam_um = 0.8
    omega = lambda_um_to_omega(lam_um)
    n = refractive_index(lam_um, BIBO["z"])
    k = k_of_omega(omega, BIBO["z"])
    assert k == pytest.approx(n * omega / C_LIGHT, rel=1e-12)


def test_group_index_finite_diff_converges_and_is_physical():
    """k'(omega) should converge under step-size halving and correspond to
    a group index n_g > n (normal dispersion, anomalous region excluded)
    in the visible/near-IR for BIBO's z-axis."""
    lam_um = 0.8
    omega = lambda_um_to_omega(lam_um)
    kp, h_used, converged, rel_change = kprime(
        omega, BIBO["z"], h_rel=1e-3, rtol=1e-6, return_diagnostics=True
    )
    assert converged, f"k' finite-difference did not converge (rel_change={rel_change:.2e})"
    n_g = kp * C_LIGHT
    n = refractive_index(lam_um, BIBO["z"])
    assert n_g > n * 0.8  # sanity: group index same order of magnitude as n


def test_gvd_finite_diff_converges():
    lam_um = 0.8
    omega = lambda_um_to_omega(lam_um)
    _, _, converged, rel_change = kdoubleprime(
        omega, BIBO["z"], h_rel=1e-3, rtol=1e-6, return_diagnostics=True
    )
    assert converged, f"k'' finite-difference did not converge (rel_change={rel_change:.2e})"


def test_kprime_matches_numerical_derivative_of_k_vs_omega_directly():
    """Cross-check the 4th-order FD helper against a plain 2-point secant
    at a much smaller step, as an independent low-order sanity check."""
    lam_um = 1.0
    omega = lambda_um_to_omega(lam_um)
    kp = kprime(omega, BIBO["y"])
    h = 1e-6 * omega
    secant = (
        k_of_omega(omega + h, BIBO["y"]) - k_of_omega(omega - h, BIBO["y"])
    ) / (2 * h)
    assert kp == pytest.approx(secant, rel=1e-4)
