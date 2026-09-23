"""Regression tests for the continuous-pump landscape module.

The landscape sweeps use vectorised re-expressions of the scalar
primitives in `phasematch`/`biaxial`/`gvm`; these tests assert that the two
agree, so that a result quoted from the landscape and a result quoted from
a point-by-point calculation are the same computation.

They also pin the structural results the manuscript relies on:
  * the type-I collapse at frequency degeneracy (k_s' = k_i', so
    theta_PM = -45 deg at every pump wavelength, and a triple-GVM point
    where all three GVM functions vanish together);
  * the three degenerate GVM points of the x-z type-II interaction;
  * the negative results for the x-y and y-z planes.
"""
import numpy as np
import pytest

from bibojsa import landscape as L
from bibojsa.gvm import group_index_at_angle
from bibojsa.materials.bibo import BIBO
from bibojsa.phasematch import ThreeWaveConfig, find_phase_matching_angle
from bibojsa.sellmeier import lambda_um_to_omega

X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]
OEO = ThreeWaveConfig(plane="xz", pol_pump="o", pol_signal="e", pol_idler="o")


# ----------------------------------------------------------- primitives

@pytest.mark.parametrize("plane,pol", [("xz", "o"), ("xz", "e"),
                                       ("yz", "e"), ("xy", "e")])
@pytest.mark.parametrize("theta_deg", [5.0, 20.0, 39.422, 60.0, 85.0])
def test_kprime_vec_matches_step_halving(plane, pol, theta_deg):
    """The fixed-step vectorised k' reproduces the step-halving scalar one."""
    theta = np.radians(theta_deg)
    lams = np.array([0.34, 0.5, 0.68, 0.81, 1.14, 1.74, 2.23, 2.81, 3.02])
    vec = L.kprime_vec(lambda_um_to_omega(lams), theta, plane, pol, X, Y, Z)
    for lam, v in zip(lams, vec):
        ref = group_index_at_angle(lambda_um_to_omega(lam), theta, plane, pol,
                                   X, Y, Z)
        assert abs(v - ref) / abs(ref) < 1e-9


@pytest.mark.parametrize("lam_p,lam_s", [(0.810, 1.1388), (0.800, 1.1090),
                                         (0.8696, 1.73919), (0.775, 1.0413)])
def test_solve_theta_matches_brentq(lam_p, lam_s):
    """The vectorised bisection reproduces the scalar brentq root."""
    lam_i = 1.0 / (1.0 / lam_p - 1.0 / lam_s)
    th_vec, li = L.solve_theta_grid(lam_p, lam_s, "xz", "oeo", X, Y, Z)
    th_ref = find_phase_matching_angle(lam_p, lam_s, lam_i, OEO, X, Y, Z)
    assert abs(np.degrees(th_vec) - np.degrees(th_ref)) < 1e-8
    assert abs(li - lam_i) < 1e-9


def test_solve_theta_grid_broadcasts():
    """Shapes are preserved, and nodes whose idler leaves the Sellmeier
    window come back as NaN rather than as a spurious root."""
    ls = np.linspace(1.10, 1.20, 7)
    th, li = L.solve_theta_grid(np.full_like(ls, 0.81), ls, "xz", "oeo",
                                X, Y, Z)
    assert th.shape == ls.shape and li.shape == ls.shape
    assert np.isfinite(th[li <= L.LAM_MAX]).all()
    assert not np.isfinite(th[li > L.LAM_MAX]).any()


# ------------------------------------------- type-I collapse at degeneracy

@pytest.mark.parametrize("seq", ["ooe", "eeo"])
@pytest.mark.parametrize("lam_p", [0.70, 0.78, 0.86, 0.94])
def test_type_i_pinned_to_anticorrelated_at_degeneracy(seq, lam_p):
    """Identically polarised daughters have identical group slowness, so
    theta_PM is -45 deg (mod 180) at every pump wavelength."""
    lam_d = 2.0 * lam_p
    th, _ = L.solve_theta_grid(lam_p, lam_d * 0.9999999, "xz", seq, X, Y, Z)
    if not np.isfinite(th):
        pytest.skip("not phase-matchable at this pump")
    g = L.gvm_functions(lam_p, lam_d, lam_d, th, "xz", seq, X, Y, Z)
    assert abs(g["ks"] - g["ki"]) / abs(g["ks"]) < 1e-12
    wrapped = (float(g["theta_pm"]) + 90.0) % 180.0 - 90.0
    assert abs(wrapped + 45.0) < 1e-6


@pytest.mark.parametrize("seq,lam_p_expected", [("ooe", 0.81870),
                                                ("eeo", 0.89676)])
def test_type_i_triple_gvm_point(seq, lam_p_expected):
    """Where the type-I curve flips sign, all three GVM functions vanish."""
    p = L.degenerate_gvm_point("xz", seq, 0, X, Y, Z, (0.60, 1.20))
    assert p is not None
    assert abs(p.lam_p_um - lam_p_expected) < 5e-4
    g = L.gvm_functions(p.lam_p_um, p.lam_s_um, p.lam_i_um,
                        np.radians(p.theta_deg), "xz", seq, X, Y, Z)
    for key in ("g0", "g45", "g90"):
        assert abs(float(g[key])) * 1e12 < 0.01   # fs/mm


# ----------------------------------------- type-II degenerate GVM points

@pytest.mark.parametrize("target,lam_p_expected,theta_expected", [
    (0,  0.682241, 45.968),
    (45, 0.869583, 43.770),
    (90, 1.116841, 46.110),
])
def test_type_ii_degenerate_points(target, lam_p_expected, theta_expected):
    """The three degenerate GVM points of the x-z oeo interaction.

    Compare with Jin et al., Quantum Eng. 2, e38 (2020), who report
    687/875/1119 nm at 46.0/43.8/46.0 deg from a different Sellmeier fit.
    """
    p = L.degenerate_gvm_point("xz", "oeo", target, X, Y, Z, (0.50, 1.40))
    assert p is not None
    assert abs(p.lam_p_um - lam_p_expected) < 1e-3
    assert abs(p.theta_deg - theta_expected) < 0.01
    assert abs(p.theta_pm_deg - target) < 1e-3
    assert p.degeneracy < 1e-6


# -------------------------------------------------------- negative results

def test_xy_ooe_never_phase_matches():
    """No collinear SPDC solution anywhere on the x-y ooe surface."""
    for lam_p in np.linspace(0.34, 1.50, 25):
        ls = np.linspace(lam_p * 1.001, 2.0 * lam_p * 0.999, 40)
        th, _ = L.solve_theta_grid(np.full_like(ls, lam_p), ls, "xy", "ooe",
                                   X, Y, Z, n_theta=601)
        assert not np.isfinite(th).any()


def test_xy_eoe_phase_matches_but_meets_no_gvm_condition():
    """x-y eoe phase-matches only for deep-blue pumps, and each GVM
    function keeps one sign over the whole admissible surface."""
    signs = {"g0": set(), "g45": set(), "g90": set()}
    pumps = []
    for lam_p in np.arange(0.34, 1.501, 0.01):
        ls = np.linspace(lam_p * 1.0005, 2.0 * lam_p * 0.999999, 120)
        th, li = L.solve_theta_grid(np.full_like(ls, lam_p), ls, "xy", "eoe",
                                    X, Y, Z, n_theta=1201)
        ok = np.isfinite(th)
        if not ok.any():
            continue
        pumps.append(lam_p)
        g = L.gvm_functions(np.full(ok.sum(), lam_p), ls[ok], li[ok], th[ok],
                            "xy", "eoe", X, Y, Z)
        for key in signs:
            signs[key] |= set(np.sign(np.asarray(g[key])).tolist())
    assert pumps, "expected x-y eoe to phase-match somewhere"
    assert 0.339 < min(pumps) < 0.345 and 0.385 < max(pumps) < 0.395
    for key, sg in signs.items():
        assert len(sg) == 1, f"{key} changes sign; a GVM crossing exists"


def test_yz_eeo_meets_no_gvm_condition():
    """y-z eeo is phase-matchable everywhere but its ridge orientation
    never reaches 0, 45 or 90 degrees."""
    lo, hi = 1e9, -1e9
    for lam_p in np.arange(0.34, 1.501, 0.02):
        ls = np.linspace(lam_p * 1.0005, 2.0 * lam_p * 0.999999, 120)
        th, li = L.solve_theta_grid(np.full_like(ls, lam_p), ls, "yz", "eeo",
                                    X, Y, Z, n_theta=1201)
        ok = np.isfinite(th)
        if not ok.any():
            continue
        g = L.gvm_functions(np.full(ok.sum(), lam_p), ls[ok], li[ok], th[ok],
                            "yz", "eeo", X, Y, Z)
        t = (np.asarray(g["theta_pm"]) + 90.0) % 180.0 - 90.0
        lo, hi = min(lo, t.min()), max(hi, t.max())
    assert -45.01 < lo and hi < -5.0


def test_yz_oeo_has_zero_deff():
    from bibojsa.deff import deff_bibo_pm_per_v
    for th_deg in np.linspace(5.0, 175.0, 25):
        d = deff_bibo_pm_per_v("yz", "o", "e", "o", np.radians(th_deg),
                               lambda_pump_um=0.8, enforce_range=False)
        assert abs(d) < 1e-12
