"""Milestone 5: JSA / heralded purity at BIBO's 1.34 um triple-GVM point.

Configuration: x-z plane, negative type-I ooe (idler o, signal o, pump e),
pump 800 nm, theta_cut = 10.9577 deg. The exact theta_PM = 45 deg
(symmetric-GVM) signal wavelength solves k_p' = (k_s' + k_i')/2 at
lambda_s = 1.3409 um (lambda_i = 1.9831 um) -- the operating point
identified in milestone 3 from LPR2010 Fig. 10b's "lambda2 ~ 1340 nm"
region.

Physics documented by these tests (the milestone's central finding):

Although theta_PM crosses 45 deg here, this is a *triple* GVM point --
all three pairwise group-delay differences are simultaneously tiny
(|Delta31| = |Delta32| ~ 2.4 fs/mm at the crossing, vs ~147 fs/mm for
the canonical KTP EPM configuration of Konig & Wong / Laudenbach). The
sinc phase-matching function is therefore orders of magnitude too
broadband for its width ever to match a realistic pump bandwidth at
0.1-10 mm crystal lengths: the JSA is a long frequency-anticorrelated
stripe whose extent is set by second-order dispersion (and ultimately
runs into the Sellmeier validity window), NOT the near-separable
round JSA of a true EPM source. Unfiltered heralded purity stays below
~0.2 across the whole (L, bandwidth) region scanned (best P ~ 0.175 at
L = 2 mm, 40 nm pump FWHM, grid-converged to <1e-3), i.e. this point is
a huge-spectral-acceptance / ultra-broadband-biphoton operating point,
not a high-purity heralded single-photon source. Positive identification
of high-purity BIBO operating points is deferred to the exploration
survey (asymmetric GVM loci and symmetric loci with sizable |Delta|).

The tests assert (a) the 45 deg crossing exists where claimed, (b) the
triple-GVM scale is fs/mm (not tens/hundreds), and (c) the low-purity
character -- so that any future dispersion-data change that silently
moved or "fixed" this conclusion would be caught.
"""
import numpy as np
import pytest
from scipy.optimize import brentq

from bibojsa.materials.bibo import BIBO
from bibojsa.phasematch import ThreeWaveConfig, find_phase_matching_angle
from bibojsa.gvm import group_index_at_angle, theta_pm_deg
from bibojsa.jsa import build_jsa, fwhm_omega_from_fwhm_wavelength
from bibojsa.sellmeier import lambda_um_to_omega

X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]
PLANE = "xz"
LAMBDA_PUMP = 0.800
CONFIG = ThreeWaveConfig(plane=PLANE, pol_pump="e", pol_signal="o", pol_idler="o")
THETA_BOUNDS = (np.radians(0.05), np.radians(24.9))

# Design point (examples/locate_epm.py)
LAMBDA_S_STAR = 1.3409
LAMBDA_I_STAR = 1.9831
THETA_STAR = np.radians(10.9577)


def _idler_um(lam_s):
    return 1.0 / (1.0 / LAMBDA_PUMP - 1.0 / lam_s)


def _kprimes(lam_s):
    lam_i = _idler_um(lam_s)
    theta = find_phase_matching_angle(LAMBDA_PUMP, lam_s, lam_i, CONFIG,
                                       X, Y, Z, theta_bounds=THETA_BOUNDS, n_scan=400)
    wp, ws, wi = (lambda_um_to_omega(l) for l in (LAMBDA_PUMP, lam_s, lam_i))
    kp = group_index_at_angle(wp, theta, PLANE, "e", X, Y, Z)
    ks = group_index_at_angle(ws, theta, PLANE, "o", X, Y, Z)
    ki = group_index_at_angle(wi, theta, PLANE, "o", X, Y, Z)
    return theta, kp, ks, ki


def test_symmetric_gvm_crossing_at_1341nm():
    """k_p' - (k_s'+k_i')/2 changes sign in [1.30, 1.40] um and the root
    is at 1.341 um (+- 2 nm), i.e. theta_PM crosses 45 deg (mod 180) there."""
    def mismatch(lam_s):
        _, kp, ks, ki = _kprimes(lam_s)
        return kp - 0.5 * (ks + ki)

    assert mismatch(1.30) * mismatch(1.40) < 0
    lam_star = brentq(mismatch, 1.30, 1.40, xtol=1e-5)
    assert lam_star == pytest.approx(LAMBDA_S_STAR, abs=0.002)

    _, kp, ks, ki = _kprimes(lam_star)
    tpm = theta_pm_deg(kp, ks, ki) % 180.0
    assert tpm == pytest.approx(45.0, abs=0.5)


def test_triple_gvm_scale_is_few_fs_per_mm():
    """At the 45 deg crossing the *individual* mismatches are only a few
    fs/mm (|Delta31| = |Delta32| ~ 2.4 fs/mm) -- roughly 60x smaller than
    KTP's EPM design (~147 fs/mm). This is what makes the sinc too
    broadband to decorrelate the JSA (triple GVM, not classic EPM)."""
    _, kp, ks, ki = _kprimes(LAMBDA_S_STAR)
    c = 2.99792458e8
    d32_fs_mm = (kp - ks) * c / 0.299792458 * 1000.0  # GD_p - GD_s in fs/mm
    d31_fs_mm = (kp - ki) * c / 0.299792458 * 1000.0
    assert abs(d32_fs_mm) < 10.0
    assert abs(d31_fs_mm) < 10.0
    assert abs(d32_fs_mm + d31_fs_mm) < 0.5  # symmetric: Delta32 ~= -Delta31


def test_unfiltered_purity_low_at_practical_length():
    """Central (negative) milestone finding: at L = 2 mm with a broad
    40 nm FWHM pump -- the best point of the scanned (L, bandwidth) grid
    -- the unfiltered heralded purity is only ~0.175 (grid-converged),
    far below the ~0.8 of a true EPM source. Guard both sides so a
    dispersion-data regression that changes this materially is caught."""
    L = 2e-3
    fw = fwhm_omega_from_fwhm_wavelength(0.040, LAMBDA_PUMP)
    ws0 = lambda_um_to_omega(LAMBDA_S_STAR)
    wi0 = lambda_um_to_omega(LAMBDA_I_STAR)
    wp0 = lambda_um_to_omega(LAMBDA_PUMP)
    # window matching the converged adaptive window of the milestone scan:
    # +-4.4e14 rad/s about each center, clipped to Sellmeier validity
    w_min = lambda_um_to_omega(3.05)
    w_max = lambda_um_to_omega(0.33)
    half = 4.4e14
    ws_axis = np.linspace(max(ws0 - half, w_min), min(ws0 + half, w_max), 160)
    wi_axis = np.linspace(max(wi0 - half, w_min), min(wi0 + half, w_max), 160)
    res = build_jsa(ws_axis, wi_axis, wp0, THETA_STAR, CONFIG, X, Y, Z,
                    L, fw, pump_shape="gaussian", include_phase=True)
    assert res.purity < 0.30, "purity unexpectedly high -- dispersion data changed?"
    assert res.purity > 0.05, "purity unexpectedly low vs the recorded ~0.17"
    assert res.schmidt_number > 3.0
