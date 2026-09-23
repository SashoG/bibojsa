"""Validation milestone 3: BIBO sanity check vs. LPR2010 Fig. 10b.

Configuration: negative type-I (ooe) collinear SPDC phase matching in the
x-z principal plane of BIBO (idler=o, signal=o, pump=e, per this project's
lambda1 >= lambda2 > lambda3 = idler,signal,pump convention), pump fixed
at 800 nm, scanning the phase-matched crystal-cut angle theta(lambda_s)
and the three pairwise group-delay (GVM) differences among pump/signal/
idler as the signal wavelength is scanned from near-degenerate (1.6 um)
down toward 1.1 um.

Target claim (LPR2010 -- V. Petrov et al., Laser & Photon. Rev. 4, 53
(2010), Fig. 10b and its caption): "Moving with the signal wavelength
towards degeneracy, one can see that BIBO always exhibits a smaller GVM
parameter Delta21 [than BBO], equivalent to a larger spectral acceptance.
There is also a fairly broad spectral range, only in the case of BIBO,
near lambda2 = 1340 nm, where the signal and idler pulses travel with
group velocities very close to that of the pump."

This is a qualitative sanity check (Fig. 10b is a plot, not a numeric
table -- LPR2010 gives no tabulated GVM values for this exact scan), so
the tests below check the qualitative shape of the claim rather than
matching literature numbers to a tolerance:

  1. The phase-matched cut angle stays near the ~11 deg region LPR2010's
     Fig. 10a text describes ("maximum polar angle slightly above 11
     deg") for negative type-I in BIBO's x-z plane.
  2. Delta21 = GD_s - GD_i shrinks in magnitude monotonically as the
     signal wavelength approaches degeneracy (trivial at exact
     degeneracy where lambda_s = lambda_i, but LPR2010's point is that it
     is already much smaller than BBO's well before degeneracy -- BBO is
     out of scope here, not implemented per the brief, so only BIBO's own
     internal trend is checked).
  3. The *three* pump-signal-idler group delays (not just Delta21) draw
     close together -- i.e., Delta31 = GD_p - GD_i and Delta32 = GD_p -
     GD_s both cross zero -- in a band bracketing lambda_s = 1.34 um,
     reproducing "signal and idler travel with group velocities very
     close to that of the pump" near lambda2 = 1340 nm. This is BIBO's
     first EPM (extended-phase-matching) candidate per the brief.

No deff/nonlinear-coefficient calculation is performed for this
milestone (not required for the GVM-only sanity check requested); the
Table 1 deff expressions and d_il coefficients remain unimplemented
pending the brief's explicit go-ahead, per "ask before adding any
physics not specified here."
"""
import numpy as np
import pytest

from bibojsa.materials.bibo import BIBO
from bibojsa.phasematch import ThreeWaveConfig, find_phase_matching_angle
from bibojsa.gvm import phase_matched_group_delays_ps_per_mm

PLANE = "xz"
LAMBDA_PUMP = 0.800  # um
CONFIG = ThreeWaveConfig(plane=PLANE, pol_pump="e", pol_signal="o", pol_idler="o")
THETA_BOUNDS = (np.radians(0.05), np.radians(24.9))  # stay inside Omega(800nm) ~= 25.9 deg


def _idler_um(lambda_signal_um: float) -> float:
    return 1.0 / (1.0 / LAMBDA_PUMP - 1.0 / lambda_signal_um)


def _scan_point(lambda_signal_um: float):
    lambda_idler_um = _idler_um(lambda_signal_um)
    theta = find_phase_matching_angle(
        LAMBDA_PUMP, lambda_signal_um, lambda_idler_um, CONFIG,
        BIBO["x"], BIBO["y"], BIBO["z"], theta_bounds=THETA_BOUNDS, n_scan=400,
    )
    GD_p, GD_s, GD_i = phase_matched_group_delays_ps_per_mm(
        LAMBDA_PUMP, lambda_signal_um, lambda_idler_um, theta, PLANE,
        "e", "o", "o", BIBO["x"], BIBO["y"], BIBO["z"],
    )
    return theta, lambda_idler_um, GD_p, GD_s, GD_i


SIGNAL_SCAN_UM = [1.10, 1.15, 1.20, 1.25, 1.30, 1.34, 1.40, 1.45, 1.50, 1.55, 1.59]


def test_phase_matched_angle_stays_near_11_degrees():
    """LPR2010's Fig. 10a text: negative type-I phase matching in BIBO's
    x-z plane near 800 nm pump reaches "a maximum polar angle slightly
    above 11 deg" before retracing. The angle should stay in a narrow
    band around 10-11.5 deg across the whole scan (not, e.g., drift to
    the eeo/other-octant regime or to the y-z / x-y planes)."""
    for lam_s in SIGNAL_SCAN_UM:
        theta, _, _, _, _ = _scan_point(lam_s)
        theta_deg = np.degrees(theta)
        assert 9.5 < theta_deg < 11.5, f"lambda_s={lam_s}: theta={theta_deg:.3f} deg out of expected band"


def test_delta21_shrinks_toward_degeneracy():
    """Delta21 = GD_s - GD_i should shrink in magnitude monotonically
    (up to small numerical wiggle) as lambda_s -> 1.6 um degenerate point,
    consistent with BIBO's comparatively small/vanishing Delta21 (BBO
    comparison out of scope; only the internal BIBO trend is checked)."""
    abs_d21 = []
    for lam_s in SIGNAL_SCAN_UM:
        _, _, _, GD_s, GD_i = _scan_point(lam_s)
        abs_d21.append(abs(GD_s - GD_i) * 1000.0)  # fs/mm
    # allow small non-monotonic wiggle near the very end (theta plateaus),
    # but require a clear, large overall decrease from short to long lambda_s
    assert abs_d21[0] > 40.0, f"expected large Delta21 near lambda_s=1.10um, got {abs_d21[0]:.1f} fs/mm"
    assert abs_d21[-1] < 2.0, f"expected near-zero Delta21 near degeneracy, got {abs_d21[-1]:.1f} fs/mm"
    # overall decreasing trend (each successive point no more than a small amount above the previous)
    for a, b in zip(abs_d21, abs_d21[1:]):
        assert b <= a + 1e-6, "Delta21 magnitude should not increase moving toward degeneracy"


def test_all_three_group_delays_converge_near_1340nm():
    """The central Fig. 10b claim: near lambda2 = 1340 nm, signal and
    idler group velocities are close to the pump's -- i.e. Delta31 = GD_p
    - GD_i and Delta32 = GD_p - GD_s both cross zero in a band bracketing
    1.34 um, so that max(|Delta31|,|Delta21|,|Delta32|) there is much
    smaller than at lambda_s=1.10um (far from the reported EPM-like
    point)."""
    _, _, GD_p_near, GD_s_near, GD_i_near = _scan_point(1.34)
    d21_near = (GD_s_near - GD_i_near) * 1000.0
    d31_near = (GD_p_near - GD_i_near) * 1000.0
    d32_near = (GD_p_near - GD_s_near) * 1000.0
    max_near = max(abs(d21_near), abs(d31_near), abs(d32_near))

    _, _, GD_p_far, GD_s_far, GD_i_far = _scan_point(1.10)
    d21_far = (GD_s_far - GD_i_far) * 1000.0
    d31_far = (GD_p_far - GD_i_far) * 1000.0
    d32_far = (GD_p_far - GD_s_far) * 1000.0
    max_far = max(abs(d21_far), abs(d31_far), abs(d32_far))

    assert max_near < 10.0, f"expected all GVM parameters small (<10 fs/mm) near 1340nm, got {max_near:.1f}"
    assert max_near < 0.2 * max_far, "expected a large collapse in GVM spread relative to lambda_s=1.10um"


def test_delta31_and_delta32_bracket_1340nm():
    """Delta31 (GD_p - GD_i) and Delta32 (GD_p - GD_s) should each change
    sign somewhere in [1.25, 1.40] um, bracketing LPR2010's reported
    lambda2 = 1340 nm -- i.e. this is a genuine simultaneous zero-crossing
    region, not just a coincidentally small value at one sampled point."""
    lo, hi = 1.25, 1.40
    _, _, GD_p_lo, GD_s_lo, GD_i_lo = _scan_point(lo)
    _, _, GD_p_hi, GD_s_hi, GD_i_hi = _scan_point(hi)
    d31_lo, d31_hi = GD_p_lo - GD_i_lo, GD_p_hi - GD_i_hi
    d32_lo, d32_hi = GD_p_lo - GD_s_lo, GD_p_hi - GD_s_hi
    assert d31_lo * d31_hi < 0, "expected Delta31 to change sign between 1.25 and 1.40 um"
    assert d32_lo * d32_hi < 0, "expected Delta32 to change sign between 1.25 and 1.40 um"
