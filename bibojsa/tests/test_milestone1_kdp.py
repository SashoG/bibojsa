"""Validation milestone 1: KDP asymmetric GVM (Mosley et al., PRL 100,
133601 (2008)).

Configuration: type-II collinear SPDC in KDP, pump near 415 nm (e-pol),
degenerate daughter photons at 830 nm (o- and e-pol), 5 mm crystal, pump
modeled as a Gaussian with 4 nm FWHM spectral bandwidth (as used in the
paper's Fig. 1 simulation). Claim to reproduce: the pump's (e, 415 nm)
group velocity matches the o-polarized daughter's (830 nm) group
velocity, giving theta_PM ~= 0 deg (asymmetric GVM) and a high heralded
purity for the unfiltered signal, with no spectral filtering required.
"""
import numpy as np
import pytest

from bibojsa.materials.kdp import KDP
from bibojsa.phasematch import ThreeWaveConfig, find_phase_matching_angle
from bibojsa.gvm import group_index_at_angle, theta_pm_deg, gvm_regime
from bibojsa.jsa import build_jsa, fwhm_omega_from_fwhm_wavelength
from bibojsa.sellmeier import lambda_um_to_omega, C_LIGHT

LAMBDA_PUMP = 0.415
LAMBDA_DEGENERATE = 0.83
CONFIG = ThreeWaveConfig(plane="xz", pol_pump="e", pol_signal="o", pol_idler="e")


@pytest.fixture(scope="module")
def phase_matching_angle():
    return find_phase_matching_angle(LAMBDA_PUMP, LAMBDA_DEGENERATE, LAMBDA_DEGENERATE,
                                      CONFIG, KDP["x"], KDP["y"], KDP["z"])


def test_phase_matching_angle_is_physical(phase_matching_angle):
    """Type-II KDP phase matching for 415/830/830 nm should land at a
    non-trivial angle strictly between 0 and 90 deg (not at a principal
    axis, and not the degenerate-optics-axis singularity)."""
    deg = np.degrees(phase_matching_angle)
    assert 30.0 < deg < 85.0


def test_pump_group_velocity_matches_ordinary_daughter(phase_matching_angle):
    """The defining GVM property (Mosley et al.): v_g(e, 415nm) = v_g(o, 830nm),
    i.e. k_p' = k_s' to within ~0.1% -- this is what makes the o-polarized
    daughter's marginal spectrum set purely by phase matching, uncorrelated
    with its twin, without any spectral filter."""
    theta = phase_matching_angle
    wp = lambda_um_to_omega(LAMBDA_PUMP)
    ws = lambda_um_to_omega(LAMBDA_DEGENERATE)
    wi = lambda_um_to_omega(LAMBDA_DEGENERATE)

    kp_prime, _, conv_p, _ = group_index_at_angle(wp, theta, "xz", "e", KDP["x"], KDP["y"], KDP["z"],
                                                    return_diagnostics=True)
    ks_prime, _, conv_s, _ = group_index_at_angle(ws, theta, "xz", "o", KDP["x"], KDP["y"], KDP["z"],
                                                    return_diagnostics=True)
    ki_prime, _, conv_i, _ = group_index_at_angle(wi, theta, "xz", "e", KDP["x"], KDP["y"], KDP["z"],
                                                    return_diagnostics=True)
    assert conv_p and conv_s and conv_i

    ng_pump, ng_signal, ng_idler = kp_prime * C_LIGHT, ks_prime * C_LIGHT, ki_prime * C_LIGHT
    rel_mismatch = abs(ng_pump - ng_signal) / ng_pump
    assert rel_mismatch < 1e-3, (
        f"n_g(pump,e,415nm)={ng_pump:.6f} vs n_g(o,830nm)={ng_signal:.6f}, "
        f"rel. mismatch {rel_mismatch:.2e} (expected < 1e-3)"
    )
    # and clearly different from the *other* (e-polarized) daughter, confirming
    # this is a genuine asymmetric match, not an accidental near-degeneracy of all three
    assert abs(ng_pump - ng_idler) / ng_pump > 1e-2


def test_theta_pm_near_zero_degrees(phase_matching_angle):
    """theta_PM = arctan[-(kp'-ks')/(kp'-ki')] should be close to 0 deg
    (asymmetric GVM locus), matching LPR2010's classification convention."""
    theta = phase_matching_angle
    wp = lambda_um_to_omega(LAMBDA_PUMP)
    ws = lambda_um_to_omega(LAMBDA_DEGENERATE)
    wi = lambda_um_to_omega(LAMBDA_DEGENERATE)
    kp_prime = group_index_at_angle(wp, theta, "xz", "e", KDP["x"], KDP["y"], KDP["z"])
    ks_prime = group_index_at_angle(ws, theta, "xz", "o", KDP["x"], KDP["y"], KDP["z"])
    ki_prime = group_index_at_angle(wi, theta, "xz", "e", KDP["x"], KDP["y"], KDP["z"])

    tpm = theta_pm_deg(kp_prime, ks_prime, ki_prime)
    assert abs(tpm) < 1.0
    assert gvm_regime(kp_prime, ks_prime, ki_prime) == "asymmetric_signal"


def test_unfiltered_purity_high_and_grid_converged(phase_matching_angle):
    """JSA for a 5 mm KDP crystal, 4 nm FWHM (intensity) pump bandwidth at
    415 nm: unfiltered heralded purity should be high (no filtering
    needed, per the paper's central claim) and stable under grid
    refinement at fixed window (grid-size convergence), per the brief's
    required convergence check. Window-size sensitivity is documented
    separately (see module docstring / milestone report) since the
    daughter photon is intentionally broadband in this GVM configuration.
    """
    theta = phase_matching_angle
    L = 5e-3  # 5 mm, per Mosley et al.
    pump_fwhm_omega = fwhm_omega_from_fwhm_wavelength(0.004, LAMBDA_PUMP)  # 4 nm FWHM

    wp0 = lambda_um_to_omega(LAMBDA_PUMP)
    ws0 = lambda_um_to_omega(LAMBDA_DEGENERATE)
    wi0 = lambda_um_to_omega(LAMBDA_DEGENERATE)

    window = 4e14  # rad/s; captures the broadband GVM-matched daughter's JSA support
    purities = []
    for n in (120, 180, 240):
        ws_axis = np.linspace(ws0 - window, ws0 + window, n)
        wi_axis = np.linspace(wi0 - window, wi0 + window, n)
        result = build_jsa(ws_axis, wi_axis, wp0, theta, CONFIG, KDP["x"], KDP["y"], KDP["z"],
                            L, pump_fwhm_omega, pump_shape="gaussian", include_phase=True)
        purities.append(result.purity)

    # grid-size convergence: purity should barely change with finer grids
    assert max(purities) - min(purities) < 1e-2, (
        f"purity not grid-converged across N=100,160,220: {purities}"
    )
    # high purity without any heralding-arm filter
    assert purities[-1] > 0.9, f"unfiltered purity {purities[-1]:.3f}, expected > 0.9"
