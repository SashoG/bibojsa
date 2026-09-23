"""Validation milestone 2: KTP symmetric GVM / extended phase matching (EPM).

Configuration: type-II collinear SPDC in KTP, propagation along the
crystallographic x-axis (walk-off-free principal axis, xz principal
plane at theta=90 deg), pump near 791 nm (o-pol) -> degenerate daughters
at 1582 nm (o- and e-pol). Target claim (Konig & Wong, Appl. Phys. Lett.
84, 1644 (2004) [papers/09.]; Evans et al., Phys. Rev. Lett. 105, 253601
(2010); Graffitti et al., Phys. Rev. A 98, 053811 (2018)): the pump's
group velocity matches the *average* of the two daughters' group
velocities (symmetric GVM / EPM, dispersion parameter D=1, theta_PM=45
deg), giving unfiltered purity ~ 0.80-0.85.

STATUS: fully validated. Precise literature numbers for this exact
configuration -- group delays GD_p=6.027, GD_s=5.880, GD_i=6.175 ps/mm,
D=0.993, QPM period 46.1-46.2 um -- are taken from Laudenbach, Jin,
Greganti, Hentschel, Walther & Hubel, "Numerical Investigation of Photon-
Pair Generation in Periodically Poled MTiOXO4", Phys. Rev. Applied 8,
024035 (2017) [arXiv:1706.04559], Fig. 1(a)/(b)/(c) caption.

HISTORY (see PROGRESS.md for the full account): this milestone was
originally closed as a "partial validation" because Konig & Wong's own
n_y Sellmeier fit was believed unsourceable (paywalled, not reproduced
verbatim in any secondary reference found by an extensive search), and a
Kato & Takaoka (2002) n_y was used as the closest available proxy. That
proxy reproduced n_y's *absolute* value well but not its *dispersion
slope* near 791 nm, which was enough to misclassify the GVM regime as
asymmetric (D~0) rather than the literature's symmetric/EPM (D~0.993).
Konig & Wong's paper (Appl. Phys. Lett. 84, 1644 (2004)) has since been
located directly (papers/09.) and its own Eq. (3) n_y coefficients
(A1=2.09930, A2=0.922683, A3=0.0467695, A4=0.0138408) are now used
verbatim in bibojsa.materials.ktp -- closing the gap. n_z was already
Fradkin et al.'s fit, which Konig & Wong themselves also adopted
unmodified, so no change was needed there.

All five headline numbers (GD_p, GD_s, GD_i, D, QPM period) now match
Laudenbach et al. 2017 to well under 1%.
"""
import numpy as np
import pytest

from bibojsa.materials.ktp import KTP
from bibojsa.phasematch import ThreeWaveConfig, find_qpm_period_um
from bibojsa.gvm import phase_matched_group_delays_ps_per_mm, theta_pm_deg

THETA = np.pi / 2  # crystallographic x-axis, xz principal plane
PLANE = "xz"
LAMBDA_PUMP = 0.791
LAMBDA_DEGENERATE = 1.582
CONFIG = ThreeWaveConfig(plane=PLANE, pol_pump="o", pol_signal="o", pol_idler="e")

# Literature reference values (Laudenbach et al. 2017, Fig. 1 caption / Fig. 3)
LIT_GD_P = 6.027  # ps/mm
LIT_GD_S = 5.880  # ps/mm
LIT_GD_I = 6.175  # ps/mm
LIT_D = 0.993
LIT_QPM_PERIOD_UM = 46.1  # 46.1-46.22 um across Evans/Graffitti/Laudenbach


def _group_delays_ps_per_mm():
    return phase_matched_group_delays_ps_per_mm(
        LAMBDA_PUMP, LAMBDA_DEGENERATE, LAMBDA_DEGENERATE, THETA, PLANE,
        "o", "o", "e", KTP["x"], KTP["y"], KTP["z"],
    )


def test_signal_and_idler_group_delay_match_literature():
    """GD_s (o, 1582nm) and GD_i (e, 1582nm) depend only on n_y and n_z
    (x-axis propagation never involves n_x). Both match Laudenbach et al.
    2017 to well under 1%."""
    _, GD_s, GD_i = _group_delays_ps_per_mm()
    assert abs(GD_s - LIT_GD_S) / LIT_GD_S < 0.01
    assert abs(GD_i - LIT_GD_I) / LIT_GD_I < 0.01


def test_pump_group_delay_matches_literature():
    """GD_p (o, 791nm) now uses Konig & Wong's own n_y (Eq. 3), sourced
    directly from their paper (papers/09.) rather than a proxy fit.
    Matches Laudenbach et al. 2017 to well under 1% (previously a
    documented ~2.9% gap with the Kato-Takaoka n_y proxy)."""
    GD_p, _, _ = _group_delays_ps_per_mm()
    rel_err = abs(GD_p - LIT_GD_P) / LIT_GD_P
    assert rel_err < 0.01, f"GD_p={GD_p:.4f} ps/mm vs literature {LIT_GD_P}; relative error {rel_err:.2%}"


def test_dispersion_parameter_reaches_symmetric_regime():
    """With Konig & Wong's own n_y, the GVM regime at the literature's
    design point is now correctly classified as symmetric/EPM (D~1,
    theta_PM near +-45 deg mod 180), matching the literature -- this was
    the central open item from the previous iteration of this milestone,
    now closed."""
    GD_p, GD_s, GD_i = _group_delays_ps_per_mm()
    D = -(GD_p - GD_s) / (GD_p - GD_i)
    tpm = theta_pm_deg(GD_p, GD_s, GD_i)
    assert abs(D - LIT_D) / LIT_D < 0.01, f"D={D:.4f} vs literature {LIT_D}"
    # theta_PM is only defined mod 180 deg; the symmetric locus is at
    # +-45 deg or, equivalently, +-135 deg.
    dist_to_45 = min(abs(abs(tpm) - 45.0), abs(abs(tpm) - 135.0))
    assert dist_to_45 < 2.0, f"theta_PM={tpm:.2f} deg not near the 45 deg symmetric locus"


def test_qpm_period_matches_literature():
    """Grating period for QPM at the literature's design point now
    matches Laudenbach et al. 2017 / Evans et al. / Graffitti et al.
    (46.1-46.22 um) to well under 1%."""
    Lambda = find_qpm_period_um(THETA, LAMBDA_PUMP, LAMBDA_DEGENERATE, LAMBDA_DEGENERATE,
                                 PLANE, "o", "o", "e", KTP["x"], KTP["y"], KTP["z"])
    rel_err = abs(Lambda - LIT_QPM_PERIOD_UM) / LIT_QPM_PERIOD_UM
    assert rel_err < 0.01, f"QPM period={Lambda:.3f} um vs literature {LIT_QPM_PERIOD_UM}"
