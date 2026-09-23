"""Group-velocity-matching (GVM) classification.

The phase-matching function phi(omega_s, omega_i) = sinc(Delta k L/2) has,
near a phase-matched point, a ridge in the (omega_s, omega_i) plane whose
orientation is set by the group velocities of the three interacting waves.
Expanding Delta k to first order in (omega_s - omega_s0), (omega_i -
omega_i0) about the design point (using omega_p = omega_s + omega_i):

    Delta k ~= (k_p' - k_s')(omega_s - omega_s0) + (k_p' - k_i')(omega_i - omega_i0)

so the ridge (Delta k = 0 locus) has slope -(k_p' - k_s')/(k_p' - k_i') in
the (omega_i, omega_s) plane; this project (per the brief) reports it as
an orientation angle::

    tan(theta_PM) = -(k_p' - k_s') / (k_p' - k_i')

Special loci (asymmetric and symmetric GVM / extended phase matching):

    theta_PM = 0 deg or 90 deg  <=>  k_p' = k_s'  or  k_p' = k_i'
        (asymmetric GVM: the pump's group velocity matches one daughter
        exactly; the ridge is parallel to one frequency axis, so that
        daughter's marginal spectrum is set purely by the phase-matching
        bandwidth and is uncorrelated with the other -- e.g. Mosley et
        al., PRL 100, 133601 (2008), KDP type-II near 415/830 nm.)

    theta_PM = 45 deg  <=>  k_p' = (k_s' + k_i')/2
        (symmetric GVM / extended phase matching (EPM); e.g. Konig &
        Wong, APL 84, 1644 (2004), KTP type-II near 1582 nm.)

Group velocities (k') at each wave's own frequency, holding the crystal
cut angle theta fixed (theta is the *phase-matching angle*, a property of
the crystal cut -- not to be confused with theta_PM above, the JSA ridge
orientation), are obtained from `biaxial.k_of_omega_principal` via the
same converged finite-difference machinery used in `sellmeier.py`.
"""
from __future__ import annotations

import numpy as np

from .biaxial import k_of_omega_principal
from .sellmeier import finite_diff_converged


def group_index_at_angle(omega: float, theta: float, plane: str, pol: str,
                          coeffs_x, coeffs_y, coeffs_z,
                          h_rel: float = 1e-3, rtol: float = 1e-6,
                          return_diagnostics: bool = False):
    """k'(omega) [s/m] for a wave of polarization `pol` at fixed crystal
    cut angle `theta`, via 4th-order central finite difference with a
    step-size convergence check."""
    def f(w):
        return k_of_omega_principal(w, theta, plane, pol, coeffs_x, coeffs_y, coeffs_z)

    h0 = h_rel * omega
    val, h_used, converged, rel_change = finite_diff_converged(f, omega, order=1, h0=h0, rtol=rtol)
    if return_diagnostics:
        return val, h_used, converged, rel_change
    return val


def theta_pm_deg(kp_prime: float, ks_prime: float, ki_prime: float) -> float:
    """theta_PM (degrees), the phase-matching-function ridge orientation:

        tan(theta_PM) = -(k_p' - k_s') / (k_p' - k_i')

    Returned in (-90, 90] via arctan (the ridge orientation is only
    defined modulo 180 deg; callers interested in the 0/45/90 deg loci
    should look at |theta_PM| and 180-|theta_PM| as needed).
    """
    num = -(kp_prime - ks_prime)
    den = (kp_prime - ki_prime)
    return np.degrees(np.arctan2(num, den))


def gvm_regime(kp_prime: float, ks_prime: float, ki_prime: float,
               rel_tol: float = 0.02) -> str:
    """Classify which GVM locus (if any) the given group velocities are
    close to, comparing relative differences of k' against `rel_tol`:

    - 'asymmetric_signal': k_p' ~= k_s'  (theta_PM ~= 0 deg)
    - 'asymmetric_idler':  k_p' ~= k_i'  (theta_PM ~= 90 deg)
    - 'symmetric':         k_p' ~= (k_s'+k_i')/2  (theta_PM ~= 45 deg, EPM)
    - 'none': none of the above within `rel_tol`
    """
    scale = abs(kp_prime) if kp_prime != 0 else 1.0

    def close(a, b):
        return abs(a - b) / scale < rel_tol

    if close(kp_prime, ks_prime):
        return "asymmetric_signal"
    if close(kp_prime, ki_prime):
        return "asymmetric_idler"
    if close(kp_prime, 0.5 * (ks_prime + ki_prime)):
        return "symmetric"
    return "none"


def phase_matched_group_delays_ps_per_mm(lambda_pump_um: float, lambda_signal_um: float,
                                          lambda_idler_um: float, theta: float, plane: str,
                                          pol_pump: str, pol_signal: str, pol_idler: str,
                                          coeffs_x, coeffs_y, coeffs_z) -> tuple[float, float, float]:
    """Group delays GD_p, GD_s, GD_i (ps/mm) of the pump/signal/idler waves
    at a given (already phase-matched) crystal-cut angle `theta` and
    wavelength triplet, via n_g = c*k'(omega) and GD [ps/mm] = n_g /
    c[mm/ps] = n_g / 0.299792458. Packages the group-delay computation
    used by the KTP (milestone 2) and BIBO (milestone 3) validation
    tests into one reusable helper."""
    from .sellmeier import lambda_um_to_omega, C_LIGHT
    wp = lambda_um_to_omega(lambda_pump_um)
    ws = lambda_um_to_omega(lambda_signal_um)
    wi = lambda_um_to_omega(lambda_idler_um)
    ng_p = group_index_at_angle(wp, theta, plane, pol_pump, coeffs_x, coeffs_y, coeffs_z) * C_LIGHT
    ng_s = group_index_at_angle(ws, theta, plane, pol_signal, coeffs_x, coeffs_y, coeffs_z) * C_LIGHT
    ng_i = group_index_at_angle(wi, theta, plane, pol_idler, coeffs_x, coeffs_y, coeffs_z) * C_LIGHT
    return ng_p / 0.299792458, ng_s / 0.299792458, ng_i / 0.299792458


def find_symmetric_gvm_pump_wavelength(theta: float, plane: str, pol_pump: str,
                                        pol_signal: str, pol_idler: str,
                                        coeffs_x, coeffs_y, coeffs_z,
                                        lambda_pump_bounds: tuple[float, float],
                                        n_scan: int = 100) -> float:
    """Root-find the pump wavelength (um) for which the *degenerate*
    (lambda_s = lambda_i = 2*lambda_pump) symmetric-GVM condition holds
    at fixed crystal-cut angle `theta`:

        k_p' = (k_s' + k_i')/2   (theta_PM = 45 deg, extended phase matching)

    Scans `lambda_pump_bounds` (um) for a sign change in
    k_p' - (k_s'+k_i')/2 and bisects with `brentq`.
    """
    from scipy.optimize import brentq
    from .sellmeier import lambda_um_to_omega

    def mismatch(lambda_pump_um):
        lambda_deg_um = 2.0 * lambda_pump_um
        wp = lambda_um_to_omega(lambda_pump_um)
        ws = lambda_um_to_omega(lambda_deg_um)
        wi = lambda_um_to_omega(lambda_deg_um)
        kp_prime = group_index_at_angle(wp, theta, plane, pol_pump, coeffs_x, coeffs_y, coeffs_z)
        ks_prime = group_index_at_angle(ws, theta, plane, pol_signal, coeffs_x, coeffs_y, coeffs_z)
        ki_prime = group_index_at_angle(wi, theta, plane, pol_idler, coeffs_x, coeffs_y, coeffs_z)
        return kp_prime - 0.5 * (ks_prime + ki_prime)

    lo, hi = lambda_pump_bounds
    lams = np.linspace(lo, hi, n_scan)
    vals = np.array([mismatch(lam) for lam in lams])
    sign_changes = np.where(np.diff(np.sign(vals)) != 0)[0]
    if len(sign_changes) == 0:
        raise ValueError(
            f"No symmetric-GVM pump wavelength found in [{lo}, {hi}] um "
            f"for plane={plane}, pols=({pol_pump},{pol_signal},{pol_idler})."
        )
    idx = sign_changes[0]
    return brentq(mismatch, lams[idx], lams[idx + 1])
