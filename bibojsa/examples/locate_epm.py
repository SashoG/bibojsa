"""Locate the exact symmetric-GVM (theta_PM = 45 deg) operating point for
BIBO x-z ooe SPDC pumped at 800 nm, near the milestone-3 candidate
lambda_s ~ 1.34 um, and print the design-point summary.

Result (recorded 2026-07-12, this Sellmeier fit):
    lambda_s = 1.3409 um, lambda_i = 1.9831 um, theta_cut = 10.9577 deg
    GD (ps/mm): p = 6.0312, s = 6.0288, i = 6.0336  (triple GVM, ~2.4 fs/mm)
    |deff| = 3.1417 pm/V; pump walk-off rho ~ 30.5 mrad (1.75 deg)
"""
import numpy as np
from scipy.optimize import brentq

from bibojsa.materials.bibo import BIBO
from bibojsa.phasematch import ThreeWaveConfig, find_phase_matching_angle
from bibojsa.gvm import group_index_at_angle, theta_pm_deg
from bibojsa.deff import deff_bibo_from_config
from bibojsa.biaxial import walkoff_angle
from bibojsa.sellmeier import lambda_um_to_omega

PLANE = "xz"
LAMBDA_PUMP = 0.800
CONFIG = ThreeWaveConfig(plane=PLANE, pol_pump="e", pol_signal="o", pol_idler="o")
THETA_BOUNDS = (np.radians(0.05), np.radians(24.9))
X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]


def idler_um(lam_s):
    return 1.0 / (1.0 / LAMBDA_PUMP - 1.0 / lam_s)


def point(lam_s):
    lam_i = idler_um(lam_s)
    theta = find_phase_matching_angle(LAMBDA_PUMP, lam_s, lam_i, CONFIG,
                                       X, Y, Z, theta_bounds=THETA_BOUNDS, n_scan=400)
    wp, ws, wi = (lambda_um_to_omega(l) for l in (LAMBDA_PUMP, lam_s, lam_i))
    kp = group_index_at_angle(wp, theta, PLANE, "e", X, Y, Z)
    ks = group_index_at_angle(ws, theta, PLANE, "o", X, Y, Z)
    ki = group_index_at_angle(wi, theta, PLANE, "o", X, Y, Z)
    return theta, lam_i, kp, ks, ki


def symmetric_mismatch(lam_s):
    """kp' - (ks'+ki')/2 = 0  <=>  theta_PM = 45 deg."""
    _, _, kp, ks, ki = point(lam_s)
    return kp - 0.5 * (ks + ki)


if __name__ == "__main__":
    for lam in (1.25, 1.30, 1.34, 1.38, 1.42):
        th, li, kp, ks, ki = point(lam)
        print(f"lam_s={lam:.3f} um  lam_i={li:.3f}  theta={np.degrees(th):7.3f} deg  "
              f"theta_PM={theta_pm_deg(kp, ks, ki):7.2f} deg  "
              f"mismatch={(kp-0.5*(ks+ki))*1e15:+.3f} fs/m")

    lam_star = brentq(symmetric_mismatch, 1.25, 1.42, xtol=1e-6)
    th, li, kp, ks, ki = point(lam_star)
    tpm = theta_pm_deg(kp, ks, ki)
    d = deff_bibo_from_config(CONFIG, th, lambda_pump_um=LAMBDA_PUMP)
    rho_p = walkoff_angle(LAMBDA_PUMP, th, PLANE, X, Y, Z)  # pump is the only e-wave
    print("\n=== symmetric-GVM (EPM) design point ===")
    print(f"lambda_s = {lam_star:.4f} um, lambda_i = {li:.4f} um, pump 800 nm")
    print(f"theta_cut = {np.degrees(th):.4f} deg")
    print(f"theta_PM = {tpm:.3f} deg (mod 180: {tpm % 180.0:.3f})")
    print(f"GD (ps/mm): p={kp*2.99792458e8/0.299792458:.4f} s={ks*2.99792458e8/0.299792458:.4f} i={ki*2.99792458e8/0.299792458:.4f}")
    print(f"deff = {d:+.4f} pm/V (|deff| = {abs(d):.4f})")
    print(f"pump walk-off rho = {rho_p*1e3:.1f} mrad ({np.degrees(rho_p):.3f} deg)")
