"""Collinear SPDC phase matching, principal-plane propagation.

omega_p = omega_s + omega_i (energy conservation)
Delta k = k_p - k_s - k_i = 0 (collinear momentum conservation)

Every validation milestone and BIBO survey target in this project uses
propagation confined to a principal dielectric plane, where each wave's
polarization is either the fixed ("o-like") or angle-dependent ("e-like")
eigenmode of `bibojsa.biaxial.principal_plane_index`. A `ThreeWaveConfig`
fixes which of the three principal planes and which polarization each of
pump/signal/idler takes; `find_phase_matching_angle` solves for the
crystal-cut angle theta that satisfies Delta k = 0 at specified
wavelengths (typically frequency degeneracy for the daughter photons).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .biaxial import k_of_omega_principal, PRINCIPAL_PLANES
from .sellmeier import lambda_um_to_omega


@dataclass(frozen=True)
class ThreeWaveConfig:
    """Collinear three-wave-mixing configuration for SPDC phase matching.

    pol_pump/signal/idler: 'o' (fixed/angle-independent index branch) or
    'e' (angle-dependent branch) in the chosen principal `plane`.

    qpm_period_um: optional quasi-phase-matching grating period (um) for
    periodically-poled crystals (e.g. PPKTP). When set, `delta_k` returns
    the *total* mismatch k_p - k_s - k_i - 2*pi/Lambda instead of the pure
    birefringent mismatch. None (default) = bulk birefringent phase
    matching only, no grating term.
    """
    plane: str
    pol_pump: str
    pol_signal: str
    pol_idler: str
    qpm_period_um: float | None = None

    def __post_init__(self):
        if self.plane not in PRINCIPAL_PLANES:
            raise ValueError(f"plane must be one of {PRINCIPAL_PLANES}, got {self.plane!r}")
        for name, pol in [("pol_pump", self.pol_pump), ("pol_signal", self.pol_signal),
                           ("pol_idler", self.pol_idler)]:
            if pol not in ("o", "e"):
                raise ValueError(f"{name} must be 'o' or 'e', got {pol!r}")
        if self.qpm_period_um is not None and self.qpm_period_um <= 0:
            raise ValueError(f"qpm_period_um must be positive, got {self.qpm_period_um!r}")


def delta_k(theta: float, omega_p: float, omega_s: float, omega_i: float,
            config: ThreeWaveConfig, coeffs_x, coeffs_y, coeffs_z) -> float:
    """Delta k(theta) = k_p - k_s - k_i [- 2*pi/Lambda if config.qpm_period_um
    is set] [rad/m] at fixed crystal-cut angle theta, for given
    pump/signal/idler angular frequencies (which need not satisfy
    omega_p = omega_s + omega_i if the caller wants an off-shell
    evaluation, e.g. for JSA phase-matching-function grids)."""
    kp = k_of_omega_principal(omega_p, theta, config.plane, config.pol_pump,
                               coeffs_x, coeffs_y, coeffs_z)
    ks = k_of_omega_principal(omega_s, theta, config.plane, config.pol_signal,
                               coeffs_x, coeffs_y, coeffs_z)
    ki = k_of_omega_principal(omega_i, theta, config.plane, config.pol_idler,
                               coeffs_x, coeffs_y, coeffs_z)
    dk = kp - ks - ki
    if config.qpm_period_um is not None:
        dk = dk - 2.0 * np.pi / (config.qpm_period_um * 1e-6)
    return dk


def find_qpm_period_um(theta: float, lambda_pump_um: float, lambda_signal_um: float,
                        lambda_idler_um: float, plane: str, pol_pump: str, pol_signal: str,
                        pol_idler: str, coeffs_x, coeffs_y, coeffs_z) -> float:
    """Grating period Lambda (um) that quasi-phase-matches the given
    (fixed) wavelengths and crystal-cut angle `theta`, i.e. solves
    Delta k_birefringent = 2*pi/Lambda for Lambda:

        Lambda = 2*pi / Delta k_birefringent

    Delta k_birefringent is k_p - k_s - k_i computed *without* any grating
    term (a plain `ThreeWaveConfig` with qpm_period_um=None). Returns a
    positive period regardless of the sign of Delta k_birefringent (only
    |Delta k| sets the required period; the sign is a grating-vector
    orientation convention, not a physical constraint on Lambda itself).
    """
    bare_config = ThreeWaveConfig(plane=plane, pol_pump=pol_pump, pol_signal=pol_signal,
                                   pol_idler=pol_idler, qpm_period_um=None)
    omega_p = lambda_um_to_omega(lambda_pump_um)
    omega_s = lambda_um_to_omega(lambda_signal_um)
    omega_i = lambda_um_to_omega(lambda_idler_um)
    dk_bir = delta_k(theta, omega_p, omega_s, omega_i, bare_config, coeffs_x, coeffs_y, coeffs_z)
    if dk_bir == 0:
        raise ValueError("Delta k_birefringent is exactly zero; no grating needed (Lambda -> infinity).")
    lambda_m = 2.0 * np.pi / abs(dk_bir)
    return lambda_m * 1e6


def find_phase_matching_angle(lambda_pump_um: float, lambda_signal_um: float,
                               lambda_idler_um: float, config: ThreeWaveConfig,
                               coeffs_x, coeffs_y, coeffs_z,
                               theta_bounds: tuple[float, float] = (1e-4, np.pi / 2 - 1e-4),
                               n_scan: int = 200) -> float:
    """Root-find the crystal-cut angle theta in `theta_bounds` (radians,
    measured per `principal_plane_index`'s convention for `config.plane`)
    at which Delta k = 0 for the given (fixed) wavelengths, by scanning
    for a sign change in Delta k(theta) and bisecting with `brentq`.

    Raises ValueError if no sign change is found (no phase-matching
    solution in the given angle range for these wavelengths/polarizations).
    """
    omega_p = lambda_um_to_omega(lambda_pump_um)
    omega_s = lambda_um_to_omega(lambda_signal_um)
    omega_i = lambda_um_to_omega(lambda_idler_um)

    def f(th):
        return delta_k(th, omega_p, omega_s, omega_i, config, coeffs_x, coeffs_y, coeffs_z)

    thetas = np.linspace(theta_bounds[0], theta_bounds[1], n_scan)
    vals = np.array([f(th) for th in thetas])
    sign_changes = np.where(np.diff(np.sign(vals)) != 0)[0]
    if len(sign_changes) == 0:
        raise ValueError(
            f"No phase-matching angle found in theta in "
            f"[{theta_bounds[0]:.4f}, {theta_bounds[1]:.4f}] rad for "
            f"lambda_pump={lambda_pump_um}, lambda_signal={lambda_signal_um}, "
            f"lambda_idler={lambda_idler_um} um, plane={config.plane}, "
            f"pols=({config.pol_pump},{config.pol_signal},{config.pol_idler})."
        )
    idx = sign_changes[0]
    return brentq(f, thetas[idx], thetas[idx + 1])
