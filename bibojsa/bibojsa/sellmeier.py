"""Sellmeier dispersion and dispersion derivatives.

Governing equations
--------------------
Two-pole Sellmeier form used throughout (Miyata-Umemura-Kato / LPR2010
convention, LPR2010 eq. in Sec. 2 and Table 2)::

    n^2(lambda) = A + B / (lambda^2 - C) - D * lambda^2

with lambda in vacuum micrometres.

Wavevector and its derivatives (user-facing API is vacuum wavelength in
micrometres; all internal frequency bookkeeping is in rad/s)::

    k(omega)  = n(omega) * omega / c                  [rad/m]
    k'(omega) = d k/d omega = n_g(omega) / c           [s/m], group index n_g
    k''(omega) = d^2 k/d omega^2                       [s^2/m], GVD

Derivatives are computed by 4th-order central finite differences directly
on k(omega) (not on n(lambda), to avoid an extra chain-rule step and its
associated round-off), with a step-size convergence check comparing the
result at h and h/2 (Richardson-style consistency check) -- see
`_finite_diff_converged`.

References
----------
[LPR2010] Petrov et al., Laser & Photonics Reviews 4, 53 (2010), Table 2
    (Sellmeier fits) and Table 3 (validation index values).
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

C_LIGHT = 2.99792458e8  # m/s, exact (SI definition)


@dataclass(frozen=True)
class SellmeierCoeffs:
    """Coefficients for one principal-axis Sellmeier fit.

    n^2(lambda_um) = A + B/(lambda_um^2 - C) - D*lambda_um^2

    Parameters
    ----------
    A : dimensionless
    B : micrometres^2
    C : micrometres^2
    D : micrometres^-2
    lambda_min, lambda_max : validity range in vacuum micrometres
    axis : which principal dielectric axis ('x', 'y', or 'z')
    source : literature citation string
    """

    A: float
    B: float
    C: float
    D: float
    lambda_min: float
    lambda_max: float
    axis: str = ""
    source: str = ""

    def n2(self, lambda_um: np.ndarray | float) -> np.ndarray | float:
        lam = np.asarray(lambda_um, dtype=float)
        return self.A + self.B / (lam**2 - self.C) - self.D * lam**2

    def n(self, lambda_um: np.ndarray | float) -> np.ndarray | float:
        n2 = self.n2(lambda_um)
        if np.any(np.asarray(n2) < 0):
            raise ValueError(
                f"Sellmeier n^2 < 0 for axis {self.axis!r} at some requested "
                f"wavelength(s); outside physical validity range "
                f"[{self.lambda_min}, {self.lambda_max}] um."
            )
        return np.sqrt(n2)

    def check_range(self, lambda_um: np.ndarray | float) -> None:
        lam = np.asarray(lambda_um, dtype=float)
        if np.any(lam < self.lambda_min) or np.any(lam > self.lambda_max):
            raise ValueError(
                f"Wavelength outside Sellmeier validity range "
                f"[{self.lambda_min}, {self.lambda_max}] um for axis "
                f"{self.axis!r} ({self.source})."
            )


def omega_to_lambda_um(omega: np.ndarray | float) -> np.ndarray | float:
    """Vacuum angular frequency (rad/s) -> vacuum wavelength (um)."""
    return 2.0 * np.pi * C_LIGHT / np.asarray(omega, dtype=float) * 1e6


def lambda_um_to_omega(lambda_um: np.ndarray | float) -> np.ndarray | float:
    """Vacuum wavelength (um) -> vacuum angular frequency (rad/s)."""
    lam_m = np.asarray(lambda_um, dtype=float) * 1e-6
    return 2.0 * np.pi * C_LIGHT / lam_m


def refractive_index(lambda_um: np.ndarray | float, coeffs: SellmeierCoeffs,
                      enforce_range: bool = True) -> np.ndarray | float:
    """Refractive index n(lambda) for one principal axis."""
    if enforce_range:
        coeffs.check_range(lambda_um)
    return coeffs.n(lambda_um)


def k_of_omega(omega: np.ndarray | float, coeffs: SellmeierCoeffs,
               enforce_range: bool = True) -> np.ndarray | float:
    """Wavevector magnitude k(omega) = n(omega) * omega / c, in rad/m."""
    lam_um = omega_to_lambda_um(omega)
    n = refractive_index(lam_um, coeffs, enforce_range=enforce_range)
    return n * np.asarray(omega, dtype=float) / C_LIGHT


def _central_diff_1st(f, x: float, h: float) -> float:
    """4th-order central finite difference for f'(x)."""
    return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h)


def _central_diff_2nd(f, x: float, h: float) -> float:
    """4th-order central finite difference for f''(x)."""
    return (
        -f(x + 2 * h) + 16 * f(x + h) - 30 * f(x) + 16 * f(x - h) - f(x - 2 * h)
    ) / (12 * h**2)


def _finite_diff_converged(f, x: float, order: int, h0: float,
                            rtol: float = 1e-6, max_halvings: int = 12):
    """Halve the step size until the derivative estimate stabilizes.

    Returns (value, h_used, converged: bool, achieved_rel_change: float).
    Comparing successive halvings guards against both truncation error
    (h too large) and cancellation/round-off error (h too small); we stop
    at the h that gives the best agreement between successive estimates.
    """
    diff_fn = _central_diff_1st if order == 1 else _central_diff_2nd
    h = h0
    prev = diff_fn(f, x, h)
    best = (prev, h, False, np.inf)
    for _ in range(max_halvings):
        h_new = h / 2
        cur = diff_fn(f, x, h_new)
        denom = abs(cur) if cur != 0 else 1.0
        rel_change = abs(cur - prev) / denom
        if rel_change < best[3]:
            best = (cur, h_new, rel_change < rtol, rel_change)
        if rel_change < rtol:
            return cur, h_new, True, rel_change
        prev = cur
        h = h_new
    return best


def group_index(omega: float, coeffs: SellmeierCoeffs, h_rel: float = 1e-3,
                 rtol: float = 1e-6) -> float:
    """n_g(omega) = c * k'(omega), group index, via converged finite diff."""
    kp, _, converged, _ = kprime(omega, coeffs, h_rel=h_rel, rtol=rtol,
                                  return_diagnostics=True)
    return kp * C_LIGHT


def kprime(omega: float, coeffs: SellmeierCoeffs, h_rel: float = 1e-3,
           rtol: float = 1e-6, return_diagnostics: bool = False):
    """k'(omega) = dk/domega [s/m], 4th-order central FD with convergence check."""
    def f(w):
        return k_of_omega(w, coeffs, enforce_range=False)

    h0 = h_rel * omega
    val, h_used, converged, rel_change = _finite_diff_converged(
        f, omega, order=1, h0=h0, rtol=rtol
    )
    if return_diagnostics:
        return val, h_used, converged, rel_change
    return val


def kdoubleprime(omega: float, coeffs: SellmeierCoeffs, h_rel: float = 1e-3,
                  rtol: float = 1e-6, return_diagnostics: bool = False):
    """k''(omega) = d^2k/domega^2 [s^2/m], GVD, 4th-order central FD."""
    def f(w):
        return k_of_omega(w, coeffs, enforce_range=False)

    h0 = h_rel * omega
    val, h_used, converged, rel_change = _finite_diff_converged(
        f, omega, order=2, h0=h0, rtol=rtol
    )
    if return_diagnostics:
        return val, h_used, converged, rel_change
    return val


def optic_axis_angle_deg(lambda_um: float, coeffs_x: SellmeierCoeffs,
                          coeffs_y: SellmeierCoeffs,
                          coeffs_z: SellmeierCoeffs) -> float:
    """Optic axis half-angle Omega (degrees) for a biaxial crystal.

    Omega = arcsin[ (nz/ny) * sqrt((ny^2 - nx^2) / (nz^2 - nx^2)) ]

    Convention n_x < n_y < n_z. Returns Omega in degrees; report 2*Omega
    (the full optic-axis angle, as tabulated in LPR2010) at the call site.
    """
    nx = refractive_index(lambda_um, coeffs_x)
    ny = refractive_index(lambda_um, coeffs_y)
    nz = refractive_index(lambda_um, coeffs_z)
    arg = (nz / ny) * np.sqrt((ny**2 - nx**2) / (nz**2 - nx**2))
    return np.degrees(np.arcsin(arg))


# Public alias: the step-size-convergence-checked finite-difference helper
# is reused by biaxial.py (walk-off, angle derivatives) and gvm.py (group
# velocities at fixed phase-matching angle).
finite_diff_converged = _finite_diff_converged


@dataclass(frozen=True)
class ZernikeSellmeierCoeffs:
    """Coefficients for the 'double-resonance' Sellmeier form used by
    Zernike (1964) for KDP/ADP and commonly for KTP-family crystals:

        n^2(lambda_um) = A + B*lambda_um^2/(lambda_um^2 - C) + D/(lambda_um^2 - E)

    This is a different functional form from `SellmeierCoeffs` (the
    LPR2010/Miyata-Umemura-Kato two-pole form used for BIBO/BBO); kept as
    a separate class rather than unified because the two forms are not
    algebraically equivalent term-by-term, and conflating them invites
    transcription errors. Duck-types the same `.n()` / `.check_range()`
    interface as `SellmeierCoeffs` so both can be used interchangeably by
    `refractive_index`, `k_of_omega`, `kprime`, `kdoubleprime`, etc.
    """

    A: float
    B: float
    C: float
    D: float
    E: float
    lambda_min: float
    lambda_max: float
    axis: str = ""
    source: str = ""

    def n2(self, lambda_um: np.ndarray | float) -> np.ndarray | float:
        lam = np.asarray(lambda_um, dtype=float)
        return self.A + self.B * lam**2 / (lam**2 - self.C) + self.D / (lam**2 - self.E)

    def n(self, lambda_um: np.ndarray | float) -> np.ndarray | float:
        n2 = self.n2(lambda_um)
        if np.any(np.asarray(n2) < 0):
            raise ValueError(
                f"Sellmeier n^2 < 0 for axis {self.axis!r} at some requested "
                f"wavelength(s); outside physical validity range "
                f"[{self.lambda_min}, {self.lambda_max}] um."
            )
        return np.sqrt(n2)

    def check_range(self, lambda_um: np.ndarray | float) -> None:
        lam = np.asarray(lambda_um, dtype=float)
        if np.any(lam < self.lambda_min) or np.any(lam > self.lambda_max):
            raise ValueError(
                f"Wavelength outside Sellmeier validity range "
                f"[{self.lambda_min}, {self.lambda_max}] um for axis "
                f"{self.axis!r} ({self.source})."
            )


@dataclass(frozen=True)
class TwoPoleLinearSellmeierCoeffs:
    """Coefficients for the 'two-pole + linear' Sellmeier form used by
    Fradkin et al., Appl. Phys. Lett. 74, 914 (1999) for KTP n_z (as
    reproduced in Heberle & Levine, J. Res. NIST 118, 375 (2013), who cite
    Konig & Wong (2004) n_y and Fradkin et al. (1999) n_z as the standard
    combination used in the GVM/EPM KTP literature; here paired with
    Kato & Takaoka (2002) n_x, n_y per the same recipe):

        n^2(lambda_um) = A + B1*lambda_um^2/(lambda_um^2 - C1)
                           + B2*lambda_um^2/(lambda_um^2 - C2)
                           - D*lambda_um^2

    Duck-types the same `.n()` / `.n2()` / `.check_range()` interface as
    `SellmeierCoeffs` / `ZernikeSellmeierCoeffs`.
    """

    A: float
    B1: float
    C1: float
    B2: float
    C2: float
    D: float
    lambda_min: float
    lambda_max: float
    axis: str = ""
    source: str = ""

    def n2(self, lambda_um: np.ndarray | float) -> np.ndarray | float:
        lam = np.asarray(lambda_um, dtype=float)
        lam2 = lam**2
        return (
            self.A
            + self.B1 * lam2 / (lam2 - self.C1)
            + self.B2 * lam2 / (lam2 - self.C2)
            - self.D * lam2
        )

    def n(self, lambda_um: np.ndarray | float) -> np.ndarray | float:
        n2 = self.n2(lambda_um)
        if np.any(np.asarray(n2) < 0):
            raise ValueError(
                f"Sellmeier n^2 < 0 for axis {self.axis!r} at some requested "
                f"wavelength(s); outside physical validity range "
                f"[{self.lambda_min}, {self.lambda_max}] um."
            )
        return np.sqrt(n2)

    def check_range(self, lambda_um: np.ndarray | float) -> None:
        lam = np.asarray(lambda_um, dtype=float)
        if np.any(lam < self.lambda_min) or np.any(lam > self.lambda_max):
            raise ValueError(
                f"Wavelength outside Sellmeier validity range "
                f"[{self.lambda_min}, {self.lambda_max}] um for axis "
                f"{self.axis!r} ({self.source})."
            )
