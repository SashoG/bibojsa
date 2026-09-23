"""Joint spectral amplitude (JSA), Schmidt decomposition, heralded purity.

    f(omega_s, omega_i) = alpha(omega_s + omega_i) * phi(omega_s, omega_i)

alpha: pump envelope (Gaussian or sech^2 amplitude, parametrized by the
       pump's *intensity* FWHM bandwidth or, equivalently, pulse duration).
phi:   phase-matching function, phi = sinc(Delta k L / 2) * exp(i Delta k
       L / 2), with the exp(...) phase included by default (flag
       `include_phase` to drop it -- the "no phase" convention some JSA
       treatments use, which affects the SVD only through an overall
       phase per (omega_s, omega_i) pair. i.e. it does affect |f|
       nowhere, but does affect the *complex* Schmidt decomposition; the
       amplitude-squared marginal/purity depend on it because the SVD is
       taken of the complex f, not |f|).

Delta k(omega_s, omega_i) is evaluated *exactly* (not Taylor-expanded)
from the real Sellmeier dispersion via `phasematch.delta_k`, at the fixed
crystal-cut angle theta found at the design (typically degenerate)
frequencies -- i.e. the crystal orientation is fixed; only the two SPDC
frequencies vary across the grid, with omega_p = omega_s + omega_i.

Schmidt decomposition: f sampled on an (Ns x Ni) grid is treated as a
matrix; its SVD f = U diag(s) V^dagger gives Schmidt coefficients
lambda_n = s_n^2 / sum(s^2), heralded purity P = sum(lambda_n^2), Schmidt
number K = 1/P.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .phasematch import ThreeWaveConfig, delta_k


def gaussian_pump_sigma(fwhm_omega: float) -> float:
    """sigma for amplitude envelope alpha(Omega) = exp(-Omega^2/(2 sigma^2))
    given the *intensity* (|alpha|^2) FWHM in angular frequency.

    |alpha(Omega)|^2 = exp(-Omega^2/sigma^2); solving exp(-Omega_h^2/sigma^2)
    = 1/2 gives Omega_h = sigma*sqrt(ln 2), so FWHM = 2*sigma*sqrt(ln 2).
    """
    return fwhm_omega / (2.0 * np.sqrt(np.log(2.0)))


def fwhm_omega_from_fwhm_wavelength(fwhm_lambda_um: float, lambda0_um: float) -> float:
    """Convert a pump *wavelength* intensity FWHM (um) at central wavelength
    lambda0_um into an angular-frequency intensity FWHM (rad/s), via the
    local Jacobian |domega/dlambda| = 2*pi*c/lambda0^2."""
    from .sellmeier import C_LIGHT
    lam0_m = lambda0_um * 1e-6
    dlam_m = fwhm_lambda_um * 1e-6
    return 2.0 * np.pi * C_LIGHT / lam0_m**2 * dlam_m


def pump_envelope_gaussian(omega_s: np.ndarray, omega_i: np.ndarray,
                            omega_p0: float, sigma: float) -> np.ndarray:
    """alpha(Omega) = exp(-(omega_s+omega_i-omega_p0)^2 / (2 sigma^2))."""
    Omega = omega_s + omega_i - omega_p0
    return np.exp(-Omega**2 / (2.0 * sigma**2))


def pump_envelope_sech2(omega_s: np.ndarray, omega_i: np.ndarray,
                         omega_p0: float, fwhm_omega: float) -> np.ndarray:
    """alpha(Omega) = sech(Omega * tau_sp), tau_sp chosen so that
    |alpha|^2 = sech^2(Omega*tau_sp) has the specified intensity FWHM.
    sech^2(x)=1/2 at x = arccosh(sqrt(2)) ~= 0.881374; FWHM_Omega =
    2*0.881374/tau_sp => tau_sp = 1.762747/FWHM_Omega.
    """
    tau_sp = 1.762747 / fwhm_omega
    Omega = omega_s + omega_i - omega_p0
    return 1.0 / np.cosh(Omega * tau_sp)


def phase_matching_function(omega_s: np.ndarray, omega_i: np.ndarray, theta: float,
                             config: ThreeWaveConfig, coeffs_x, coeffs_y, coeffs_z,
                             length_m: float, include_phase: bool = True) -> np.ndarray:
    """phi(omega_s, omega_i) = sinc(Delta k L / 2) [* exp(i Delta k L / 2)].

    Evaluated fully vectorized: `phasematch.delta_k` (and the whole
    Sellmeier / principal-plane index chain under it) is numpy-array-safe,
    so the (Ns x Ni) grid is computed in one shot rather than per point.
    """
    dk = delta_k(theta, omega_s + omega_i, omega_s, omega_i, config,
                 coeffs_x, coeffs_y, coeffs_z)
    x = dk * length_m / 2.0
    sinc = np.sinc(x / np.pi)  # np.sinc(y) = sin(pi y)/(pi y); y=x/pi -> sin(x)/x
    if include_phase:
        return sinc * np.exp(1j * x)
    return sinc.astype(complex)


@dataclass
class JSAResult:
    omega_s: np.ndarray
    omega_i: np.ndarray
    f: np.ndarray
    schmidt_lambdas: np.ndarray
    purity: float
    schmidt_number: float


def build_jsa(omega_s_axis: np.ndarray, omega_i_axis: np.ndarray, omega_p0: float,
              theta: float, config: ThreeWaveConfig, coeffs_x, coeffs_y, coeffs_z,
              length_m: float, pump_fwhm_omega: float, pump_shape: str = "gaussian",
              include_phase: bool = True) -> JSAResult:
    """Grid the JSA f(omega_s, omega_i) = alpha * phi over the given axes,
    and compute its Schmidt decomposition. Grids must be wide/fine enough
    for a converged purity -- see `purity_convergence_scan`."""
    Ws, Wi = np.meshgrid(omega_s_axis, omega_i_axis, indexing="ij")
    if pump_shape == "gaussian":
        sigma = gaussian_pump_sigma(pump_fwhm_omega)
        alpha = pump_envelope_gaussian(Ws, Wi, omega_p0, sigma)
    elif pump_shape == "sech2":
        alpha = pump_envelope_sech2(Ws, Wi, omega_p0, pump_fwhm_omega)
    else:
        raise ValueError(f"pump_shape must be 'gaussian' or 'sech2', got {pump_shape!r}")

    phi = phase_matching_function(Ws, Wi, theta, config, coeffs_x, coeffs_y, coeffs_z,
                                   length_m, include_phase=include_phase)
    f = alpha.astype(complex) * phi

    lambdas, purity, K = schmidt_decompose(f)
    return JSAResult(omega_s_axis, omega_i_axis, f, lambdas, purity, K)


def schmidt_decompose(f: np.ndarray) -> tuple[np.ndarray, float, float]:
    """SVD-based Schmidt decomposition of a gridded JSA matrix f.

    Returns (lambda_n, purity, schmidt_number). lambda_n are normalized
    (sum to 1) Schmidt eigenvalues from the singular values of f;
    purity = sum(lambda_n^2), schmidt_number K = 1/purity.
    Note: this treats the grid as a plain matrix (unweighted SVD), which
    is the standard approximation for a uniform grid; non-uniform grids
    would need singular-value decomposition of a suitably weighted
    kernel.
    """
    s = np.linalg.svd(f, compute_uv=False)
    weights = s**2
    total = np.sum(weights)
    lambdas = weights / total
    purity = float(np.sum(lambdas**2))
    K = 1.0 / purity
    return lambdas, purity, K


def apply_heralding_filter_gaussian(f: np.ndarray, omega_axis: np.ndarray,
                                     omega0: float, fwhm_omega: float,
                                     herald_axis: int) -> np.ndarray:
    """Multiply the JSA by a Gaussian *amplitude* filter (intensity FWHM
    `fwhm_omega`) on the heralding arm (axis 0 = signal/omega_s, axis 1 =
    idler/omega_i), modeling a spectral filter placed before the herald
    detector. Returns the filtered JSA (not renormalized)."""
    sigma = gaussian_pump_sigma(fwhm_omega)
    filt = np.exp(-(omega_axis - omega0) ** 2 / (2.0 * sigma**2))
    if herald_axis == 0:
        return f * filt[:, None]
    elif herald_axis == 1:
        return f * filt[None, :]
    else:
        raise ValueError("herald_axis must be 0 or 1")


def purity_convergence_scan(build_fn, window_scales: list[float], n_points: list[int]) -> list[dict]:
    """Convergence check: recompute purity over a set of (window, n_points)
    choices and report how it changes, per the brief's requirement to add
    a convergence check on purity vs grid size and window (higher-order
    curvature of the phase-matching band can matter near degeneracy).

    `build_fn(window_scale, n)` must return a JSAResult-like object with
    a `.purity` attribute for given window scale and grid size.
    """
    records = []
    for scale in window_scales:
        for n in n_points:
            result = build_fn(scale, n)
            records.append({"window_scale": scale, "n_points": n, "purity": result.purity})
    return records
