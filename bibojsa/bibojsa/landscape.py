"""Continuous-pump GVM landscape: the phase-matching surface and its
group-velocity-matching level curves.

Motivation
----------
Collinear SPDC in a fixed interaction (principal plane + o/e sequence) is a
*two*-parameter family: choose a pump wavelength and a signal wavelength and
the crystal-cut angle theta follows from Delta k = 0.  The joint-spectral
ridge orientation theta_PM is therefore a scalar field on the
(lambda_p, lambda_s) phase-matching surface, and the engineering targets

    theta_PM =  0 deg   (asymmetric GVM, k_p' = k_s')
    theta_PM = 45 deg   (symmetric GVM / extended phase matching)
    theta_PM = 90 deg   (asymmetric GVM, k_p' = k_i')

are *level curves* of that field, not isolated points.  Surveying at a fixed
list of pump wavelengths samples those curves at arbitrary places and can
miss the parts of them that matter -- in particular the frequency-degenerate
point, where each curve crosses lambda_s = 2*lambda_p.  This module traces
the curves themselves.

Everything here is a vectorised re-expression of the scalar primitives in
`bibojsa.phasematch`, `bibojsa.biaxial` and `bibojsa.gvm`; the agreement is
asserted in `tests/test_landscape.py` at the 1e-9 level, so the landscape
results and the point-by-point results in the rest of the package are the
same computation at different scales.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .biaxial import k_of_omega_principal
from .phasematch import ThreeWaveConfig
from .sellmeier import (C_LIGHT, lambda_um_to_omega, omega_to_lambda_um,
                        optic_axis_angle_deg)

# Sellmeier validity of the BIBO fit is 0.32-3.09 um; keep a margin so that
# the finite-difference stencils below never step outside it.
LAM_MIN, LAM_MAX = 0.33, 3.05

# Relative finite-difference step for the vectorised k'.  Chosen by the
# convergence study in tests/test_landscape.py: it reproduces the
# step-halving result of `gvm.group_index_at_angle` to better than 1e-10
# relative over the whole 0.33-3.05 um range.
KPRIME_H_REL = 5e-4


# ---------------------------------------------------------------- primitives

def kprime_vec(omega, theta, plane: str, pol: str, cx, cy, cz,
               h_rel: float = KPRIME_H_REL):
    """k'(omega) [s/m] at fixed crystal-cut angle, vectorised over arrays.

    Fourth-order central difference on k(omega) at a fixed relative step,
    the array-safe counterpart of `gvm.group_index_at_angle` (which halves
    the step until the estimate stabilises; that loop cannot vectorise
    because its stopping test is per-element).
    """
    omega = np.asarray(omega, dtype=float)
    h = h_rel * omega

    def k(w):
        return k_of_omega_principal(w, theta, plane, pol, cx, cy, cz)

    return (-k(omega + 2 * h) + 8 * k(omega + h)
            - 8 * k(omega - h) + k(omega - 2 * h)) / (12 * h)


def delta_k_vec(theta, omega_p, omega_s, omega_i, plane: str, seq: str,
                cx, cy, cz):
    """Delta k = k_p - k_s - k_i [rad/m], broadcasting over theta and omegas.

    `seq` is the LPR2010 o/e sequence idler-signal-pump.
    """
    kp = k_of_omega_principal(omega_p, theta, plane, seq[2], cx, cy, cz)
    ks = k_of_omega_principal(omega_s, theta, plane, seq[1], cx, cy, cz)
    ki = k_of_omega_principal(omega_i, theta, plane, seq[0], cx, cy, cz)
    return kp - ks - ki


def theta_pm_deg_vec(kp, ks, ki):
    """JSA ridge orientation in the (omega_s, omega_i) plane, degrees.

    tan(theta_PM) = d omega_i / d omega_s along the Delta k = 0 locus
                  = -(k_p' - k_s') / (k_p' - k_i')

    i.e. the angle is measured from the *signal* frequency axis.  Returned
    on (-90, 90]; 0 deg and 90 deg are the two asymmetric-GVM orientations
    and +45 deg is extended phase matching.
    """
    return np.degrees(np.arctan2(-(kp - ks), (kp - ki)))


def angle_bounds(plane: str):
    """Allowed crystal-cut angle range (radians) for a principal plane."""
    if plane == "xy":
        return np.radians(0.05), np.radians(89.95)
    return np.radians(0.05), np.radians(179.95)


def region_ok(plane: str, seq: str, theta_deg, lam_p_um, cx, cy, cz):
    """Whether a root lies in the angular region where this interaction is
    the physically realised one.

    In the x-z plane the eigenmode ordering swaps at the optic axes: the
    "+" interactions (eeo, oeo) exist between them, the "-" interactions
    (ooe, eoe) outside.  The y-z and x-y planes have no such split.
    """
    if plane != "xz":
        return np.ones_like(np.asarray(theta_deg, dtype=float), dtype=bool)
    om = optic_axis_angle_deg(lam_p_um, cx, cy, cz)
    inside = (np.asarray(theta_deg) > om) & (np.asarray(theta_deg) < 180.0 - om)
    return inside if seq in ("eeo", "oeo") else ~inside


# ------------------------------------------------------- phase-matching roots

def solve_theta_grid(lam_p, lam_s, plane: str, seq: str, cx, cy, cz,
                     n_theta: int = 1801, branch: str = "first",
                     chunk: int = 4000):
    """Crystal-cut angles solving Delta k = 0, over arrays of (lam_p, lam_s).

    Scans Delta k on a theta grid, brackets every sign change, and refines
    by 60 rounds of vectorised bisection (machine precision).  Returns an
    array of theta in radians, NaN where no root exists in the allowed
    range or where the wavelengths fall outside the Sellmeier window.

    `branch` selects which root to keep when the scan finds several:
    "first" (smallest theta) or "second".  In the x-z plane the mirror
    partner at 180-theta always exists and is physically equivalent, so
    only theta <= 90 deg is retained there.
    """
    scalar = np.isscalar(lam_p) and np.isscalar(lam_s)
    lam_p = np.atleast_1d(np.asarray(lam_p, dtype=float))
    lam_s = np.atleast_1d(np.asarray(lam_s, dtype=float))
    lam_p, lam_s = np.broadcast_arrays(lam_p, lam_s)
    shape = lam_p.shape
    lam_p, lam_s = lam_p.ravel(), lam_s.ravel()

    lam_i = np.full_like(lam_p, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = 1.0 / lam_p - 1.0 / lam_s
        lam_i = np.where(denom > 0, 1.0 / denom, np.nan)

    good = (np.isfinite(lam_i) & (lam_i >= LAM_MIN) & (lam_i <= LAM_MAX)
            & (lam_s >= LAM_MIN) & (lam_s <= LAM_MAX)
            & (lam_p >= LAM_MIN) & (lam_p <= LAM_MAX)
            & (lam_i >= lam_s))
    out = np.full(lam_p.shape, np.nan)
    idx_all = np.flatnonzero(good)
    if idx_all.size == 0:
        if scalar:
            return float(out[0]), float(lam_i[0])
        return out.reshape(shape), lam_i.reshape(shape)

    lo, hi = angle_bounds(plane)
    if plane == "xz":
        hi = np.radians(90.0)          # mirror octant is equivalent
    grid = np.linspace(lo, hi, n_theta)

    for start in range(0, idx_all.size, chunk):
        idx = idx_all[start:start + chunk]
        wp = lambda_um_to_omega(lam_p[idx])[:, None]
        ws = lambda_um_to_omega(lam_s[idx])[:, None]
        wi = lambda_um_to_omega(lam_i[idx])[:, None]
        vals = delta_k_vec(grid[None, :], wp, ws, wi, plane, seq, cx, cy, cz)

        sgn = np.signbit(vals)
        cross = sgn[:, :-1] != sgn[:, 1:]
        has = cross.any(axis=1)
        if not has.any():
            continue
        which = np.argmax(cross, axis=1)
        if branch == "second":
            second = np.where(cross.sum(axis=1) >= 2)
            for r in second[0]:
                cols = np.flatnonzero(cross[r])
                which[r] = cols[1]
            has = has & (cross.sum(axis=1) >= 2)

        a = grid[which].copy()
        b = grid[which + 1].copy()
        fa = np.take_along_axis(vals, which[:, None], 1).ravel()
        for _ in range(60):
            m = 0.5 * (a + b)
            fm = delta_k_vec(m, wp.ravel(), ws.ravel(), wi.ravel(),
                             plane, seq, cx, cy, cz)
            left = np.signbit(fm) == np.signbit(fa)
            a = np.where(left, m, a)
            fa = np.where(left, fm, fa)
            b = np.where(left, b, m)
        theta = 0.5 * (a + b)
        theta = np.where(has, theta, np.nan)

        ok = region_ok(plane, seq, np.degrees(theta), lam_p[idx], cx, cy, cz)
        out[idx] = np.where(ok, theta, np.nan)

    if scalar:
        return float(out[0]), float(lam_i[0])
    return out.reshape(shape), lam_i.reshape(shape)


def gvm_functions(lam_p, lam_s, lam_i, theta, plane: str, seq: str,
                  cx, cy, cz):
    """The three GVM mismatch functions and theta_PM at given roots.

    Returns a dict with entries 'g0', 'g45', 'g90' (units s/m; each is zero
    exactly on its GVM condition) plus 'theta_pm' in degrees and the three
    group slownesses.
    """
    wp = lambda_um_to_omega(np.asarray(lam_p, dtype=float))
    ws = lambda_um_to_omega(np.asarray(lam_s, dtype=float))
    wi = lambda_um_to_omega(np.asarray(lam_i, dtype=float))
    kp = kprime_vec(wp, theta, plane, seq[2], cx, cy, cz)
    ks = kprime_vec(ws, theta, plane, seq[1], cx, cy, cz)
    ki = kprime_vec(wi, theta, plane, seq[0], cx, cy, cz)
    return {
        "kp": kp, "ks": ks, "ki": ki,
        "g0": kp - ks,
        "g90": kp - ki,
        "g45": kp - 0.5 * (ks + ki),
        "theta_pm": theta_pm_deg_vec(kp, ks, ki),
    }


# --------------------------------------------------------- level-curve tracing

@dataclass
class GVMPoint:
    plane: str
    seq: str
    target: int            # 0, 45 or 90 degrees
    lam_p_um: float
    lam_s_um: float
    lam_i_um: float
    theta_deg: float
    theta_pm_deg: float
    branch: str
    kp: float
    ks: float
    ki: float

    @property
    def degeneracy(self) -> float:
        """|lambda_i - lambda_s| / lambda_s; 0 at frequency degeneracy."""
        return abs(self.lam_i_um - self.lam_s_um) / self.lam_s_um


def trace_gvm_curve(plane: str, seq: str, target: int, cx, cy, cz,
                    lam_p_grid, n_ls: int = 260, branch: str = "first",
                    n_theta: int = 1801):
    """Trace one theta_PM = `target` level curve across a pump-wavelength grid.

    For each pump wavelength the signal axis is scanned, the phase-matching
    angle solved at every node, the appropriate GVM function evaluated, and
    every sign change refined by `brentq` on the scalar composition
    (solve theta, then evaluate g).  Returns a list of `GVMPoint`.
    """
    key = {0: "g0", 45: "g45", 90: "g90"}[target]
    points: list[GVMPoint] = []

    for lam_p in np.atleast_1d(lam_p_grid):
        ls_lo = max(lam_p * 1.0005, LAM_MIN)
        ls_hi = min(2.0 * lam_p, LAM_MAX)
        if ls_hi <= ls_lo:
            continue
        ls = np.linspace(ls_lo, ls_hi, n_ls)
        th, li = solve_theta_grid(np.full_like(ls, lam_p), ls, plane, seq,
                                  cx, cy, cz, n_theta=n_theta, branch=branch)
        ok = np.isfinite(th)
        if ok.sum() < 2:
            continue
        g = np.full_like(ls, np.nan)
        g[ok] = gvm_functions(np.full(ok.sum(), lam_p), ls[ok], li[ok],
                              th[ok], plane, seq, cx, cy, cz)[key]

        def scalar_g(x):
            t, lii = solve_theta_grid(float(lam_p), float(x), plane, seq,
                                      cx, cy, cz, n_theta=n_theta,
                                      branch=branch)
            if not np.isfinite(t):
                return np.nan
            return float(gvm_functions(lam_p, x, lii, t, plane, seq,
                                       cx, cy, cz)[key])

        for j in range(len(ls) - 1):
            ga, gb = g[j], g[j + 1]
            if not (np.isfinite(ga) and np.isfinite(gb)) or ga * gb >= 0:
                continue
            try:
                ls_star = brentq(scalar_g, ls[j], ls[j + 1], xtol=1e-9)
            except (ValueError, RuntimeError):
                continue
            t, lii = solve_theta_grid(float(lam_p), float(ls_star), plane,
                                      seq, cx, cy, cz, n_theta=n_theta,
                                      branch=branch)
            if not np.isfinite(t):
                continue
            gf = gvm_functions(lam_p, ls_star, lii, t, plane, seq, cx, cy, cz)
            points.append(GVMPoint(
                plane=plane, seq=seq, target=target, lam_p_um=float(lam_p),
                lam_s_um=float(ls_star), lam_i_um=lii,
                theta_deg=float(np.degrees(t)),
                theta_pm_deg=float(gf["theta_pm"]), branch=branch,
                kp=float(gf["kp"]), ks=float(gf["ks"]), ki=float(gf["ki"])))
    return points


def degenerate_gvm_point(plane: str, seq: str, target: int, cx, cy, cz,
                         lam_p_bounds, n_scan: int = 90, n_theta: int = 1801,
                         branch: str = "first"):
    """Pump wavelength at which a GVM curve crosses frequency degeneracy.

    Solves g_target(lambda_p) = 0 along lambda_s = lambda_i = 2*lambda_p.
    Returns a `GVMPoint`, or None if no crossing lies in `lam_p_bounds`.
    """
    key = {0: "g0", 45: "g45", 90: "g90"}[target]

    def g_of_pump(lp):
        lam_d = 2.0 * lp
        t, _ = solve_theta_grid(float(lp), float(lam_d) * 0.9999999, plane,
                                seq, cx, cy, cz, n_theta=n_theta,
                                branch=branch)
        if not np.isfinite(t):
            return np.nan
        return float(gvm_functions(lp, lam_d, lam_d, t, plane, seq,
                                   cx, cy, cz)[key])

    lps = np.linspace(lam_p_bounds[0], lam_p_bounds[1], n_scan)
    vals = np.array([g_of_pump(lp) for lp in lps])
    fin = np.isfinite(vals)
    for j in range(len(lps) - 1):
        if not (fin[j] and fin[j + 1]) or vals[j] * vals[j + 1] >= 0:
            continue
        try:
            lp_star = brentq(g_of_pump, lps[j], lps[j + 1], xtol=1e-10)
        except (ValueError, RuntimeError):
            continue
        lam_d = 2.0 * lp_star
        t, _ = solve_theta_grid(float(lp_star), float(lam_d) * 0.9999999,
                                plane, seq, cx, cy, cz, n_theta=n_theta,
                                branch=branch)
        gf = gvm_functions(lp_star, lam_d, lam_d, t, plane, seq, cx, cy, cz)
        return GVMPoint(plane=plane, seq=seq, target=target,
                        lam_p_um=float(lp_star), lam_s_um=lam_d,
                        lam_i_um=lam_d, theta_deg=float(np.degrees(t)),
                        theta_pm_deg=float(gf["theta_pm"]), branch=branch,
                        kp=float(gf["kp"]), ks=float(gf["ks"]),
                        ki=float(gf["ki"]))
    return None


def gd_ps_per_mm(kprime):
    """k' [s/m] -> group delay [ps/mm].  GD = k' * 1e-3 m/mm * 1e12 ps/s."""
    return np.asarray(kprime) * 1e9
