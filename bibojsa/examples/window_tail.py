"""Collection-window dependence of the unfiltered heralded purity.

An unapodized sinc phase-matching function has 1/x amplitude tails.  Along
the pump stripe (delta_omega_s = -delta_omega_i) the pump envelope does not
suppress them, so the weight the tails carry inside a square collection
window of half-width W falls only as 1/W, whichever way the phase-matching
ridge points.  The Schmidt purity computed on that window therefore obeys

    P(W) = P_inf + C / W

to high accuracy, and *every* unfiltered sinc-limited purity is a property
of the collection band, not of the crystal alone.  This script measures
P(W) at each featured operating point on square windows of half-width
W = 6...24 pump sigma at a fixed pixel size (so grid resolution is not a
variable), fits P_inf and C, and records the paper's own adaptive window
in the same units.  It also reports the filtered purities at two windows to
show that they are window-independent.

Outputs
  survey/window_tail.csv          per point: P at each W, P_inf, C, fit residual,
                                  the paper's window in sigma and in nm, and the
                                  filtered purities at 6 and 16 sigma
  survey/length_extension.csv     P (paper window) and P_inf for A, B, C at
                                  L = 10...30 mm over pump bandwidths
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from landscape_headline import FEATURED
from landscape_points import W_MAX, W_MIN, adaptive_axes, cfg_for

from bibojsa.jsa import (apply_heralding_filter_gaussian, build_jsa,
                         heralding_axis,
                         fwhm_omega_from_fwhm_wavelength, gaussian_pump_sigma,
                         schmidt_decompose)
from bibojsa.materials.bibo import BIBO
from bibojsa.sellmeier import C_LIGHT, lambda_um_to_omega

X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "survey")

W_SIGMAS = (6.0, 8.0, 12.0, 16.0, 24.0)
PIXEL_RAD_PS = 1.0          # checked against 0.5 rad/ps: |dP| < 3e-5


def axes_for(ws0, wi0, half, pixel):
    n = int(2 * half * 1e-12 / pixel) + 1
    ws = np.linspace(max(ws0 - half, W_MIN), min(ws0 + half, W_MAX), n)
    wi = np.linspace(max(wi0 - half, W_MIN), min(wi0 + half, W_MAX), n)
    clipped = (ws0 - half < W_MIN) or (ws0 + half > W_MAX) \
        or (wi0 - half < W_MIN) or (wi0 + half > W_MAX)
    return ws, wi, clipped


def purity_on(pt, L_mm, bw_nm, ws, wi):
    tag, plane, seq, _, pump_nm, lam_s, lam_i, theta_deg, *_ = pt
    cfg = cfg_for(plane, seq)
    wp0 = lambda_um_to_omega(pump_nm / 1000.0)
    fw = fwhm_omega_from_fwhm_wavelength(bw_nm * 1e-3, pump_nm / 1000.0)
    res = build_jsa(ws, wi, wp0, np.radians(theta_deg), cfg, X, Y, Z,
                    L_mm * 1e-3, fw, pump_shape="gaussian", include_phase=True)
    return res


def transmitted_purity(res, ws, wi, ws0, wi0, fw, target):
    """Gaussian filter on the narrow (heralding) arm transmitting `target`
    of the pairs: idler at theta_PM = 0, signal at theta_PM = 90."""
    from scipy.optimize import brentq
    tot = float(np.sum(np.abs(res.f) ** 2))
    hax = heralding_axis(res.f, ws, wi)
    h_axis, h0 = (ws, ws0) if hax == 0 else (wi, wi0)

    def trans(k):
        ff = apply_heralding_filter_gaussian(res.f, h_axis, h0, k * fw, hax)
        return float(np.sum(np.abs(ff) ** 2)) / tot

    k = brentq(lambda k: trans(k) - target, 0.01, 60.0)
    ff = apply_heralding_filter_gaussian(res.f, h_axis, h0, k * fw, hax)
    return schmidt_decompose(ff)[1], k * fw, h0


def nm_halfwidth(omega0, half):
    lam0 = 2 * np.pi * C_LIGHT / omega0
    return lam0 ** 2 * half / (2 * np.pi * C_LIGHT) * 1e9


def tail_rows():
    rows = []
    for pt in FEATURED:
        tag, plane, seq, _, pump_nm, lam_s, lam_i, theta_deg, _, L_mm, bw_nm, _ = pt
        lam_p = pump_nm / 1000.0
        wp0, ws0, wi0 = (lambda_um_to_omega(v) for v in (lam_p, lam_s, lam_i))
        fw = fwhm_omega_from_fwhm_wavelength(bw_nm * 1e-3, lam_p)
        sig = gaussian_pump_sigma(fw)
        row = {"tag": tag, "L_mm": L_mm, "pump_fwhm_nm": bw_nm,
               "pump_sigma_rad_ps": round(sig * 1e-12, 3)}
        Ps, Ws = [], []
        for w in W_SIGMAS:
            ws, wi, clipped = axes_for(ws0, wi0, w * sig, PIXEL_RAD_PS)
            P = purity_on(pt, L_mm, bw_nm, ws, wi).purity
            row[f"P_{int(w)}sig"] = round(P, 5)
            row[f"clipped_{int(w)}sig"] = int(clipped)
            if not clipped:
                Ps.append(P)
                Ws.append(w)
        if len(Ws) >= 3:
            A = np.vstack([np.ones(len(Ws)), 1.0 / np.array(Ws)]).T
            (p_inf, c), *_ = np.linalg.lstsq(A, np.array(Ps), rcond=None)
            resid = float(np.max(np.abs(A @ np.array([p_inf, c]) - np.array(Ps))))
            row.update({"P_inf": round(p_inf, 4), "C_per_sigma": round(c, 4),
                        "fit_max_residual": round(resid, 5),
                        "fit_points": len(Ws)})
        else:
            row.update({"P_inf": "", "C_per_sigma": "", "fit_max_residual": "",
                        "fit_points": len(Ws)})
        # the paper's own adaptive window, in the same units
        wsa, wia = adaptive_axes(np.radians(theta_deg), cfg_for(plane, seq),
                                 wp0, ws0, wi0, L_mm * 1e-3, bw_nm, lam_p, 200)
        hs = 0.5 * (wsa[-1] - wsa[0])
        hi = 0.5 * (wia[-1] - wia[0])
        row.update({"paper_halfwindow_s_sigma": round(hs / sig, 2),
                    "paper_halfwindow_i_sigma": round(hi / sig, 2),
                    "paper_halfwindow_s_nm": round(nm_halfwidth(ws0, hs), 1),
                    "paper_halfwindow_i_nm": round(nm_halfwidth(wi0, hi), 1)})
        # filtered purities are window independent: show at 6 and 16 sigma
        for w in (6.0, 16.0):
            ws, wi, _ = axes_for(ws0, wi0, w * sig, PIXEL_RAD_PS)
            res = purity_on(pt, L_mm, bw_nm, ws, wi)
            for target in (0.9, 0.7):
                P, fwf, h0 = transmitted_purity(res, ws, wi, ws0, wi0, fw, target)
                row[f"P_T{int(target*100)}_{int(w)}sig"] = round(P, 5)
                if w == 6.0:
                    row[f"filter_fwhm_T{int(target*100)}_nm"] = round(
                        2 * nm_halfwidth(h0, fwf / 2), 1)
        rows.append(row)
        print(tag, {k: v for k, v in row.items() if k.startswith(("P_", "C_", "paper"))})
    return rows


def length_rows():
    rows = []
    for pt in FEATURED[:3]:
        tag, plane, seq, _, pump_nm, lam_s, lam_i, theta_deg, _, _, _, _ = pt
        lam_p = pump_nm / 1000.0
        wp0, ws0, wi0 = (lambda_um_to_omega(v) for v in (lam_p, lam_s, lam_i))
        for L_mm in (10.0, 15.0, 20.0, 30.0):
            for bw in (2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0):
                fw = fwhm_omega_from_fwhm_wavelength(bw * 1e-3, lam_p)
                sig = gaussian_pump_sigma(fw)
                Ps, Ws = [], []
                for w in (6.0, 12.0, 24.0):
                    ws, wi, clipped = axes_for(ws0, wi0, w * sig, PIXEL_RAD_PS)
                    if clipped:
                        continue
                    Ps.append(purity_on(pt, L_mm, bw, ws, wi).purity)
                    Ws.append(w)
                if len(Ws) < 2:
                    continue
                A = np.vstack([np.ones(len(Ws)), 1.0 / np.array(Ws)]).T
                (p_inf, c), *_ = np.linalg.lstsq(A, np.array(Ps), rcond=None)
                rows.append({"tag": tag, "L_mm": L_mm, "pump_fwhm_nm": bw,
                             "P_6sig": round(Ps[0], 5), "P_inf": round(p_inf, 4),
                             "C_per_sigma": round(c, 4)})
                print(rows[-1])
    return rows


def write(name, rows):
    with open(os.path.join(OUT, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("tail", "all"):
        write("window_tail.csv", tail_rows())
    if what in ("length", "all"):
        write("length_extension.csv", length_rows())
