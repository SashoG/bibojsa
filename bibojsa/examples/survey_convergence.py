"""Purity convergence tables for the top-3 survey operating points
(per the brief's requirement: convergence check on purity vs grid size
AND window, since higher-order curvature of the phase-matching band can
matter). Milestone-5 style, but two-dimensional: grid size N x window
scale s, where s multiplies the adaptively-probed window half-widths
before clipping to the Sellmeier validity range (0.32-3.09 um).

Points and their best low-truncation (L, pump-FWHM) cases
(survey_summary.csv):
  #45  810 nm  xz/oeo t45   L = 5 mm,  bw = 8 nm
  #31  800 nm  xz/oeo t45   L = 10 mm, bw = 4 nm
  #42  810 nm  xz/eeo t45   L = 10 mm, bw = 16 nm

Interpretation guide (recorded in each CSV):
- Convergence in N at fixed s: numerical grid-resolution convergence.
- Convergence in s: window-size convergence. For the oeo points the
  idler window is clipped by the validity limit, so purity may keep
  drifting with s and/or the edge metric stays high: that quantifies the
  documented truncation systematic rather than removing it.
- purity_filt3x (Gaussian idler-arm heralding filter, intensity FWHM =
  3x the pump's in angular frequency) is the window-robust quantity: the
  filter suppresses the slowly-decaying sinc tails (and the clipped
  region near 3.05 um), so it converges in BOTH N and s.

RESULTS (2026-07-13 run, this Sellmeier fit):
- N-convergence: |dP| < 8e-4 from N=80 to N=300 at every window, all
  three points. Grid resolution is a non-issue.
- Unfiltered purity drifts slowly downward with window scale (sinc-tail
  Schmidt weight): #45: 0.819->0.808, #31: 0.804->0.788,
  #42: 0.806->0.770 over s = 0.7->1.3 (still decreasing at 1.3).
  Unfiltered numbers must therefore be quoted with a defined window (or
  treated as upper bounds).
- Filtered purity is converged to < 2e-3 across the whole s range and
  flat to < 5e-5 for s >= 1.15:
      #45: P_filt = 0.8696;  #31: 0.8622;  #42: 0.8540.
  These are the paper-grade numbers.

Output: survey/convergence/point<id>_convergence.csv
        survey/convergence/convergence_all.csv
"""
import csv
import os

import numpy as np

from bibojsa.materials.bibo import BIBO
from bibojsa.phasematch import ThreeWaveConfig, delta_k
from bibojsa.jsa import (build_jsa, fwhm_omega_from_fwhm_wavelength,
                          apply_heralding_filter_gaussian, schmidt_decompose)
from bibojsa.sellmeier import lambda_um_to_omega

X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]
W_MAX = lambda_um_to_omega(0.33)
W_MIN = lambda_um_to_omega(3.05)

CASES = [  # point_id, pump_um, plane, seq, lambda_s, lambda_i, theta_deg, L_mm, bw_nm
    (45, 0.810, "xz", "oeo", 1.1388, 2.8053, 39.422, 5.0, 8.0),
    (31, 0.800, "xz", "oeo", 1.1090, 2.8710, 39.302, 10.0, 4.0),
    (42, 0.810, "xz", "eeo", 1.1735, 2.6151, 35.234, 10.0, 16.0),
]
N_GRID = [80, 120, 160, 220, 300]
SCALES = [0.7, 0.85, 1.0, 1.15, 1.3]


def axes_for(theta, cfg, wp0, ws0, wi0, L_m, fwhm_nm, lam_p, n, scale):
    def first_ext(deltas):
        ws, wi = ws0 + deltas, wi0 - deltas
        valid = (ws <= W_MAX) & (ws >= W_MIN) & (wi <= W_MAX) & (wi >= W_MIN)
        x = np.where(valid,
                     delta_k(theta, np.where(valid, ws + wi, wp0),
                             np.where(valid, ws, ws0), np.where(valid, wi, wi0),
                             cfg, X, Y, Z) * L_m / 2.0, 0.0)
        dead = (~valid) | (np.abs(np.sinc(x / np.pi)) < 0.05)
        return abs(deltas[np.argmax(dead)]) if dead.any() else abs(deltas[-1])

    probe = np.linspace(0, 6e14, 400)
    ep, em = first_ext(probe), first_ext(-probe)
    sigma = fwhm_omega_from_fwhm_wavelength(fwhm_nm * 1e-3, lam_p)
    pad = 3.0 * sigma
    ws_lo = max(ws0 - scale * (1.15 * em + pad), W_MIN)
    ws_hi = min(ws0 + scale * (1.15 * ep + pad), W_MAX)
    wi_lo = max(wi0 - scale * (1.15 * ep + pad), W_MIN)
    wi_hi = min(wi0 + scale * (1.15 * em + pad), W_MAX)
    return np.linspace(ws_lo, ws_hi, n), np.linspace(wi_lo, wi_hi, n)


def main():
    os.makedirs("survey/convergence", exist_ok=True)
    allrows = []
    for pid, lam_p, plane, seq, lam_s, lam_i, th_deg, L_mm, bw in CASES:
        cfg = ThreeWaveConfig(plane=plane, pol_pump=seq[2], pol_signal=seq[1],
                               pol_idler=seq[0])
        theta = np.radians(th_deg)
        wp0, ws0, wi0 = (lambda_um_to_omega(l) for l in (lam_p, lam_s, lam_i))
        L_m = L_mm * 1e-3
        fw = fwhm_omega_from_fwhm_wavelength(bw * 1e-3, lam_p)
        rows = []
        for s in SCALES:
            for n in N_GRID:
                ws_axis, wi_axis = axes_for(theta, cfg, wp0, ws0, wi0, L_m,
                                             bw, lam_p, n, s)
                res = build_jsa(ws_axis, wi_axis, wp0, theta, cfg, X, Y, Z,
                                L_m, fw, pump_shape="gaussian", include_phase=True)
                A = np.abs(res.f); m = A.max()
                edge = max(A[0].max(), A[-1].max(), A[:, 0].max(),
                           A[:, -1].max()) / m if m else 0.0
                pf = schmidt_decompose(apply_heralding_filter_gaussian(
                    res.f, wi_axis, wi0, 3.0 * fw, 1))[1]
                r = {"point_id": pid, "L_mm": L_mm, "pump_fwhm_nm": bw,
                     "window_scale": s, "N": n,
                     "purity": f"{res.purity:.6f}", "K": f"{res.schmidt_number:.4f}",
                     "purity_filt3x": f"{pf:.6f}",
                     "edge": f"{edge:.4f}",
                     "ws_lo_um": f"{2*np.pi*2.99792458e8/ws_axis[-1]*1e6:.4f}",
                     "ws_hi_um": f"{2*np.pi*2.99792458e8/ws_axis[0]*1e6:.4f}",
                     "wi_lo_um": f"{2*np.pi*2.99792458e8/wi_axis[-1]*1e6:.4f}",
                     "wi_hi_um": f"{2*np.pi*2.99792458e8/wi_axis[0]*1e6:.4f}"}
                rows.append(r); allrows.append(r)
            print(f"#{pid} s={s}: P(N) = " +
                  " ".join(f"{float(r['purity']):.5f}" for r in rows[-len(N_GRID):]),
                  f"edge={rows[-1]['edge']}", flush=True)
        with open(f"survey/convergence/point{pid}_convergence.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
    with open("survey/convergence/convergence_all.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(allrows[0].keys()))
        w.writeheader(); w.writerows(allrows)
    print("convergence done", flush=True)


if __name__ == "__main__":
    main()
