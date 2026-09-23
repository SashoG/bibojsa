"""Per-operating-point data sets for the BIBO survey (CSV only, no figures).

Selects unique operating points from survey/survey_crossings.csv
(|deff| >= 0.1 pm/V; x-z mirror octant duplicates removed -- theta and
180-theta are equivalent in the x-z plane, |deff| and dispersion both
symmetric; y-z octants would be kept separately, but all y-z crossings
found have deff = 0 exactly and are excluded).

For each point <id> writes survey/points/point<id>_*.csv:
  _info.csv       one-row header: pump, plane, type, wavelengths, angle,
                  theta_PM, deff, walk-offs, GDs
  _purity.csv     rows: L_mm, pump_fwhm_nm, purity_unfiltered, K,
                  purity_filt1x, purity_filt3x, edge
                  (filtered = Gaussian heralding filter on the idler arm,
                  intensity-FWHM = 1x / 3x the pump's intensity FWHM in
                  angular frequency)
  _schmidt.csv    first 30 Schmidt eigenvalues at the best unfiltered
                  (L, bw) of the grid
  _jsa_map.csv    |f| (normalized) matrix at the best (L, bw); first row =
                  idler wavelength axis (um), first column = signal
                  wavelength axis (um)
  _brightness.csv waist_um, rho_max_mrad, L_eff_mm, and the brightness
                  proxy deff^2 * min(L, L_eff)^2 [pm^2/V^2 mm^2] for
                  L = 0.5, 1, 2, 5, 10 mm
"""
import csv
import json
import os
import sys

import numpy as np

from bibojsa.materials.bibo import BIBO
from bibojsa.phasematch import ThreeWaveConfig, delta_k
from bibojsa.jsa import (build_jsa, fwhm_omega_from_fwhm_wavelength,
                          apply_heralding_filter_gaussian, schmidt_decompose)
from bibojsa.sellmeier import lambda_um_to_omega, omega_to_lambda_um

X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]
W_MAX = lambda_um_to_omega(0.33)
W_MIN = lambda_um_to_omega(3.05)

L_GRID_MM = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
BW_GRID_NM = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
PROXY_L_MM = [0.5, 1.0, 2.0, 5.0, 10.0]
WAISTS_UM = [20, 40, 60, 80, 100, 120, 140, 160, 180, 200]
N_GRID = 140


def load_points():
    rows = list(csv.DictReader(open("survey/survey_crossings.csv")))
    pts, seen = [], set()
    for r in rows:
        if float(r["abs_deff_pm_V"]) < 0.1:
            continue
        th = float(r["angle_deg"])
        th_canon = min(th, 180.0 - th) if r["plane"] == "xz" else th
        key = (r["pump_nm"], r["plane"], r["interaction"],
               r["crossing_type_deg"], round(float(r["lambda_s_um"]), 4),
               round(th_canon, 2))
        if key in seen:
            continue
        seen.add(key)
        if r["plane"] == "xz" and th > 90.0:
            continue  # mirror duplicate
        pts.append(r)
    return pts


def adaptive_axes(theta, cfg, wp0, ws0, wi0, L_m, fwhm_nm, lam_p, n):
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
    ws_lo = max(ws0 - 1.15 * em - pad, W_MIN)
    ws_hi = min(ws0 + 1.15 * ep + pad, W_MAX)
    wi_lo = max(wi0 - 1.15 * ep - pad, W_MIN)
    wi_hi = min(wi0 + 1.15 * em + pad, W_MAX)
    return np.linspace(ws_lo, ws_hi, n), np.linspace(wi_lo, wi_hi, n)


def run_case(pt, L_mm, fwhm_nm, n=N_GRID):
    plane, seq = pt["plane"], pt["interaction"]
    cfg = ThreeWaveConfig(plane=plane, pol_pump=seq[2], pol_signal=seq[1],
                           pol_idler=seq[0])
    lam_p = float(pt["pump_nm"]) / 1000.0
    theta = np.radians(float(pt["angle_deg"]))
    wp0 = lambda_um_to_omega(lam_p)
    ws0 = lambda_um_to_omega(float(pt["lambda_s_um"]))
    wi0 = lambda_um_to_omega(float(pt["lambda_i_um"]))
    L_m = L_mm * 1e-3
    fw = fwhm_omega_from_fwhm_wavelength(fwhm_nm * 1e-3, lam_p)
    ws_axis, wi_axis = adaptive_axes(theta, cfg, wp0, ws0, wi0, L_m, fwhm_nm, lam_p, n)
    res = build_jsa(ws_axis, wi_axis, wp0, theta, cfg, X, Y, Z, L_m, fw,
                    pump_shape="gaussian", include_phase=True)
    A = np.abs(res.f)
    m = A.max()
    edge = max(A[0].max(), A[-1].max(), A[:, 0].max(), A[:, -1].max()) / m if m else 0.0
    # heralding filter on the idler arm (axis 1), centered at wi0
    p1 = schmidt_decompose(apply_heralding_filter_gaussian(res.f, wi_axis, wi0, 1.0 * fw, 1))[1]
    p3 = schmidt_decompose(apply_heralding_filter_gaussian(res.f, wi_axis, wi0, 3.0 * fw, 1))[1]
    return res, edge, p1, p3


def main():
    os.makedirs("survey/points", exist_ok=True)
    pts = load_points()
    print(f"{len(pts)} unique operating points", flush=True)
    for pt in pts:
        pid = pt["point_id"]
        base = f"survey/points/point{pid}"
        if os.path.exists(base + "_brightness.csv"):
            continue  # done
        # info
        with open(base + "_info.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(pt.keys()))
            w.writeheader(); w.writerow(pt)
        # purity table (resumable per point via temp jsonl)
        cache = base + "_purity.jsonl"
        done = {}
        if os.path.exists(cache):
            for line in open(cache):
                r = json.loads(line)
                done[(r["L_mm"], r["fwhm_nm"])] = r
        for L in L_GRID_MM:
            for bw in BW_GRID_NM:
                if (L, bw) in done:
                    continue
                res, edge, p1, p3 = run_case(pt, L, bw)
                r = {"L_mm": L, "fwhm_nm": bw, "purity_unfiltered": res.purity,
                     "K": res.schmidt_number, "purity_filt1x": p1,
                     "purity_filt3x": p3, "edge": edge}
                done[(L, bw)] = r
                with open(cache, "a") as fh:
                    fh.write(json.dumps(r) + "\n")
        with open(base + "_purity.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["L_mm", "fwhm_nm",
                "purity_unfiltered", "K", "purity_filt1x", "purity_filt3x", "edge"])
            w.writeheader()
            for L in L_GRID_MM:
                for bw in BW_GRID_NM:
                    w.writerow(done[(L, bw)])
        # best unfiltered case -> schmidt + jsa map
        best = max(done.values(), key=lambda r: r["purity_unfiltered"])
        res, edge, _, _ = run_case(pt, best["L_mm"], best["fwhm_nm"])
        with open(base + "_schmidt.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["n", "lambda_n"])
            w.writerow(["# L_mm", best["L_mm"], "fwhm_nm", best["fwhm_nm"],
                        "purity", f"{res.purity:.5f}"])
            for i, lam in enumerate(res.schmidt_lambdas[:30], 1):
                w.writerow([i, f"{lam:.6e}"])
        A = np.abs(res.f); A = A / A.max()
        ls = omega_to_lambda_um(res.omega_s)
        li = omega_to_lambda_um(res.omega_i)
        with open(base + "_jsa_map.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["lambda_s_um\\lambda_i_um"] + [f"{v:.6f}" for v in li])
            for irow in range(A.shape[0]):
                w.writerow([f"{ls[irow]:.6f}"] + [f"{A[irow, j]:.5f}" for j in range(A.shape[1])])
        # walk-off / brightness proxy
        rho_max = max(abs(float(pt["rho_pump_mrad"])),
                      abs(float(pt["rho_signal_mrad"])),
                      abs(float(pt["rho_idler_mrad"])))  # mrad
        d2 = float(pt["abs_deff_pm_V"]) ** 2
        with open(base + "_brightness.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["waist_um", "rho_max_mrad", "L_eff_mm"]
                       + [f"proxy_L{L:g}mm_pm2V-2mm2" for L in PROXY_L_MM])
            for w0 in WAISTS_UM:
                L_eff_mm = np.inf if rho_max == 0 else (w0 * 1e-3) / (rho_max * 1e-3)
                row = [w0, f"{rho_max:.3f}",
                       "inf" if np.isinf(L_eff_mm) else f"{L_eff_mm:.3f}"]
                for L in PROXY_L_MM:
                    row.append(f"{d2 * min(L, L_eff_mm) ** 2:.4f}")
                w.writerow(row)
        print(f"point {pid} done: {pt['pump_nm']}nm {pt['plane']}/{pt['interaction']} "
              f"type={pt['crossing_type_deg']} bestP={best['purity_unfiltered']:.3f} "
              f"(L={best['L_mm']}mm bw={best['fwhm_nm']}nm) "
              f"filt3x={best['purity_filt3x']:.3f}", flush=True)
    print("points done", flush=True)


if __name__ == "__main__":
    main()
