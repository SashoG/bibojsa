"""Joint-spectral analysis of the operating points selected from the
continuous-pump landscape (`landscape_survey.py`).

For every selected point this writes, under survey/points2/:

  point<id>_info.csv        one-row header (interaction, wavelengths, angle,
                            theta_PM, deff, group-delay differences, walk-off)
  point<id>_purity.csv      L_mm, pump_fwhm_nm, purity, K, edge, and for each
                            heralding-filter width the purity *and* the
                            fraction of pairs the filter transmits
  point<id>_convergence.csv purity vs (window scale, grid size)
  point<id>_jsa.csv         |f| on the best gently-filtered configuration,
                            in angular-frequency detuning (rad/s) so that the
                            45 deg ridge orientation is visible as 45 deg
  point<id>_schmidt.csv     Schmidt spectrum
  survey/points2_summary.csv  one row per point

Filter widths are quoted as multiples of the pump's own intensity FWHM in
angular frequency; the transmitted fraction is reported alongside every
filtered purity, because a purity bought with a filter is only meaningful
together with the flux it costs.
"""
import csv
import os
import sys

import numpy as np

from bibojsa import landscape as L
from bibojsa.biaxial import walkoff_angle
from bibojsa.deff import deff_bibo_pm_per_v
from bibojsa.jsa import (apply_heralding_filter_gaussian, build_jsa,
                         fwhm_omega_from_fwhm_wavelength, schmidt_decompose)
from bibojsa.materials.bibo import BIBO
from bibojsa.phasematch import ThreeWaveConfig, delta_k
from bibojsa.sellmeier import lambda_um_to_omega

X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]
W_MAX = lambda_um_to_omega(L.LAM_MIN)
W_MIN = lambda_um_to_omega(L.LAM_MAX)

L_GRID_MM = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
BW_GRID_NM = [1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 32.0]
FILTERS = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
N_GRID = 160
EDGE_MAX = 0.15
WAISTS_UM = [20, 40, 60, 80, 100, 140, 200, 300, 420]


def cfg_for(plane, seq):
    return ThreeWaveConfig(plane=plane, pol_pump=seq[2], pol_signal=seq[1],
                           pol_idler=seq[0])


def adaptive_axes(theta, cfg, wp0, ws0, wi0, L_m, fwhm_nm, lam_p, n, scale=1.0):
    """Grid windows from the actual decay of |f| along the antidiagonal,
    clipped to the Sellmeier validity range (same recipe as the original
    survey, with an explicit window-scale factor for convergence tests)."""
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
    ws_lo = max(ws0 - scale * 1.15 * em - pad, W_MIN)
    ws_hi = min(ws0 + scale * 1.15 * ep + pad, W_MAX)
    wi_lo = max(wi0 - scale * 1.15 * ep - pad, W_MIN)
    wi_hi = min(wi0 + scale * 1.15 * em + pad, W_MAX)
    return np.linspace(ws_lo, ws_hi, n), np.linspace(wi_lo, wi_hi, n)


def run_case(pt, L_mm, fwhm_nm, n=N_GRID, scale=1.0):
    plane, seq = pt["plane"], pt["interaction"]
    cfg = cfg_for(plane, seq)
    lam_p = float(pt["pump_nm"]) / 1000.0
    theta = np.radians(float(pt["angle_deg"]))
    wp0 = lambda_um_to_omega(lam_p)
    ws0 = lambda_um_to_omega(float(pt["lambda_s_um"]))
    wi0 = lambda_um_to_omega(float(pt["lambda_i_um"]))
    L_m = L_mm * 1e-3
    fw = fwhm_omega_from_fwhm_wavelength(fwhm_nm * 1e-3, lam_p)
    wsa, wia = adaptive_axes(theta, cfg, wp0, ws0, wi0, L_m, fwhm_nm, lam_p,
                             n, scale=scale)
    res = build_jsa(wsa, wia, wp0, theta, cfg, X, Y, Z, L_m, fw,
                    pump_shape="gaussian", include_phase=True)
    A = np.abs(res.f)
    m = A.max()
    edge = (max(A[0].max(), A[-1].max(), A[:, 0].max(), A[:, -1].max()) / m
            if m else 0.0)
    tot = float(np.sum(A ** 2))
    filt = {}
    for k in FILTERS:
        ff = apply_heralding_filter_gaussian(res.f, wia, wi0, k * fw, 1)
        filt[k] = (schmidt_decompose(ff)[1],
                   float(np.sum(np.abs(ff) ** 2)) / tot if tot else 0.0)
    return res, edge, filt, wsa, wia, wp0, ws0, wi0


def brightness_rows(pt):
    """Walk-off-limited effective length and brightness proxy vs waist.

    L_eff is the aperture length l_a = sqrt(pi) * w0 / rho (Boyd's
    convention); the proxy is deff^2 * min(L, L_eff)^2 in pm^2/V^2 mm^2.
    """
    plane, seq = pt["plane"], pt["interaction"]
    theta = np.radians(float(pt["angle_deg"]))
    rho = 0.0
    for lam, pol in ((float(pt["pump_nm"]) / 1000.0, seq[2]),
                     (float(pt["lambda_s_um"]), seq[1]),
                     (float(pt["lambda_i_um"]), seq[0])):
        if pol == "e":
            rho = max(rho, abs(walkoff_angle(lam, theta, plane, X, Y, Z)))
    d = abs(float(pt["abs_deff_pm_V"]))
    rows = []
    for w0 in WAISTS_UM:
        leff = np.sqrt(np.pi) * (w0 * 1e-6) / rho * 1e3 if rho > 0 else np.inf
        for Lm in [0.5, 1.0, 2.0, 5.0, 10.0]:
            rows.append({"waist_um": w0, "rho_max_mrad": round(rho * 1e3, 3),
                         "L_eff_mm": round(leff, 4) if np.isfinite(leff) else "",
                         "L_mm": Lm,
                         "brightness_proxy": round(d ** 2 * min(Lm, leff) ** 2, 5)})
    return rows, rho


def analyse(pt, pid, outdir):
    rows = []
    best_unf = (0.0, None)
    best_gentle = (0.0, None)
    for Lmm in L_GRID_MM:
        for bw in BW_GRID_NM:
            res, edge, filt, *_ = run_case(pt, Lmm, bw)
            row = {"L_mm": Lmm, "pump_fwhm_nm": bw,
                   "purity": round(res.purity, 6),
                   "K": round(res.schmidt_number, 6),
                   "edge": round(edge, 5)}
            for k, (p, t) in filt.items():
                row[f"purity_filt{k:g}x"] = round(p, 6)
                row[f"trans_filt{k:g}x"] = round(t, 6)
            rows.append(row)
            if edge < EDGE_MAX:
                if res.purity > best_unf[0]:
                    best_unf = (res.purity, row)
                if filt[1.0][0] > best_gentle[0]:
                    best_gentle = (filt[1.0][0], row)

    fields = (["L_mm", "pump_fwhm_nm", "purity", "K", "edge"]
              + [f"{a}_filt{k:g}x" for k in FILTERS for a in ("purity", "trans")])
    with open(f"{outdir}/point{pid}_purity.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    # convergence at the best gently-filtered configuration
    conv = []
    if best_gentle[1]:
        Lb, bwb = best_gentle[1]["L_mm"], best_gentle[1]["pump_fwhm_nm"]
        for scale in (0.7, 0.85, 1.0, 1.15, 1.3):
            for n in (100, 140, 200, 300):
                res, edge, filt, *_ = run_case(pt, Lb, bwb, n=n, scale=scale)
                conv.append({"window_scale": scale, "N": n,
                             "purity": round(res.purity, 6),
                             "purity_filt1x": round(filt[1.0][0], 6),
                             "purity_filt3x": round(filt[3.0][0], 6),
                             "edge": round(edge, 5)})
        with open(f"{outdir}/point{pid}_convergence.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["window_scale", "N", "purity",
                                               "purity_filt1x",
                                               "purity_filt3x", "edge"])
            w.writeheader()
            w.writerows(conv)

        # JSA map in angular-frequency detuning, plus the Schmidt spectrum
        res, edge, filt, wsa, wia, wp0, ws0, wi0 = run_case(pt, Lb, bwb, n=140)
        A = np.abs(res.f)
        A = A / A.max()
        with open(f"{outdir}/point{pid}_jsa.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["dws_rad_ps", "dwi_rad_ps", "amp"])
            for a, ws in enumerate(wsa):
                for b, wi in enumerate(wia):
                    w.writerow([round((ws - ws0) * 1e-12, 6),
                                round((wi - wi0) * 1e-12, 6),
                                round(float(A[a, b]), 6)])
        with open(f"{outdir}/point{pid}_schmidt.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["n", "lambda_n"])
            for i, lam in enumerate(res.schmidt_lambdas[:30]):
                w.writerow([i + 1, round(float(lam), 8)])

    br, rho = brightness_rows(pt)
    with open(f"{outdir}/point{pid}_brightness.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["waist_um", "rho_max_mrad",
                                           "L_eff_mm", "L_mm",
                                           "brightness_proxy"])
        w.writeheader()
        w.writerows(br)

    summ = dict(pt)
    summ["point_id"] = pid
    summ["rho_max_mrad"] = round(rho * 1e3, 3)
    if best_unf[1]:
        summ["best_unfilt_purity"] = best_unf[1]["purity"]
        summ["best_unfilt_L_mm"] = best_unf[1]["L_mm"]
        summ["best_unfilt_fwhm_nm"] = best_unf[1]["pump_fwhm_nm"]
        summ["best_unfilt_edge"] = best_unf[1]["edge"]
    if best_gentle[1]:
        g = best_gentle[1]
        summ["best_gentle_purity"] = g["purity_filt1x"]
        summ["best_gentle_trans"] = g["trans_filt1x"]
        summ["best_gentle_unfilt"] = g["purity"]
        summ["best_gentle_L_mm"] = g["L_mm"]
        summ["best_gentle_fwhm_nm"] = g["pump_fwhm_nm"]
        summ["best_gentle_purity_2x"] = g["purity_filt2x"]
        summ["best_gentle_trans_2x"] = g["trans_filt2x"]
        summ["best_gentle_purity_3x"] = g["purity_filt3x"]
        summ["best_gentle_trans_3x"] = g["trans_filt3x"]
        summ["best_gentle_edge"] = g["edge"]
    return summ


def select_points():
    """Operating points to analyse, drawn from the traced landscape."""
    pts = []
    seen = set()

    def add(row, tag):
        key = (row["plane"], row["interaction"], row["target_deg"],
               round(float(row["pump_nm"]), 1))
        if key in seen or float(row["abs_deff_pm_V"]) < 0.1:
            return
        seen.add(key)
        r = dict(row)
        r["selection"] = tag
        pts.append(r)

    with open("survey/landscape_degenerate.csv") as fh:
        for r in csv.DictReader(fh):
            add(r, "degenerate")

    with open("survey/landscape_curves.csv") as fh:
        curves = list(csv.DictReader(fh))
    import collections
    fam = collections.defaultdict(list)
    for r in curves:
        if float(r["abs_deff_pm_V"]) < 0.1:
            continue
        fam[(r["plane"], r["interaction"], r["target_deg"])].append(r)
    for key, rs in fam.items():
        rs.sort(key=lambda r: float(r["pump_nm"]))
        inside = [r for r in rs if r["lambda_i_in_transparency"] == "1"]
        # the member closest to degeneracy inside the transparency window,
        # the most nondegenerate one inside it, and the family's endpoints
        if inside:
            add(min(inside, key=lambda r: float(r["degeneracy"])), "in-window")
            add(max(inside, key=lambda r: float(r["degeneracy"])),
                "in-window-nondeg")
        add(rs[0], "endpoint-short-pump")
        add(rs[-1], "endpoint-long-pump")
    return pts


if __name__ == "__main__":
    outdir = "survey/points2"
    os.makedirs(outdir, exist_ok=True)
    pts = select_points()
    only = sys.argv[1] if len(sys.argv) > 1 else None
    print(f"{len(pts)} operating points selected", flush=True)
    summaries = []
    for i, pt in enumerate(pts, start=1):
        pid = i
        if only and only != str(pid):
            continue
        s = analyse(pt, pid, outdir)
        summaries.append(s)
        print(f"#{pid} {s['plane']}/{s['interaction']} t={s['target_deg']} "
              f"{s['pump_nm']}nm ({s['selection']}): "
              f"P={s.get('best_unfilt_purity')} "
              f"P1x={s.get('best_gentle_purity')} "
              f"(T={s.get('best_gentle_trans')})", flush=True)
    if summaries:
        keys = []
        for s in summaries:
            for k in s:
                if k not in keys:
                    keys.append(k)
        with open("survey/points2_summary.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(summaries)
