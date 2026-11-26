"""Numerical checks added for the first revision of the manuscript
(Optics Express 622047, November 2026; bibojsa v1.2.0).

Each stage answers one request of the reviewers and writes one or two CSVs
under survey/:

  resolution   survey/resolution_check.csv
               Folds and ends of the five type-II GVM curves, and the two
               weak x-z eoe features, retraced on a 0.1-nm pump grid at the
               survey's node counts (320 signal nodes, 1201 angles) and at
               1280 / 2401.  The survey itself uses a 2-nm pump grid.
  negatives    survey/negative_results_check.csv
               The four interactions the survey reports as never meeting a
               GVM condition, re-scanned at 1 nm / 640 signal nodes / 2401
               angles (survey field: 5 nm / 220 / 1201).  theta_PM is mapped
               to (-90, 90] deg.
  altfit       survey/xy_altfit_curves.csv, survey/xy_altfit_H.csv
               The two x-y type-II GVM curves traced with the Hellwig (2000)
               and Umemura (2007) Sellmeier fits, and the member pumped at
               464 nm (point H) with its purity on 6- and 12-sigma windows.
               The alternative fits are evaluated outside their stated
               validity ranges where the curve requires it; every row says
               whether its three wavelengths lie inside the fit's own range.
  window       survey/window_prescription.csv
               The adaptive collection window against the prescription
               R = 10 x (mean of the signal and idler central-peak FWHMs)
               of Cai, Baek, Jin & Kaneda, Adv. Phys. Res. 4, e00167 (2025),
               at A, B, H (asymmetric GVM) and D (EPM).
  tailfit      survey/window_tail_dense.csv, survey/window_tail_fits.csv
               P(W) on every integer window half-width the Sellmeier range
               allows, and the (P_inf, C) fits of P = P_inf + C/W over
               sub-ranges of W -- the sensitivity of the extrapolation to
               the fitted range, including the three-window fits behind
               points B and C in window_tail.csv.
  conditions   survey/window_tail_conditions.csv
               Tail coefficient C at matched crystal length: the EPM point D
               at 10 mm against its tabulated 2 mm.
  deff         survey/deff_benchmarks.csv
               |d_eff| at the three degenerate x-z benchmarks with the
               Kleinman-averaged coefficient (3.2 pm/V) and with the measured
               d26 = 3.48 pm/V, against Jin et al., Quantum Eng. 2, e38
               (2020).  The last columns reconstruct those values as
               d26 cos(theta - rho) times a Miller factor referred to
               second-harmonic generation of 1079.5 nm; that reconstruction
               is an inference, not a statement of the cited paper.

Run from bibojsa/:   python examples/revision_checks.py [stage ...]
(no argument = all stages; the tailfit stage takes about half an hour.)
"""
import csv
import dataclasses
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from landscape_headline import FEATURED
from landscape_points import W_MAX, W_MIN, adaptive_axes, cfg_for
from window_tail import PIXEL_RAD_PS, axes_for, purity_on

from bibojsa import landscape as L
from bibojsa.biaxial import principal_plane_index, walkoff_angle
from bibojsa.deff import deff_bibo_pm_per_v
from bibojsa.jsa import (build_jsa, fwhm_omega_from_fwhm_wavelength,
                         gaussian_pump_sigma)
from bibojsa.materials.bibo import BIBO, BIBO_HELLWIG, BIBO_UMEMURA
from bibojsa.phasematch import ThreeWaveConfig
from bibojsa.sellmeier import C_LIGHT, lambda_um_to_omega, refractive_index

X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "survey")
F = {f[0]: f for f in FEATURED}
FITS = {"Miyata-Umemura-Kato": BIBO, "Hellwig 2000": BIBO_HELLWIG,
        "Umemura 2007": BIBO_UMEMURA}


def write(name, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(OUT, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"  survey/{name}: {len(rows)} rows", flush=True)


# ------------------------------------------------------------- resolution

RESOLUTION_JOBS = [
    # label, plane, sequence, target, pump windows (um)
    ("x-z k_p'=k_e' fold", "xz", "oeo", 0, [(0.7180, 0.7240)]),
    ("x-z EPM fold", "xz", "oeo", 45, [(0.8780, 0.8840)]),
    ("x-z k_p'=k_e' short-pump end", "xz", "eoo", 90, [(0.4900, 0.4980)]),
    ("x-z EPM short-pump end", "xz", "eoo", 45, [(0.5520, 0.5600)]),
    ("x-z k_p'=k_o' short-pump end", "xz", "oeo", 90, [(0.8240, 0.8320)]),
    ("x-z k_p'=k_o' long-pump end", "xz", "eoo", 0, [(1.2160, 1.2260)]),
    ("x-y k_p'=k_o' both ends", "xy", "oee", 90,
     [(0.4340, 0.4420), (0.4880, 0.4960)]),
    ("x-y EPM both ends", "xy", "oee", 45,
     [(0.4800, 0.4880), (0.5480, 0.5560)]),
    ("x-z eoe theta_PM = 45 deg feature", "xz", "eoe", 45, [(0.7700, 0.7950)]),
    ("x-z eoe theta_PM = 0 both ends", "xz", "eoe", 0,
     [(0.6000, 0.6100), (0.6380, 0.6480)]),
]
RESOLUTION_STEP_UM = 0.0001


def stage_resolution():
    rows = []
    for label, plane, seq, target, wins in RESOLUTION_JOBS:
        pumps = np.concatenate([np.arange(a, b, RESOLUTION_STEP_UM)
                                for a, b in wins])
        for n_ls, n_th in ((320, 1201), (1280, 2401)):
            pts = []
            for branch in ("first", "second"):
                pts += L.trace_gvm_curve(plane, seq, target, X, Y, Z, pumps,
                                         n_ls=n_ls, branch=branch,
                                         n_theta=n_th)
            lp = np.array([p.lam_p_um for p in pts]) * 1e3
            ls = np.array([p.lam_s_um for p in pts])
            li = np.array([p.lam_i_um for p in pts])
            th = np.array([p.theta_deg for p in pts])
            rows.append({
                "feature": label, "plane": plane, "interaction": seq,
                "target_deg": target, "pump_step_nm": RESOLUTION_STEP_UM * 1e3,
                "n_signal_nodes": n_ls, "n_theta_nodes": n_th,
                "pump_scanned_nm": "; ".join(f"{a*1e3:.1f}-{b*1e3:.1f}"
                                             for a, b in wins),
                "n_points": len(pts),
                "pump_min_nm": round(lp.min(), 2), "pump_max_nm": round(lp.max(), 2),
                "lambda_s_min_um": round(ls.min(), 4), "lambda_s_max_um": round(ls.max(), 4),
                "lambda_i_min_um": round(li.min(), 4), "lambda_i_max_um": round(li.max(), 4),
                "angle_min_deg": round(th.min(), 3), "angle_max_deg": round(th.max(), 3)})
            print(f"  {label} [{n_ls}, {n_th}]: pump "
                  f"{lp.min():.1f}-{lp.max():.1f} nm ({len(pts)} points)",
                  flush=True)
    write("resolution_check.csv", rows)


# -------------------------------------------------------------- negatives

NEGATIVES = [("xy", "ooe"), ("xy", "eoe"), ("xz", "oee"), ("yz", "eeo")]


def stage_negatives():
    rows = []
    for n_ls, n_th, step in ((220, 1201, 0.005), (640, 2401, 0.001)):
        pumps = np.arange(0.34, 1.52 + 1e-9, step)
        for plane, seq in NEGATIVES:
            tested = matched = 0
            tmin, tmax, pm = np.inf, -np.inf, []
            for branch in ("first", "second"):
                for lp in pumps:
                    ls = np.linspace(max(lp * 1.0005, L.LAM_MIN),
                                     min(2 * lp * 0.999999, L.LAM_MAX), n_ls)
                    if ls[-1] <= ls[0]:
                        continue
                    tested += n_ls
                    th, li = L.solve_theta_grid(np.full_like(ls, lp), ls,
                                                plane, seq, X, Y, Z,
                                                n_theta=n_th, branch=branch)
                    ok = np.isfinite(th)
                    if not ok.any():
                        continue
                    g = L.gvm_functions(np.full(ok.sum(), lp), ls[ok], li[ok],
                                        th[ok], plane, seq, X, Y, Z)
                    t = ((g["theta_pm"] + 90.0) % 180.0) - 90.0
                    matched += int(ok.sum())
                    tmin, tmax = min(tmin, float(t.min())), max(tmax, float(t.max()))
                    pm.append(lp)
            rows.append({
                "plane": plane, "interaction": seq,
                "pump_range_um": "0.34-1.52", "pump_step_nm": step * 1e3,
                "n_signal_nodes": n_ls, "n_theta_nodes": n_th,
                "roots_per_node_searched": 2,
                "nodes_tested": tested, "phase_matched_nodes": matched,
                "pump_min_nm": round(min(pm) * 1e3, 1) if pm else "",
                "pump_max_nm": round(max(pm) * 1e3, 1) if pm else "",
                "theta_pm_min_deg": round(tmin, 2) if matched else "",
                "theta_pm_max_deg": round(tmax, 2) if matched else ""})
            print(f"  {plane}/{seq} [{step*1e3:.0f} nm, {n_ls}, {n_th}]: "
                  f"{matched} of {tested} nodes phase matched", flush=True)
    write("negative_results_check.csv", rows)


# ----------------------------------------------------------------- altfit

def _wide(fit):
    """The fit with its range check opened to the survey's 0.33-3.05 um, so
    that a curve can be followed across the edge of a fit's stated range."""
    return {a: dataclasses.replace(c, lambda_min=0.30, lambda_max=3.20)
            for a, c in fit.items()}


def stage_altfit():
    curves, hrows = [], []
    cfg = ThreeWaveConfig(plane="xy", pol_pump="e", pol_signal="e",
                          pol_idler="o")
    for name, fit in FITS.items():
        lo, hi = fit["x"].lambda_min, fit["x"].lambda_max
        c = _wide(fit)
        for target, cond in ((90, "k_p'=k_o'"), (45, "EPM")):
            pts = []
            for branch in ("first", "second"):
                pts += L.trace_gvm_curve("xy", "oee", target, c["x"], c["y"],
                                         c["z"], np.arange(0.40, 0.60 + 1e-9, 0.001),
                                         n_ls=320, branch=branch, n_theta=1201)
            for p in sorted(pts, key=lambda q: q.lam_p_um):
                d = abs(deff_bibo_pm_per_v("xy", "e", "e", "o",
                                           np.radians(p.theta_deg),
                                           enforce_range=False))
                curves.append({
                    "fit": name, "fit_range_um": f"{lo}-{hi}",
                    "condition": cond, "target_deg": target,
                    "pump_nm": round(p.lam_p_um * 1e3, 1),
                    "lambda_e_um": round(p.lam_s_um, 5),
                    "lambda_o_um": round(p.lam_i_um, 5),
                    "phi_deg": round(p.theta_deg, 3),
                    "abs_deff_pm_V": round(d, 4),
                    "gvm_scale_fs_mm": round(max(abs(p.kp - p.ks),
                                                 abs(p.kp - p.ki)) * 1e12, 2),
                    "inside_fit_range": int(p.lam_p_um >= lo and p.lam_i_um <= hi)})
            lp = [p.lam_p_um * 1e3 for p in pts]
            print(f"  {name}, {cond}: pump {min(lp):.0f}-{max(lp):.0f} nm "
                  f"({len(pts)} points)", flush=True)
        # the member pumped at 464 nm (point H) and its purity
        p = L.trace_gvm_curve("xy", "oee", 90, c["x"], c["y"], c["z"], [0.464],
                              n_ls=320, n_theta=1201)[0]
        wp0, ws0, wi0 = (lambda_um_to_omega(v)
                         for v in (0.464, p.lam_s_um, p.lam_i_um))
        fw = fwhm_omega_from_fwhm_wavelength(2e-3, 0.464)
        sig = gaussian_pump_sigma(fw)
        row = {"fit": name, "pump_nm": 464.0, "L_mm": 10.0, "pump_fwhm_nm": 2.0,
               "lambda_e_um": round(p.lam_s_um, 5),
               "lambda_o_um": round(p.lam_i_um, 5),
               "phi_deg": round(p.theta_deg, 3),
               "gvm_scale_fs_mm": round(abs(p.kp - p.ks) * 1e12, 2)}
        for w in (6, 12):
            n = int(2 * w * sig * 1e-12 / PIXEL_RAD_PS) + 1
            ws = np.linspace(ws0 - w * sig, ws0 + w * sig, n)
            wi = np.linspace(wi0 - w * sig, wi0 + w * sig, n)
            row[f"P_{w}sig"] = round(build_jsa(
                ws, wi, wp0, np.radians(p.theta_deg), cfg, c["x"], c["y"],
                c["z"], 10e-3, fw).purity, 5)
        hrows.append(row)
        print("  ", row, flush=True)
    write("xy_altfit_curves.csv", curves)
    write("xy_altfit_H.csv", hrows)


# ------------------------------------------------------- window utilities

def _setup(tag, L_mm=None, bw=None):
    _, plane, seq, _, pnm, ls, li, th, _, L0, bw0, _ = F[tag]
    L_mm, bw = L_mm or L0, bw or bw0
    lam_p = pnm / 1e3
    wp0, ws0, wi0 = (lambda_um_to_omega(v) for v in (lam_p, ls, li))
    fw = fwhm_omega_from_fwhm_wavelength(bw * 1e-3, lam_p)
    return dict(tag=tag, pt=F[tag], plane=plane, seq=seq, th=np.radians(th),
                lam_p=lam_p, wp0=wp0, ws0=ws0, wi0=wi0, fw=fw,
                sig=gaussian_pump_sigma(fw), L=L_mm, bw=bw)


def _clipped(s, half):
    return ((s["ws0"] - half < W_MIN) or (s["ws0"] + half > W_MAX)
            or (s["wi0"] - half < W_MIN) or (s["wi0"] + half > W_MAX))


def _purity_square(s, half, pixel=PIXEL_RAD_PS):
    ws, wi, _ = axes_for(s["ws0"], s["wi0"], half, pixel)
    return purity_on(s["pt"], s["L"], s["bw"], ws, wi)


def _scan(s, widths):
    Ws, Ps = [], []
    for w in widths:
        if _clipped(s, w * s["sig"]):
            continue
        Ws.append(w)
        Ps.append(_purity_square(s, w * s["sig"]).purity)
    return Ws, Ps


def _fit(Ws, Ps):
    A = np.vstack([np.ones(len(Ws)), 1.0 / np.array(Ws, dtype=float)]).T
    (p, c), *_ = np.linalg.lstsq(A, np.array(Ps), rcond=None)
    res = float(np.max(np.abs(A @ np.array([p, c]) - np.array(Ps))))
    return float(p), float(c), res


def _fwhm(x, y):
    """FWHM of the central peak of y(x), linearly interpolated."""
    y = y / y.max()
    k = int(np.argmax(y))
    lft = k
    while lft > 0 and y[lft - 1] >= 0.5:
        lft -= 1
    rgt = k
    while rgt < len(y) - 1 and y[rgt + 1] >= 0.5:
        rgt += 1
    xl = x[lft] if lft == 0 else np.interp(0.5, [y[lft - 1], y[lft]],
                                           [x[lft - 1], x[lft]])
    xr = x[rgt] if rgt == len(y) - 1 else np.interp(0.5, [y[rgt + 1], y[rgt]],
                                                    [x[rgt + 1], x[rgt]])
    return abs(xr - xl)


def _nm(omega0, half):
    lam0 = 2 * np.pi * C_LIGHT / omega0
    return lam0 ** 2 * half / (2 * np.pi * C_LIGHT) * 1e9


# ----------------------------------------------------------------- window

def stage_window():
    rows = []
    for tag in ("A", "B", "H", "D"):
        s = _setup(tag)
        ws, wi, _ = axes_for(s["ws0"], s["wi0"], 12 * s["sig"], 0.1)
        inten = np.abs(purity_on(s["pt"], s["L"], s["bw"], ws, wi).f) ** 2
        ks = int(np.argmin(np.abs(ws - s["ws0"])))
        ki = int(np.argmin(np.abs(wi - s["wi0"])))
        widths = {
            "cuts through the peak": (_fwhm(ws, inten[:, ki]), _fwhm(wi, inten[ks, :])),
            "marginal spectra": (_fwhm(ws, inten.sum(1)), _fwhm(wi, inten.sum(0)))}
        wsa, wia = adaptive_axes(s["th"], cfg_for(s["plane"], s["seq"]),
                                 s["wp0"], s["ws0"], s["wi0"], s["L"] * 1e-3,
                                 s["bw"], s["lam_p"], 200)
        p_ad = build_jsa(wsa, wia, s["wp0"], s["th"],
                         cfg_for(s["plane"], s["seq"]), X, Y, Z,
                         s["L"] * 1e-3, s["fw"]).purity
        for name, (a, b) in widths.items():
            dw = 0.5 * (a + b)
            half = 5.0 * dw                      # R = 10 dw is the full width
            p = _purity_square(s, half, pixel=min(PIXEL_RAD_PS,
                                                  dw * 1e-12 / 20)).purity
            rows.append({
                "tag": tag, "L_mm": s["L"], "pump_fwhm_nm": s["bw"],
                "fwhm_definition": name,
                "fwhm_signal_rad_ps": round(a * 1e-12, 2),
                "fwhm_idler_rad_ps": round(b * 1e-12, 2),
                "R_full_width_rad_ps": round(2 * half * 1e-12, 1),
                "halfwidth_pump_sigma": round(half / s["sig"], 2),
                "halfwidth_signal_nm": round(_nm(s["ws0"], half), 1),
                "halfwidth_idler_nm": round(_nm(s["wi0"], half), 1),
                "P_prescription": round(p, 5),
                "adaptive_halfwidth_pump_sigma": round(
                    0.5 * (wsa[-1] - wsa[0]) / s["sig"], 2),
                "P_adaptive": round(p_ad, 5)})
            print("  ", rows[-1], flush=True)
    write("window_prescription.csv", rows)


# ---------------------------------------------------------------- tailfit

DENSE = {  # tag: window half-widths in pump sigma (clipped ones are skipped)
    "A": (5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 20, 24, 32),
    "G": (5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 20, 24),
    "H": (5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 20, 24, 32),
    "D": (5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 20, 24),
    "B": (4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16),
    "C": (4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16),
}
FIT_RANGES = ((6, 8), (6, 12), (8, 12), (5, 12), (4, 12), (9, 16), (6, 16),
              (6, 24), (12, 32), (6, 32))


def stage_tailfit():
    dense, fits = [], []
    for tag, widths in DENSE.items():
        t0 = time.time()
        s = _setup(tag)
        Ws, Ps = _scan(s, widths)
        for w, p in zip(Ws, Ps):
            dense.append({"tag": tag, "L_mm": s["L"], "pump_fwhm_nm": s["bw"],
                          "W_pump_sigma": w, "purity": round(p, 5)})
        sets = [("{6, 8, 12}", [6, 8, 12])]
        sets += [(f"{lo}-{hi}", [w for w in Ws if lo <= w <= hi])
                 for lo, hi in FIT_RANGES]
        seen = set()
        for label, sel in sets:
            sel = [w for w in sel if w in Ws]
            if len(sel) < 3 or tuple(sel) in seen:
                continue
            seen.add(tuple(sel))
            p, c, res = _fit(sel, [Ps[Ws.index(w)] for w in sel])
            fits.append({"tag": tag, "windows_fitted_sigma": label,
                         "n_windows": len(sel), "W_min": min(sel),
                         "W_max": max(sel), "P_inf": round(p, 4),
                         "C_per_sigma": round(c, 4),
                         "fit_max_residual": round(res, 5)})
        print(f"  {tag}: {len(Ws)} windows up to {max(Ws)} sigma "
              f"({time.time()-t0:.0f} s)", flush=True)
    write("window_tail_dense.csv", dense)
    write("window_tail_fits.csv", fits)


# ------------------------------------------------------------- conditions

CONDITIONS = [("D", 2.0, 12.0), ("D", 10.0, 2.4), ("D", 10.0, 2.0),
              ("D", 10.0, 3.0)]


def stage_conditions():
    rows = []
    for tag, L_mm, bw in CONDITIONS:
        s = _setup(tag, L_mm, bw)
        Ws, Ps = _scan(s, (6, 8, 12, 16, 24))
        p, c, res = _fit(Ws, Ps)
        rows.append({"tag": tag, "L_mm": L_mm, "pump_fwhm_nm": bw,
                     "n_windows": len(Ws), "P_6sig": round(Ps[0], 5),
                     "P_inf": round(p, 4), "C_per_sigma": round(c, 4),
                     "fit_max_residual": round(res, 5)})
        print("  ", rows[-1], flush=True)
    write("window_tail_conditions.csv", rows)


# ------------------------------------------------------------------- deff

JIN2020 = {"A": 2.50, "D": 2.48, "C": 2.30}   # Quantum Eng. 2, e38, Table III
D26_MEASURED = 3.48                           # pm/V (Hellwig 2000 via Ghotbi 2004)
D_KLEINMAN = 3.2
MILLER_REFERENCE_UM = 1.0795                  # fundamental of the Maker-fringe SHG


def _chi(n):
    return n * n - 1.0


def stage_deff():
    rows = []
    for tag, jin in JIN2020.items():
        _, _, _, _, pnm, ls, _, th, *_ = F[tag]
        lam_p, t = pnm / 1e3, np.radians(th)
        rho = abs(walkoff_angle(ls, t, "xz", X, Y, Z))

        def miller_product(lp, ld):
            return (_chi(refractive_index(lp, Y)) * _chi(refractive_index(ld, Y))
                    * _chi(principal_plane_index(ld, t, "xz", X, Y, Z)[1]))

        m = (miller_product(lam_p, ls)
             / miller_product(MILLER_REFERENCE_UM / 2, MILLER_REFERENCE_UM))
        ours, meas = D_KLEINMAN * np.cos(t), D26_MEASURED * np.cos(t)
        rows.append({
            "tag": tag, "pump_nm": pnm, "theta_deg": th,
            "deff_kleinman_3p2_pm_V": round(ours, 3),
            "deff_d26_3p48_pm_V": round(meas, 3),
            "deff_jin2020_pm_V": jin,
            "jin_minus_kleinman_percent": round(100 * (jin - ours) / ours, 1),
            "jin_minus_d26_percent": round(100 * (jin - meas) / meas, 1),
            "rho_mrad": round(rho * 1e3, 1),
            "miller_factor": round(m, 4),
            "reconstruction_d26_cos_theta_minus_rho_miller_pm_V": round(
                D26_MEASURED * np.cos(t - rho) * m, 3)})
        print("  ", rows[-1], flush=True)
    write("deff_benchmarks.csv", rows)


STAGES = {"deff": stage_deff, "altfit": stage_altfit,
          "resolution": stage_resolution, "negatives": stage_negatives,
          "window": stage_window, "conditions": stage_conditions,
          "tailfit": stage_tailfit}

if __name__ == "__main__":
    todo = sys.argv[1:] or list(STAGES)
    for name in todo:
        print(name, flush=True)
        STAGES[name]()
