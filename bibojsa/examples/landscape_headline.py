"""Manuscript-grade tables for the featured operating points.

Consumes the traced landscape (`landscape_survey.py`) and produces the CSVs
the paper's tables and figures are drawn from:

  survey/headline_points.csv      one row per featured point (Table 2; A-F, R
                                  and, since 17 Oct 2026, G, H, J, K, M)
  survey/headline_filter.csv      purity and transmitted heralding flux vs
                                  filter width, for the featured points
  survey/headline_convergence.csv purity vs (window scale, grid size) at each
                                  featured point's quoted configuration
  survey/family_purity.csv        best converged purity along each GVM family,
                                  as a function of pump wavelength
  survey/headline_jsa_<tag>.csv   |f| in angular-frequency detuning

Every purity here is reported together with the spectral window it was
computed on (through the boundary-amplitude "edge" metric) and, when a
heralding filter is applied, with the fraction of pairs that filter
transmits.
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from landscape_points import (FILTERS, adaptive_axes, brightness_rows,
                              cfg_for, run_case)

from bibojsa import landscape as L
from bibojsa.jsa import (apply_heralding_filter_gaussian, build_jsa,
                         heralding_axis,
                         fwhm_omega_from_fwhm_wavelength, schmidt_decompose)
from bibojsa.materials.bibo import BIBO
from bibojsa.sellmeier import lambda_um_to_omega

X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]

# tag, plane, seq, target, pump_nm, lam_s, lam_i, theta_deg, |deff|, L_mm, bw_nm, note
FEATURED = [
    ("A", "xz", "oeo",  0,  682.241, 1.364482, 1.364482, 45.96827, 2.224, 10.0,  8.0,
     "asymmetric GVM at frequency degeneracy"),
    ("B", "xz", "oeo",  0,  722.000, 1.155982, 1.923690, 41.35980, 2.404, 10.0, 12.0,
     "asymmetric GVM, Ti:sapphire-accessible member"),
    ("C", "xz", "oeo", 90, 1116.841, 2.233682, 2.233682, 46.10995, 2.218, 10.0, 16.0,
     "asymmetric GVM at frequency degeneracy, long-pump family"),
    ("D", "xz", "oeo", 45,  869.583, 1.739166, 1.739166, 43.77035, 2.311,  2.0, 12.0,
     "extended phase matching at frequency degeneracy"),
    ("E", "xz", "ooe",  0,  818.696, 1.637392, 1.637392, 11.03604, 3.141, 10.0, 32.0,
     "triple-GVM degenerate point, largest deff"),
    ("F", "xz", "eeo",  0,  896.760, 1.793520, 1.793520, 34.90800, 1.558, 10.0, 32.0,
     "triple-GVM degenerate point, eeo"),
    ("R", "xz", "oeo", 45,  810.000, 1.138800, 2.805300, 39.42211, 2.472,  5.0,  8.0,
     "nondegenerate EPM member (idler in the infrared absorption band)"),
    # --- second half of the type-II surfaces (17 Oct 2026); seq is still
    # idler-signal-pump with lambda_i >= lambda_s, so "eoo" = e daughter longer,
    # "oee" (x-y) = z-polarized o daughter longer.  Cut angle = theta (x-z) or
    # phi (x-y).  Selected by examples/landscape_headline_full.py (now folded
    # in here): (L, bandwidth) = best unfiltered purity on the adaptive window.
    ("G", "xz", "eoo", 90,  532.000, 0.773630, 1.703290, 67.74340, 1.212, 10.0,  4.0,
     "continuation of A's k_p'=k_e' curve, e daughter longer; 532-nm pump"),
    ("H", "xy", "oee", 90,  464.000, 0.622290, 1.824160, 44.91930, 1.660, 10.0,  2.0,
     "x-y plane asymmetric GVM (pump matched to the z-polarized wave) at the d14 maximum"),
    ("J", "xy", "oee", 45,  532.000, 0.673500, 2.532160, 31.47900, 1.479,  2.0,  2.0,
     "x-y plane EPM member, 532-nm pump"),
    ("K", "xz", "eoo", 45,  800.000, 1.335830, 1.994420, 48.22960, 2.132,  2.0,  8.0,
     "x-z EPM member with both daughters below 2.6 um, 800-nm pump"),
    ("M", "xz", "eoo",  0, 1150.000, 2.168230, 2.448820, 48.08110, 2.138, 10.0, 24.0,
     "continuation of C's k_p'=k_o' curve, e daughter longer"),
]

TRANSPARENCY_HI_UM = 2.60   # 1 cm^-1 absorption onset (Ghotbi 2006 OE; LPR2010 Table 4), not a transmission cutoff
# Log-spaced so that the transmitted-flux axis is covered from a hard
# filter (a few per cent of the heralding flux) to effectively none.
FILTER_RATIOS = np.geomspace(0.01, 24.0, 46)


def as_pt(f):
    tag, plane, seq, target, pnm, ls, li, th, d, Lmm, bw, note = f
    return {"plane": plane, "interaction": seq, "target_deg": target,
            "pump_nm": pnm, "lambda_s_um": ls, "lambda_i_um": li,
            "angle_deg": th, "abs_deff_pm_V": d}


def filter_curve(f):
    """Purity and transmitted heralding flux vs Gaussian filter width."""
    tag, plane, seq, target, pnm, ls, li, th, d, Lmm, bw, note = f
    pt = as_pt(f)
    cfg = cfg_for(plane, seq)
    lam_p = pnm / 1000.0
    theta = np.radians(th)
    wp0 = lambda_um_to_omega(lam_p)
    ws0 = lambda_um_to_omega(ls)
    wi0 = lambda_um_to_omega(li)
    L_m = Lmm * 1e-3
    fw = fwhm_omega_from_fwhm_wavelength(bw * 1e-3, lam_p)
    wsa, wia = adaptive_axes(theta, cfg, wp0, ws0, wi0, L_m, bw, lam_p, 200)
    res = build_jsa(wsa, wia, wp0, theta, cfg, X, Y, Z, L_m, fw,
                    pump_shape="gaussian", include_phase=True)
    tot = float(np.sum(np.abs(res.f) ** 2))
    # heralding filter on the narrow arm: idler (axis 1) at theta_PM = 0,
    # signal (axis 0) at theta_PM = 90; symmetric points keep the idler
    hax = heralding_axis(res.f, wsa, wia)
    h_axis, h0, lam_h = (wsa, ws0, ls) if hax == 0 else (wia, wi0, li)
    rows = []
    for k in FILTER_RATIOS:
        ff = apply_heralding_filter_gaussian(res.f, h_axis, h0, k * fw, hax)
        p = schmidt_decompose(ff)[1]
        t = float(np.sum(np.abs(ff) ** 2)) / tot
        # filter intensity FWHM expressed as a wavelength width on the heralding arm
        dlam_nm = (k * fw) * (lam_h * 1e-6) ** 2 / (2 * np.pi * 2.99792458e8) * 1e9
        rows.append({"tag": tag, "herald_arm": "signal" if hax == 0 else "idler",
                     "ratio": round(float(k), 4),
                     "filter_fwhm_nm": round(dlam_nm, 3),
                     "transmission": round(t, 5), "purity": round(p, 5)})
    rows.append({"tag": tag, "herald_arm": "signal" if hax == 0 else "idler",
                 "ratio": "", "filter_fwhm_nm": "",
                 "transmission": 1.0, "purity": round(res.purity, 5)})
    return rows, res


def convergence(f):
    tag = f[0]
    pt = as_pt(f)
    out = []
    for scale in (0.7, 0.85, 1.0, 1.15, 1.3):
        for n in (120, 180, 260, 340):
            res, edge, filt, *_ = run_case(pt, f[9], f[10], n=n, scale=scale)
            out.append({"tag": tag, "window_scale": scale, "N": n,
                        "purity": round(res.purity, 6),
                        "purity_filt1x": round(filt[1.0][0], 6),
                        "trans_filt1x": round(filt[1.0][1], 5),
                        "edge": round(edge, 5)})
    return out


def family_purity():
    """Best converged unfiltered purity along each GVM family vs pump."""
    rows = list(csv.DictReader(open("survey/landscape_curves.csv")))
    import collections
    fam = collections.defaultdict(list)
    for r in rows:
        if float(r["abs_deff_pm_V"]) < 0.5:
            continue
        fam[(r["plane"], r["interaction"], r["target_deg"])].append(r)
    out = []
    for key, rs in sorted(fam.items()):
        rs.sort(key=lambda r: float(r["pump_nm"]))
        idx = np.unique(np.linspace(0, len(rs) - 1, min(len(rs), 18)).astype(int))
        for i in idx:
            r = rs[i]
            pt = {"plane": r["plane"], "interaction": r["interaction"],
                  "pump_nm": r["pump_nm"], "lambda_s_um": r["lambda_s_um"],
                  "lambda_i_um": r["lambda_i_um"], "angle_deg": r["angle_deg"],
                  "abs_deff_pm_V": r["abs_deff_pm_V"]}
            best = (0.0, None, None, None)
            for Lmm in (1.0, 2.0, 5.0, 10.0):
                for bw in (2.0, 4.0, 8.0, 12.0, 16.0, 24.0):
                    try:
                        res, edge, filt, *_ = run_case(pt, Lmm, bw, n=120)
                    except Exception:
                        continue
                    if edge < 0.15 and res.purity > best[0]:
                        best = (res.purity, Lmm, bw, (edge, filt[1.0]))
            if best[1] is None:
                continue
            out.append({
                "plane": r["plane"], "interaction": r["interaction"],
                "target_deg": r["target_deg"], "pump_nm": float(r["pump_nm"]),
                "lambda_s_um": float(r["lambda_s_um"]),
                "lambda_i_um": float(r["lambda_i_um"]),
                "angle_deg": float(r["angle_deg"]),
                "abs_deff_pm_V": float(r["abs_deff_pm_V"]),
                "gvm_scale_fs_mm": float(r["gvm_scale_fs_mm"]),
                "degeneracy": float(r["degeneracy"]),
                "lambda_i_in_transparency": int(r["lambda_i_in_transparency"]),
                "best_purity": round(best[0], 5), "best_L_mm": best[1],
                "best_fwhm_nm": best[2], "edge": round(best[3][0], 4),
                "purity_filt1x": round(best[3][1][0], 5),
                "trans_filt1x": round(best[3][1][1], 4)})
            print(f"  {r['plane']}/{r['interaction']} t={r['target_deg']:>2} "
                  f"{float(r['pump_nm']):7.1f} nm: P={best[0]:.4f} "
                  f"(L={best[1]}, {best[2]} nm)", flush=True)
    return out


def main():
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"

    if stage in ("points", "all"):
        prows, frows, crows = [], [], []
        for f in FEATURED:
            tag, plane, seq, target, pnm, ls, li, th, d, Lmm, bw, note = f
            pt = as_pt(f)
            res, edge, filt, wsa, wia, wp0, ws0, wi0 = run_case(pt, Lmm, bw, n=200)
            br, rho = brightness_rows(pt)
            leff60 = [b["L_eff_mm"] for b in br if b["waist_um"] == 60][0]
            gf = L.gvm_functions(pnm / 1000.0, ls, li, np.radians(th), plane,
                                 seq, X, Y, Z)
            fc, _ = filter_curve(f)
            frows += fc
            crows += convergence(f)
            # purity at the filter width that keeps 90 % / 70 % of the flux
            tt = np.array([r["transmission"] for r in fc if r["ratio"] != ""])
            pp = np.array([r["purity"] for r in fc if r["ratio"] != ""])
            o = np.argsort(tt)
            def at(t):
                if t < tt[o][0] or t > tt[o][-1]:
                    return float("nan")
                return float(np.interp(t, tt[o], pp[o]))
            p90, p70 = at(0.90), at(0.70)
            prows.append({
                "tag": tag, "plane": plane, "interaction": seq,
                "target_deg": target, "pump_nm": pnm,
                "lambda_s_um": round(ls, 5), "lambda_i_um": round(li, 5),
                "angle_deg": round(th, 4), "abs_deff_pm_V": d,
                "d32_fs_mm": round(float(gf["g0"]) * 1e12, 3),
                "d31_fs_mm": round(float(gf["g90"]) * 1e12, 3),
                "theta_pm_deg": round(float(gf["theta_pm"]), 3),
                "L_mm": Lmm, "pump_fwhm_nm": bw,
                "purity": round(res.purity, 5), "edge": round(edge, 5),
                "herald_arm": fc[0]["herald_arm"],
                "purity_T90": round(p90, 5), "purity_T70": round(p70, 5),
                "rho_max_mrad": round(rho * 1e3, 2),
                "L_eff_60um_mm": leff60,
                "lambda_i_in_transparency": int(li <= TRANSPARENCY_HI_UM),
                "note": note})
            A = np.abs(res.f); A = A / A.max()
            with open(f"survey/headline_jsa_{tag}.csv", "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["dws_rad_ps", "dwi_rad_ps", "amp"])
                step = max(1, len(wsa) // 140)
                for a in range(0, len(wsa), step):
                    for b in range(0, len(wia), step):
                        w.writerow([round((wsa[a] - ws0) * 1e-12, 6),
                                    round((wia[b] - wi0) * 1e-12, 6),
                                    round(float(A[a, b]), 6)])
            print(f"{tag}: P={res.purity:.4f} edge={edge:.4f} "
                  f"P(T=90%)={p90:.4f} P(T=70%)={p70:.4f} rho={rho*1e3:.1f} mrad",
                  flush=True)
        for name, rows, keys in [
            ("survey/headline_points.csv", prows, list(prows[0])),
            ("survey/headline_filter.csv", frows, list(frows[0])),
            ("survey/headline_convergence.csv", crows, list(crows[0]))]:
            with open(name, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=keys)
                w.writeheader()
                w.writerows(rows)

    if stage in ("family", "all"):
        rows = family_purity()
        with open("survey/family_purity.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
