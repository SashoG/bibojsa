"""BIBO exploration survey (project brief, exploration targets).

For pump wavelengths 405/515/775/780/800/810 nm, scan all phase-matchable
principal-plane interaction types of BIBO (LPR2010 Table 1, both relevant
octants: theta in (0, 180 deg) for the y-z and x-z planes, phi in
(0, 90 deg) for x-y), solving collinear Delta k = 0 at each signal
wavelength; at every phase-matched point record theta_PM (JSA ridge
orientation), signed deff and |deff| (LPR2010 Table 1), spatial walk-off
of every e-polarized wave, and the three group delays.

All outputs are plain CSV (OriginLab-ready), no figures:

  survey/survey_loci.csv        every phase-matched point on every locus
  survey/survey_crossings.csv   refined theta_PM = 0/45/90 deg crossings

Convention: lambda1 >= lambda2 > lambda3 (idler, signal, pump); o/e
sequence strings are idler-signal-pump, matching bibojsa.deff.
theta_PM = arctan2(-(kp'-ks'), (kp'-ki')) per the brief; crossings sought:
  0  deg : kp' = ks' (asymmetric, signal)
  90 deg : kp' = ki' (asymmetric, idler)
  +45 deg: kp' = (ks'+ki')/2 (symmetric / EPM)
(-45/135 deg is the trivial anticorrelated orientation and is not a target.)
"""
import csv
import json
import os
import sys

import numpy as np
from scipy.optimize import brentq

from bibojsa.materials.bibo import BIBO
from bibojsa.phasematch import ThreeWaveConfig, delta_k
from bibojsa.gvm import group_index_at_angle
from bibojsa.deff import deff_bibo_pm_per_v
from bibojsa.biaxial import walkoff_angle
from bibojsa.sellmeier import lambda_um_to_omega, optic_axis_angle_deg, C_LIGHT

X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]
PUMPS_UM = [0.405, 0.515, 0.775, 0.780, 0.800, 0.810]
LAM_MIN, LAM_MAX = 0.33, 3.05  # Sellmeier validity with small margin

# (plane, seq) with seq = idler-signal-pump o/e labels (LPR2010 Table 1)
INTERACTIONS = [
    ("xy", "ooe"), ("xy", "eoe"),
    ("yz", "eeo"), ("yz", "oeo"),
    ("xz", "eeo"), ("xz", "oeo"), ("xz", "ooe"), ("xz", "eoe"),
]

N_LS = 41          # signal-wavelength grid points per locus
N_ANGLE = 1441     # angle-grid resolution for root bracketing


def angle_range(plane):
    if plane == "xy":
        return np.radians(0.2), np.radians(89.8)
    return np.radians(0.2), np.radians(179.8)


def config_for(plane, seq):
    return ThreeWaveConfig(plane=plane, pol_pump=seq[2], pol_signal=seq[1],
                            pol_idler=seq[0])


def pm_roots(lam_p, lam_s, lam_i, plane, seq):
    """All angles in the allowed range with Delta k = 0 (vectorized scan +
    brentq refinement)."""
    cfg = config_for(plane, seq)
    wp, ws, wi = (lambda_um_to_omega(l) for l in (lam_p, lam_s, lam_i))
    lo, hi = angle_range(plane)
    grid = np.linspace(lo, hi, N_ANGLE)
    vals = delta_k(grid, wp, ws, wi, cfg, X, Y, Z)
    s = np.sign(vals)
    idx = np.where(np.diff(s) != 0)[0]
    roots = []
    for i in idx:
        try:
            r = brentq(lambda th: delta_k(th, wp, ws, wi, cfg, X, Y, Z),
                       grid[i], grid[i + 1], xtol=1e-10)
            roots.append(r)
        except ValueError:
            pass
    return roots


def kprimes(lam_p, lam_s, lam_i, theta, plane, seq):
    wp, ws, wi = (lambda_um_to_omega(l) for l in (lam_p, lam_s, lam_i))
    kp = group_index_at_angle(wp, theta, plane, seq[2], X, Y, Z)
    ks = group_index_at_angle(ws, theta, plane, seq[1], X, Y, Z)
    ki = group_index_at_angle(wi, theta, plane, seq[0], X, Y, Z)
    return kp, ks, ki


def theta_pm_from_kprimes(kp, ks, ki):
    return np.degrees(np.arctan2(-(kp - ks), (kp - ki)))


def walkoffs_mrad(lam_p, lam_s, lam_i, theta, plane, seq):
    out = []
    for lam, pol in ((lam_p, seq[2]), (lam_s, seq[1]), (lam_i, seq[0])):
        if pol == "e":
            out.append(1e3 * walkoff_angle(lam, theta, plane, X, Y, Z))
        else:
            out.append(0.0)
    return out  # pump, signal, idler in mrad


def row_at(lam_p, lam_s, plane, seq, theta):
    lam_i = 1.0 / (1.0 / lam_p - 1.0 / lam_s)
    kp, ks, ki = kprimes(lam_p, lam_s, lam_i, theta, plane, seq)
    tpm = theta_pm_from_kprimes(kp, ks, ki)
    d = deff_bibo_pm_per_v(plane, seq[2], seq[1], seq[0], theta,
                            lambda_pump_um=lam_p, enforce_range=False)
    rp, rs, ri = walkoffs_mrad(lam_p, lam_s, lam_i, theta, plane, seq)
    region_ok = 1
    if plane == "xz":
        om = optic_axis_angle_deg(lam_p, X, Y, Z)
        inside = om < np.degrees(theta) < 180.0 - om
        expected = seq in ("eeo", "oeo")
        region_ok = int(inside == expected)
    gd = lambda k: k * C_LIGHT / 0.299792458  # ps/mm
    return {
        "pump_nm": round(lam_p * 1000, 1), "plane": plane, "interaction": seq,
        "lambda_s_um": lam_s, "lambda_i_um": lam_i,
        "angle_deg": np.degrees(theta), "theta_pm_deg": tpm,
        "deff_pm_V": d, "abs_deff_pm_V": abs(d), "region_ok": region_ok,
        "rho_pump_mrad": rp, "rho_signal_mrad": rs, "rho_idler_mrad": ri,
        "GD_p_ps_mm": gd(kp), "GD_s_ps_mm": gd(ks), "GD_i_ps_mm": gd(ki),
    }


def lam_s_grid(lam_p):
    lam_s_min = max(1.0 / (1.0 / lam_p - 1.0 / LAM_MAX), lam_p * 1.02, LAM_MIN)
    lam_s_max = 2.0 * lam_p * 0.9995
    if lam_s_min >= lam_s_max:
        return np.array([])
    return np.linspace(lam_s_min, lam_s_max, N_LS)


FIELDS = ["pump_nm", "plane", "interaction", "lambda_s_um", "lambda_i_um",
          "angle_deg", "theta_pm_deg", "deff_pm_V", "abs_deff_pm_V",
          "region_ok", "rho_pump_mrad", "rho_signal_mrad", "rho_idler_mrad",
          "GD_p_ps_mm", "GD_s_ps_mm", "GD_i_ps_mm"]


def stage_loci():
    done = set()
    path = "survey/survey_loci.csv"
    if os.path.exists(path):
        with open(path) as fh:
            for r in csv.DictReader(fh):
                done.add((r["pump_nm"], r["plane"], r["interaction"]))
    new_file = not os.path.exists(path)
    fh = open(path, "a", newline="")
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new_file:
        w.writeheader()
    for lam_p in PUMPS_UM:
        for plane, seq in INTERACTIONS:
            key = (str(round(lam_p * 1000, 1)), plane, seq)
            if key in done:
                continue
            grid = lam_s_grid(lam_p)
            nrows = 0
            for lam_s in grid:
                lam_i = 1.0 / (1.0 / lam_p - 1.0 / lam_s)
                if not (LAM_MIN <= lam_i <= LAM_MAX):
                    continue
                for th in pm_roots(lam_p, lam_s, lam_i, plane, seq):
                    w.writerow(row_at(lam_p, lam_s, plane, seq, th))
                    nrows += 1
            fh.flush()
            print(f"loci {lam_p*1000:.0f}nm {plane}/{seq}: {nrows} rows", flush=True)
    fh.close()
    print("loci done", flush=True)


def stage_crossings():
    import collections
    rows = list(csv.DictReader(open("survey/survey_loci.csv")))
    for r in rows:
        for k in ("lambda_s_um", "angle_deg", "theta_pm_deg", "abs_deff_pm_V"):
            r[k] = float(r[k])
    groups = collections.defaultdict(list)
    for r in rows:
        groups[(r["pump_nm"], r["plane"], r["interaction"])].append(r)

    out = open("survey/survey_crossings.csv", "w", newline="")
    w = csv.DictWriter(out, fieldnames=["point_id", "crossing_type_deg"] + FIELDS)
    w.writeheader()
    pid = 0

    def gfuncs(lam_p, plane, seq, lam_s, theta):
        lam_i = 1.0 / (1.0 / lam_p - 1.0 / lam_s)
        kp, ks, ki = kprimes(lam_p, lam_s, lam_i, theta, plane, seq)
        return {0: kp - ks, 90: kp - ki, 45: kp - 0.5 * (ks + ki)}

    for (pump_nm, plane, seq), rs in sorted(groups.items()):
        lam_p = float(pump_nm) / 1000.0
        rs.sort(key=lambda r: (r["lambda_s_um"], r["angle_deg"]))
        # split into branches: consecutive lambda_s, root index by angle order
        by_ls = collections.defaultdict(list)
        for r in rs:
            by_ls[r["lambda_s_um"]].append(r)
        ls_sorted = sorted(by_ls)
        branches = collections.defaultdict(list)
        for ls in ls_sorted:
            pts = sorted(by_ls[ls], key=lambda r: r["angle_deg"])
            for j, r in enumerate(pts):
                branches[(len(pts), j)].append(r)
        for bkey, br in branches.items():
            if len(br) < 2:
                continue
            for target in (0, 90, 45):
                vals = []
                for r in br:
                    g = gfuncs(lam_p, plane, seq, r["lambda_s_um"], np.radians(r["angle_deg"]))
                    vals.append(g[target])
                for a, b, va, vb in zip(br, br[1:], vals, vals[1:]):
                    if va * vb >= 0 or va == 0:
                        continue
                    lo, hi = a["lambda_s_um"], b["lambda_s_um"]

                    def f(ls):
                        roots = pm_roots(lam_p, ls, 1.0/(1.0/lam_p-1.0/ls), plane, seq)
                        if not roots:
                            return np.nan
                        th = min(roots, key=lambda t: abs(np.degrees(t) - 0.5*(a["angle_deg"]+b["angle_deg"])))
                        return gfuncs(lam_p, plane, seq, ls, th)[target]

                    try:
                        ls_star = brentq(f, lo, hi, xtol=1e-6)
                    except ValueError:
                        continue
                    roots = pm_roots(lam_p, ls_star, 1.0/(1.0/lam_p-1.0/ls_star), plane, seq)
                    th_star = min(roots, key=lambda t: abs(np.degrees(t) - 0.5*(a["angle_deg"]+b["angle_deg"])))
                    row = row_at(lam_p, ls_star, plane, seq, th_star)
                    pid += 1
                    w.writerow(dict({"point_id": pid, "crossing_type_deg": target}, **row))
                    print(f"crossing #{pid}: {pump_nm}nm {plane}/{seq} "
                          f"type={target} lam_s={ls_star:.4f} th={row['angle_deg']:.2f} "
                          f"|deff|={row['abs_deff_pm_V']:.2f}", flush=True)
    out.close()
    print("crossings done", flush=True)


if __name__ == "__main__":
    os.makedirs("survey", exist_ok=True)
    stage = sys.argv[1] if len(sys.argv) > 1 else "loci"
    if stage == "loci":
        stage_loci()
    elif stage == "crossings":
        stage_crossings()
