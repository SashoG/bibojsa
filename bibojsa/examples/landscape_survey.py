"""Continuous-pump GVM landscape survey of BIBO.

Supersedes `bibo_survey.py`, which sampled the phase-matching surface at six
fixed pump wavelengths.  Here the pump wavelength is a *survey variable*:
for every principal-plane interaction of LPR2010 Table 1 the whole
(lambda_p, lambda_s) phase-matching surface is covered, and the three
group-velocity-matching conditions

    theta_PM =  0 deg  <->  k_p' = k_s'          (asymmetric, signal)
    theta_PM = 45 deg  <->  k_p' = (k_s'+k_i')/2 (symmetric / EPM)
    theta_PM = 90 deg  <->  k_p' = k_i'          (asymmetric, idler)

are traced as continuous level curves, together with the pump wavelength at
which each curve crosses frequency degeneracy (lambda_s = lambda_i = 2 lambda_p).

Outputs (all CSV, all under survey/):

  landscape_curves.csv      every point on every traced GVM level curve
  landscape_degenerate.csv  the frequency-degenerate member of each curve
  landscape_field.csv       theta_PM sampled on a coarse (lambda_p, lambda_s)
                            grid, for the landscape contour figure
  landscape_coverage.csv    per-interaction record of what was found,
                            including the negative results

Run:  python examples/landscape_survey.py [field|curves|degenerate|all]
"""
import csv
import os
import sys
import time

import numpy as np

from bibojsa import landscape as L
from bibojsa.biaxial import walkoff_angle
from bibojsa.deff import deff_bibo_pm_per_v
from bibojsa.materials.bibo import BIBO

X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]

# (plane, seq) with seq = idler-signal-pump o/e labels (LPR2010 Table 1)
INTERACTIONS = [
    ("xy", "ooe"), ("xy", "eoe"),
    ("yz", "eeo"), ("yz", "oeo"),
    ("xz", "eeo"), ("xz", "oeo"), ("xz", "ooe"), ("xz", "eoe"),
]
TARGETS = (0, 45, 90)
BRANCHES = ("first", "second")

# Pump range: the short end is the Sellmeier floor, the long end is where
# even the degenerate pair reaches the 3.05 um ceiling (2 * 1.525 um).
PUMP_LO, PUMP_HI = 0.34, 1.52
PUMP_STEP_TRACE = 0.002      # 2 nm, curve tracing
PUMP_STEP_FIELD = 0.005      # 5 nm, contour field
N_LS_TRACE = 320
N_LS_FIELD = 220
N_THETA = 1201

CURVE_FIELDS = [
    "plane", "interaction", "target_deg", "branch", "pump_nm",
    "lambda_s_um", "lambda_i_um", "angle_deg", "theta_pm_deg",
    "degeneracy", "deff_pm_V", "abs_deff_pm_V",
    "d31_fs_mm", "d32_fs_mm", "gvm_scale_fs_mm",
    "rho_pump_mrad", "rho_signal_mrad", "rho_idler_mrad",
    "GD_p_ps_mm", "GD_s_ps_mm", "GD_i_ps_mm",
    "lambda_i_in_transparency",
]

# BIBO's measured transmission window (Hellwig/Liebertz/Bohaty 2000).  The
# Sellmeier fit is quoted to 3.09 um, which is a *fit* range, not a
# transparency range; every point is flagged against the physical window.
TRANSPARENCY_HI_UM = 2.60   # 1 cm^-1 absorption onset (Ghotbi 2006 OE; LPR2010 Table 4), not a transmission cutoff


def enrich(p: L.GVMPoint) -> dict:
    """Full record for one traced GVM point."""
    th = np.radians(p.theta_deg)
    seq, plane = p.seq, p.plane
    d = deff_bibo_pm_per_v(plane, seq[2], seq[1], seq[0], th,
                           lambda_pump_um=p.lam_p_um, enforce_range=False)
    rho = []
    for lam, pol in ((p.lam_p_um, seq[2]), (p.lam_s_um, seq[1]),
                     (p.lam_i_um, seq[0])):
        rho.append(1e3 * walkoff_angle(lam, th, plane, X, Y, Z)
                   if pol == "e" else 0.0)
    d31 = (p.kp - p.ki) * 1e15 * 1e-3      # fs/mm
    d32 = (p.kp - p.ks) * 1e15 * 1e-3
    return {
        "plane": plane, "interaction": seq, "target_deg": p.target,
        "branch": p.branch, "pump_nm": round(p.lam_p_um * 1000, 3),
        "lambda_s_um": round(p.lam_s_um, 6),
        "lambda_i_um": round(p.lam_i_um, 6),
        "angle_deg": round(p.theta_deg, 5),
        "theta_pm_deg": round(p.theta_pm_deg, 5),
        "degeneracy": round(p.degeneracy, 6),
        "deff_pm_V": round(d, 5), "abs_deff_pm_V": round(abs(d), 5),
        "d31_fs_mm": round(d31, 4), "d32_fs_mm": round(d32, 4),
        "gvm_scale_fs_mm": round(max(abs(d31), abs(d32)), 4),
        "rho_pump_mrad": round(rho[0], 3),
        "rho_signal_mrad": round(rho[1], 3),
        "rho_idler_mrad": round(rho[2], 3),
        "GD_p_ps_mm": round(float(L.gd_ps_per_mm(p.kp)), 6),
        "GD_s_ps_mm": round(float(L.gd_ps_per_mm(p.ks)), 6),
        "GD_i_ps_mm": round(float(L.gd_ps_per_mm(p.ki)), 6),
        "lambda_i_in_transparency": int(p.lam_i_um <= TRANSPARENCY_HI_UM),
    }


def stage_field():
    """theta_PM sampled over the phase-matching surface, for the contour map."""
    pumps = np.arange(PUMP_LO, PUMP_HI + 1e-9, PUMP_STEP_FIELD)
    with open("survey/landscape_field.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["plane", "interaction", "branch", "pump_nm",
                    "lambda_s_um", "lambda_i_um", "angle_deg",
                    "theta_pm_deg", "abs_deff_pm_V"])
        for plane, seq in INTERACTIONS:
            for branch in BRANCHES:
                n = 0
                t0 = time.time()
                for lam_p in pumps:
                    ls = np.linspace(max(lam_p * 1.0005, L.LAM_MIN),
                                     min(2.0 * lam_p * 0.999999, L.LAM_MAX),
                                     N_LS_FIELD)
                    if ls[-1] <= ls[0]:
                        continue
                    th, li = L.solve_theta_grid(np.full_like(ls, lam_p), ls,
                                                plane, seq, X, Y, Z,
                                                n_theta=N_THETA, branch=branch)
                    ok = np.isfinite(th)
                    if not ok.any():
                        continue
                    g = L.gvm_functions(np.full(ok.sum(), lam_p), ls[ok],
                                        li[ok], th[ok], plane, seq, X, Y, Z)
                    for k in range(ok.sum()):
                        thk = th[ok][k]
                        d = deff_bibo_pm_per_v(plane, seq[2], seq[1], seq[0],
                                               thk, lambda_pump_um=lam_p,
                                               enforce_range=False)
                        w.writerow([plane, seq, branch,
                                    round(lam_p * 1000, 2),
                                    round(float(ls[ok][k]), 5),
                                    round(float(li[ok][k]), 5),
                                    round(float(np.degrees(thk)), 4),
                                    round(float(g["theta_pm"][k]), 4),
                                    round(abs(d), 5)])
                        n += 1
                print(f"field {plane}/{seq}/{branch}: {n} nodes "
                      f"({time.time()-t0:.1f} s)", flush=True)


def stage_curves():
    """Trace every theta_PM = 0/45/90 level curve."""
    pumps = np.arange(PUMP_LO, PUMP_HI + 1e-9, PUMP_STEP_TRACE)
    cov = []
    with open("survey/landscape_curves.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CURVE_FIELDS)
        w.writeheader()
        for plane, seq in INTERACTIONS:
            for target in TARGETS:
                total = 0
                dmax = 0.0
                t0 = time.time()
                for branch in BRANCHES:
                    pts = L.trace_gvm_curve(plane, seq, target, X, Y, Z,
                                            pumps, n_ls=N_LS_TRACE,
                                            branch=branch, n_theta=N_THETA)
                    for p in pts:
                        row = enrich(p)
                        w.writerow(row)
                        total += 1
                        dmax = max(dmax, row["abs_deff_pm_V"])
                fh.flush()
                cov.append({"plane": plane, "interaction": seq,
                            "target_deg": target, "n_points": total,
                            "max_abs_deff_pm_V": round(dmax, 5)})
                print(f"curve {plane}/{seq} target={target:>2}: {total} points, "
                      f"max|deff|={dmax:.3f} pm/V ({time.time()-t0:.1f} s)",
                      flush=True)
    with open("survey/landscape_coverage.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["plane", "interaction",
                                           "target_deg", "n_points",
                                           "max_abs_deff_pm_V"])
        w.writeheader()
        w.writerows(cov)


def stage_degenerate():
    """Frequency-degenerate member of each GVM family."""
    rows = []
    for plane, seq in INTERACTIONS:
        for target in TARGETS:
            for branch in BRANCHES:
                p = L.degenerate_gvm_point(plane, seq, target, X, Y, Z,
                                           (PUMP_LO, PUMP_HI), n_scan=140,
                                           n_theta=N_THETA, branch=branch)
                if p is None:
                    continue
                row = enrich(p)
                if any(abs(row["pump_nm"] - r["pump_nm"]) < 0.5
                       and r["interaction"] == seq and r["plane"] == plane
                       for r in rows):
                    continue
                rows.append(row)
                print(f"degenerate {plane}/{seq} target={target:>2} "
                      f"({branch}): pump {row['pump_nm']:.2f} nm -> "
                      f"{row['lambda_s_um']*1000:.1f} nm, theta="
                      f"{row['angle_deg']:.3f}, |deff|="
                      f"{row['abs_deff_pm_V']:.3f}", flush=True)
    with open("survey/landscape_degenerate.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CURVE_FIELDS)
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    os.makedirs("survey", exist_ok=True)
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage in ("degenerate", "all"):
        stage_degenerate()
    if stage in ("curves", "all"):
        stage_curves()
    if stage in ("field", "all"):
        stage_field()
