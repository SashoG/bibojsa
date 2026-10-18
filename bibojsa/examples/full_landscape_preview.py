"""Preview of the complete type-II landscape (both halves of each surface),
drawn from survey_full/landscape_curves.csv and the family_purity files.
Also writes merged per-condition CSVs for the manuscript figures:

  survey_full/figdata/gvmcurve_xz_typeII_<cond>.csv   cond = kp_eq_ke, epm, kp_eq_ko
  survey_full/figdata/gvmcurve_xy_typeII_<cond>.csv   cond = epm, kp_eq_ko

with a column `half` ("e_shorter" = the sequence traced before 17 Oct 2026,
"e_longer" = the partner sequence) and the usual columns of
figures/data/gvmcurve_*.csv (pump_nm, lambda_s_um, lambda_i_um, deff, scale).

The physical conditions map onto the label convention as follows (x-z):
    k_p' = k_e'   = oeo target 0  U  eoo target 90
    EPM           = oeo target 45 U  eoo target 45
    k_p' = k_o'   = oeo target 90 U  eoo target 0
and (x-y, where the o wave is z-polarized):
    k_p' = k_o'   = oee target 90   (eoe has no crossing)
    EPM           = oee target 45
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "survey_full"
os.makedirs(f"{OUT}/figdata", exist_ok=True)
rows = list(csv.DictReader(open(f"{OUT}/landscape_curves.csv")))

COND = {
    ("xz", "oeo", "0"): ("xz", "kp_eq_ke", "e_shorter"), ("xz", "eoo", "90"): ("xz", "kp_eq_ke", "e_longer"),
    ("xz", "oeo", "45"): ("xz", "epm", "e_shorter"), ("xz", "eoo", "45"): ("xz", "epm", "e_longer"),
    ("xz", "oeo", "90"): ("xz", "kp_eq_ko", "e_shorter"), ("xz", "eoo", "0"): ("xz", "kp_eq_ko", "e_longer"),
    ("xy", "oee", "90"): ("xy", "kp_eq_ko", "o_longer"), ("xy", "oee", "45"): ("xy", "epm", "o_longer"),
}
merged = {}
for r in rows:
    key = (r["plane"], r["interaction"], r["target_deg"])
    if key not in COND:
        continue
    plane, cond, half = COND[key]
    # e-wave / o-wave wavelengths regardless of the signal/idler labels
    seq = r["interaction"]
    ls, li = float(r["lambda_s_um"]), float(r["lambda_i_um"])
    pol_s, pol_i = seq[1], seq[0]
    lam_e = ls if pol_s == "e" else li
    lam_o = li if pol_s == "e" else ls
    if plane == "xy":   # x-y: 'o' is the z-polarized wave, 'e' in-plane
        lam_e = ls if pol_s == "e" else li
        lam_o = li if pol_s == "e" else ls
    merged.setdefault((plane, cond), []).append({
        "half": half, "pump_nm": float(r["pump_nm"]), "lambda_s_um": ls, "lambda_i_um": li,
        "lambda_e_um": lam_e, "lambda_o_um": lam_o, "angle_deg": float(r["angle_deg"]),
        "abs_deff_pm_V": float(r["abs_deff_pm_V"]), "gvm_scale_fs_mm": float(r["gvm_scale_fs_mm"]),
        "rho_max_mrad": max(abs(float(r["rho_pump_mrad"])), abs(float(r["rho_signal_mrad"])), abs(float(r["rho_idler_mrad"]))),
        "lambda_i_in_transparency": int(r["lambda_i_in_transparency"])})
for (plane, cond), rs in merged.items():
    rs.sort(key=lambda x: (x["pump_nm"], x["lambda_s_um"]))
    with open(f"{OUT}/figdata/gvmcurve_{plane}_typeII_{cond}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rs[0])); w.writeheader(); w.writerows(rs)
    p = np.array([x["pump_nm"] for x in rs])
    print(f"{plane} {cond}: {len(rs)} points, pump {p.min():.0f}-{p.max():.0f} nm, "
          f"|deff| {min(x['abs_deff_pm_V'] for x in rs):.2f}-{max(x['abs_deff_pm_V'] for x in rs):.2f}, "
          f"scale {min(x['gvm_scale_fs_mm'] for x in rs):.0f}-{max(x['gvm_scale_fs_mm'] for x in rs):.0f} fs/mm, "
          f"rho {min(x['rho_max_mrad'] for x in rs):.0f}-{max(x['rho_max_mrad'] for x in rs):.0f} mrad")

# ------------------------------------------------------------- preview figure
fp_old = list(csv.DictReader(open("survey/family_purity.csv")))
fp_new = list(csv.DictReader(open(f"{OUT}/family_purity_new.csv"))) if os.path.exists(f"{OUT}/family_purity_new.csv") else []
fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
col = {"kp_eq_ke": "#1f77b4", "epm": "#ff7f0e", "kp_eq_ko": "#2ca02c"}
name = {"kp_eq_ke": "k_p' = k_e'", "epm": "EPM", "kp_eq_ko": "k_p' = k_o'"}
for (plane, cond), rs in merged.items():
    a = ax[0, 0] if plane == "xz" else ax[0, 1]
    for half, mk in (("e_shorter", "."), ("e_longer", "x"), ("o_longer", "x")):
        q = [x for x in rs if x["half"] == half]
        if not q:
            continue
        p = [x["pump_nm"] for x in q]
        a.plot(p, [x["lambda_e_um"] for x in q], mk, ms=3, color=col[cond],
               label=f"{name[cond]} ({'traced before' if half == 'e_shorter' else 'new'}) - e wave")
        a.plot(p, [x["lambda_o_um"] for x in q], mk, ms=3, color=col[cond], alpha=0.45)
    b = ax[1, 0] if plane == "xz" else ax[1, 1]
    for half, mk in (("e_shorter", "."), ("e_longer", "x"), ("o_longer", "x")):
        q = [x for x in rs if x["half"] == half]
        if q:
            b.plot([x["pump_nm"] for x in q], [x["abs_deff_pm_V"] / x["rho_max_mrad"] for x in q], mk, ms=3, color=col[cond])
for a, t in ((ax[0, 0], "x-z plane, oeo = eoo (dots: in paper; crosses: untraced half)"), (ax[0, 1], "x-y plane, eoe = oee (all new)")):
    a.set_title(t, fontsize=10); a.set_xlabel("pump (nm)"); a.set_ylabel("daughter wavelengths (um): e wave dark, o wave light")
    a.axhline(2.6, ls="-.", c="k", lw=0.6); a.legend(fontsize=6)
for b in (ax[1, 0], ax[1, 1]):
    b.set_xlabel("pump (nm)"); b.set_ylabel("d_eff / rho_max  (pm/V/mrad)")
    b.axhline(0.102, ls=":", c="k", lw=0.6); b.text(b.get_xlim()[0], 0.104, "point E (triple GVM)", fontsize=7)
    b.axhline(0.026, ls=":", c="gray", lw=0.6); b.text(b.get_xlim()[0], 0.0275, "point A", fontsize=7, color="gray")
# purity overlay on the lower panels (right axis)
for plane, b in (("xz", ax[1, 0]), ("xy", ax[1, 1])):
    c = b.twinx()
    for src, mk, lab in ((fp_old, ".", "best unfiltered P (in paper)"), (fp_new, "x", "best unfiltered P (new half)")):
        q = [r for r in src if r["plane"] == plane and r["interaction"] in ("oeo", "eoo", "oee")]
        if q:
            c.plot([float(r["pump_nm"]) for r in q], [float(r["best_purity"]) for r in q], mk, ms=4, color="purple", alpha=0.6, label=lab)
    c.set_ylabel("best unfiltered purity", color="purple"); c.set_ylim(0.6, 1.0); c.legend(fontsize=6, loc="lower right")
fig.suptitle("BIBO type-II GVM landscape on both halves of the phase-matching surface (survey_full, 17 Oct 2026)")
fig.tight_layout()
fig.savefig(f"{OUT}/fig_full_landscape_preview.png", dpi=130)
print("wrote", f"{OUT}/fig_full_landscape_preview.png")
