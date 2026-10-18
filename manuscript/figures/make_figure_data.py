"""Derive the plot-ready CSVs for the manuscript figures from the survey.

Reads   ../../bibojsa/survey/*.csv   (relative to manuscript/figures/)
Writes  data/*.csv

  fig1  gvmcurve_xz_<cond>.csv     x-z type-II level curves (both halves of the
        gvmcurve_xy_<cond>.csv     oeo = eoo surface merged into one continuous
                                   curve per physical condition; x-y eoe = oee
                                   likewise), columns lam_e_um / lam_o_um
        gvmcurve_<seq>_<target>.csv  type-I families (ooe, eeo), as before
  fig2  deff_xz_<cond>.csv, deff_xy_<cond>.csv, deff_<seq>_<target>.csv
                                   |deff| and group-delay scale along each curve
  fig3  purity_xz_<cond>.csv, purity_xy_<cond>.csv, purity_<seq>_<target>.csv
                                   best converged purity along each curve
  fig4  jsa_<tag>.csv              |f| in angular-frequency detuning, cropped
        jsa_mesh.tex               mesh/cols, mesh/rows per panel
  fig5  filter_<tag>.csv           purity vs transmitted heralding flux

cond:  ke  = pump matched to the e wave (k_p' = k_e'), x-z: oeo/0 + eoo/90
       epm = extended phase matching,                 x-z: oeo/45 + eoo/45
       ko  = pump matched to the o wave (k_p' = k_o'), x-z: oeo/90 + eoo/0,
                                                      x-y: oee/90 (o = z-polarized)
Each merged curve is ordered along the curve by the e-wave wavelength, which
is monotone through the degenerate member (17 Oct 2026, v0.11).

Run from manuscript/figures/:   python3 make_figure_data.py
"""
import csv
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SURVEY = os.path.normpath(os.path.join(HERE, "..", "..", "bibojsa", "survey"))
OUT = os.path.join(HERE, "data")

JSA_TAGS = ["A", "H", "G", "D", "K", "E"]           # panels of Fig. 4
FILTER_TAGS = ["A", "B", "C", "D", "G", "H", "K", "E"]  # curves of Fig. 5
TYPE_I = ["ooe", "eeo"]
TARGETS = [0, 45, 90]
COND = {("xz", "oeo", 0): ("xz", "ke"), ("xz", "eoo", 90): ("xz", "ke"),
        ("xz", "oeo", 45): ("xz", "epm"), ("xz", "eoo", 45): ("xz", "epm"),
        ("xz", "oeo", 90): ("xz", "ko"), ("xz", "eoo", 0): ("xz", "ko"),
        ("xy", "oee", 90): ("xy", "ko"), ("xy", "oee", 45): ("xy", "epm"),
        ("xy", "eoe", 90): ("xy", "ko"), ("xy", "eoe", 45): ("xy", "epm"), ("xy", "eoe", 0): ("xy", "ke")}


def rd(name):
    with open(os.path.join(SURVEY, name)) as fh:
        return list(csv.DictReader(fh))


def wr(name, header, rows):
    with open(os.path.join(OUT, name), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"  data/{name}: {len(rows)} rows")


def e_and_o(r):
    """(lambda_e, lambda_o) of a row from the idler-signal-pump sequence."""
    seq = r["interaction"]
    ls, li = float(r["lambda_s_um"]), float(r["lambda_i_um"])
    return (ls, li) if seq[1] == "e" else (li, ls)


def split_branches(rows, keyfun, gap_nm=10.0):
    """Order a traced level curve along the curve by `keyfun` (monotone along
    it), breaking only on a genuine gap in pump wavelength; never drop points."""
    rows = sorted(rows, key=keyfun)
    out, cur = [], []
    for r in rows:
        if cur and abs(float(r["pump_nm"]) - float(cur[-1]["pump_nm"])) > gap_nm:
            out.append(cur)
            cur = []
        cur.append(r)
    if cur:
        out.append(cur)
    return out


def curves_and_deff():
    curves = [r for r in rd("landscape_curves.csv") if float(r["abs_deff_pm_V"]) >= 0.5]
    # --- type II, merged by physical condition
    fam = {}
    for r in curves:
        key = (r["plane"], r["interaction"], int(r["target_deg"]))
        if key in COND:
            fam.setdefault(COND[key], []).append(r)
    for (plane, cond), rows in sorted(fam.items()):
        branches = split_branches(rows, lambda r: e_and_o(r)[0])
        out, d = [], []
        for b in branches:
            for r in b:
                le, lo = e_and_o(r)
                out.append([round(float(r["pump_nm"]) / 1000.0, 6), round(le, 6), round(lo, 6),
                            r["lambda_i_in_transparency"], r["interaction"]])
                d.append([round(float(r["pump_nm"]) / 1000.0, 6), r["abs_deff_pm_V"], r["gvm_scale_fs_mm"],
                          round(max(abs(float(r["rho_pump_mrad"])), abs(float(r["rho_signal_mrad"])),
                                    abs(float(r["rho_idler_mrad"]))), 3)])
            out.append(["nan", "nan", "nan", "0", ""])
            d.append(["nan", "nan", "nan", "nan"])
        wr(f"gvmcurve_{plane}_{cond}.csv", ["lam_p_um", "lam_e_um", "lam_o_um", "inwin", "seq"], out)
        wr(f"deff_{plane}_{cond}.csv", ["lam_p_um", "abs_deff", "gvm_fs_mm", "rho_mrad"], d)
    # --- type I, as before (daughters interchangeable, one half is the whole surface)
    fam1 = {}
    for r in curves:
        if r["plane"] == "xz" and r["interaction"] in TYPE_I:
            fam1.setdefault((r["interaction"], int(r["target_deg"])), []).append(r)
    for (seq, tgt), rows in sorted(fam1.items()):
        branches = split_branches(rows, lambda r: float(r["lambda_s_um"]))
        out, d = [], []
        for b in branches:
            for r in b:
                out.append([round(float(r["pump_nm"]) / 1000.0, 6), r["lambda_s_um"], r["lambda_i_um"],
                            r["lambda_i_in_transparency"]])
                d.append([round(float(r["pump_nm"]) / 1000.0, 6), r["abs_deff_pm_V"], r["gvm_scale_fs_mm"]])
            out.append(["nan", "nan", "nan", "0"])
            d.append(["nan", "nan", "nan"])
        wr(f"gvmcurve_{seq}_{tgt}.csv", ["lam_p_um", "lam_s_um", "lam_i_um", "inwin"], out)
        wr(f"deff_{seq}_{tgt}.csv", ["lam_p_um", "abs_deff", "gvm_fs_mm"], d)


def purity():
    pur = rd("family_purity.csv")
    fam, fam1 = {}, {}
    for r in pur:
        key = (r["plane"], r["interaction"], int(r["target_deg"]))
        if key in COND:
            fam.setdefault(COND[key], []).append(r)
        elif r["plane"] == "xz" and r["interaction"] in TYPE_I:
            fam1.setdefault((r["interaction"], int(r["target_deg"])), []).append(r)
    for (plane, cond), rows in sorted(fam.items()):
        rows.sort(key=lambda r: e_and_o(r)[0])
        wr(f"purity_{plane}_{cond}.csv", ["lam_p_um", "purity", "purity_filt1x", "trans", "inwin"],
           [[round(float(r["pump_nm"]) / 1000.0, 6), r["best_purity"], r["purity_filt1x"],
             r["trans_filt1x"], r["lambda_i_in_transparency"]] for r in rows])
    for (seq, tgt), rows in sorted(fam1.items()):
        rows.sort(key=lambda r: float(r["lambda_s_um"]))
        wr(f"purity_{seq}_{tgt}.csv", ["lam_p_um", "purity", "purity_filt1x", "trans", "inwin"],
           [[round(float(r["pump_nm"]) / 1000.0, 6), r["best_purity"], r["purity_filt1x"],
             r["trans_filt1x"], r["lambda_i_in_transparency"]] for r in rows])


JSA_COLS = 150
# Axis window of each panel of fig5_jsa_maps.tex (symmetric, rad/ps); keep in
# step with xmin/xmax there.  The survey grids are cropped to these windows
# BEFORE subsampling.
JSA_WINDOW = {"A": 70, "H": 65, "G": 89, "D": 115, "K": 90, "E": 509}


def jsa_and_filter():
    cols = {}
    for tag in JSA_TAGS:
        with open(os.path.join(SURVEY, f"headline_jsa_{tag}.csv")) as fh:
            rows = list(csv.reader(fh))
        head, body = rows[0], rows[1:]
        win = JSA_WINDOW[tag]
        xs_all = sorted({float(r[0]) for r in body})
        ys_all = sorted({float(r[1]) for r in body})
        dx, dy = xs_all[1] - xs_all[0], ys_all[1] - ys_all[0]
        body = [r for r in body if abs(float(r[0])) <= win + dx and abs(float(r[1])) <= win + dy]
        xs = sorted({r[0] for r in body}, key=float)
        ys = sorted({r[1] for r in body}, key=float)
        keepx = {xs[i] for i in np.unique(np.linspace(0, len(xs) - 1, min(len(xs), JSA_COLS)).astype(int))}
        keepy = {ys[i] for i in np.unique(np.linspace(0, len(ys) - 1, min(len(ys), JSA_COLS)).astype(int))}
        out = [r for r in body if r[0] in keepx and r[1] in keepy]
        cols[tag] = (len(keepx), len(keepy))
        wr(f"jsa_{tag}.csv", head, out)
    with open(os.path.join(OUT, "jsa_mesh.tex"), "w") as fh:
        fh.write("% written by make_figure_data.py -- mesh/cols and\n")
        fh.write("% mesh/rows of data/jsa_<tag>.csv for fig5_jsa_maps.tex\n")
        for tag, (nx, ny) in cols.items():
            fh.write(f"\\expandafter\\def\\csname jsacols{tag}\\endcsname{{{nx}}}\n")
            fh.write(f"\\expandafter\\def\\csname jsarows{tag}\\endcsname{{{ny}}}\n")
    print("  data/jsa_mesh.tex:", cols)
    by = {}
    for r in rd("headline_filter.csv"):
        by.setdefault(r["tag"], []).append(r)
    for tag in FILTER_TAGS:
        rows = sorted(by[tag], key=lambda r: float(r["transmission"]))
        wr(f"filter_{tag}.csv", ["transmission", "purity", "ratio", "filter_fwhm_nm"],
           [[r["transmission"], r["purity"], r["ratio"], r["filter_fwhm_nm"]] for r in rows])


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for f in os.listdir(OUT):
        os.remove(os.path.join(OUT, f))
    print("fig1 / fig2 data"); curves_and_deff()
    print("fig3 data"); purity()
    print("fig4 / fig5 data"); jsa_and_filter()
    print("done")
