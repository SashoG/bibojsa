"""Derive the plot-ready CSVs for the manuscript figures from the survey.

Reads   ../../bibojsa/survey/*.csv   (relative to manuscript/figures/)
Writes  data/*.csv

  fig1  gvmcurve_<seq>_<target>.csv   GVM level curves in the (lambda_p,
                                      lambda_s) plane (the 2.60/3.05 um idler
                                      limits are drawn as constants in the .tex)
  fig2  tpmdeg_<seq>.csv              theta_PM at frequency degeneracy vs pump
  fig3  deff_<seq>_<target>.csv       |deff| along each family vs pump
  fig4  purity_<seq>_<target>.csv     best converged purity along each family
  fig5  jsa_<tag>.csv                 |f| in angular-frequency detuning,
                                      cropped to the panel window (A-E, R)
        jsa_mesh.tex                  mesh/cols, mesh/rows per panel
  fig6  filter_<tag>.csv              purity vs transmitted heralding flux
                                      (A-E, R)

Only files that a figure reads are written.  Point F is in Table 2 but in
neither Fig. 5 nor Fig. 6 (qualitatively identical to E), so no jsa_F/filter_F.

Run from manuscript/figures/:   python3 make_figure_data.py
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SURVEY = os.path.normpath(os.path.join(HERE, "..", "..", "bibojsa", "survey"))
OUT = os.path.join(HERE, "data")

FIGURE_TAGS = ["A", "B", "C", "D", "E", "R"]   # panels of Fig. 5 / curves of Fig. 6
XZ_SEQS = ["oeo", "ooe", "eeo"]
TARGETS = [0, 45, 90]


def rd(name):
    with open(os.path.join(SURVEY, name)) as fh:
        return list(csv.DictReader(fh))


def wr(name, header, rows):
    with open(os.path.join(OUT, name), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"  data/{name}: {len(rows)} rows")


def split_branches(rows, gap_nm=10.0):
    """Order a traced level curve along the curve.

    landscape_curves.csv is sorted by pump wavelength, and every x-z family
    folds back in pump, so the two arms of a family are interleaved there
    (two rows per pump).  The signal wavelength, by contrast, increases
    monotonically along each traced family -- ordered by lambda_s, consecutive
    points are always one 2-nm pump step apart -- so sorting on lambda_s
    recovers the curve.  Break only on a genuine gap in pump wavelength, and
    never drop points.
    """
    rows = sorted(rows, key=lambda r: float(r["lambda_s_um"]))
    out, cur = [], []
    for r in rows:
        if cur and abs(float(r["pump_nm"])
                       - float(cur[-1]["pump_nm"])) > gap_nm:
            out.append(cur)
            cur = []
        cur.append(r)
    if cur:
        out.append(cur)
    return out


def fig1_and_3_and_4():
    curves = [r for r in rd("landscape_curves.csv")
              if r["plane"] == "xz" and float(r["abs_deff_pm_V"]) >= 0.5]
    fam = {}
    for r in curves:
        fam.setdefault((r["interaction"], int(r["target_deg"])), []).append(r)

    for (seq, tgt), rows in sorted(fam.items()):
        branches = split_branches(rows)
        out = []
        for b in branches:
            for r in b:
                out.append([round(float(r["pump_nm"]) / 1000.0, 6),
                            r["lambda_s_um"], r["lambda_i_um"],
                            r["lambda_i_in_transparency"]])
            out.append(["nan", "nan", "nan", "0"])   # break the line
        wr(f"gvmcurve_{seq}_{tgt}.csv",
           ["lam_p_um", "lam_s_um", "lam_i_um", "inwin"], out)

        d = []
        for b in branches:
            for r in b:
                d.append([round(float(r["pump_nm"]) / 1000.0, 6),
                          r["abs_deff_pm_V"], r["gvm_scale_fs_mm"]])
            d.append(["nan", "nan", "nan"])
        wr(f"deff_{seq}_{tgt}.csv",
           ["lam_p_um", "abs_deff", "gvm_fs_mm"], d)

    # fig4: purity along each family
    pur = rd("family_purity.csv")
    fam4 = {}
    for r in pur:
        fam4.setdefault((r["interaction"], int(r["target_deg"])), []).append(r)
    for (seq, tgt), rows in sorted(fam4.items()):
        # order along the family (see split_branches), not by pump
        rows.sort(key=lambda r: float(r["lambda_s_um"]))
        wr(f"purity_{seq}_{tgt}.csv",
           ["lam_p_um", "purity", "purity_filt1x", "trans", "inwin"],
           [[round(float(r["pump_nm"]) / 1000.0, 6), r["best_purity"],
             r["purity_filt1x"], r["trans_filt1x"],
             r["lambda_i_in_transparency"]] for r in rows])


def fig2():
    """theta_PM along the frequency-degenerate line, per interaction."""
    sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..",
                                                     "bibojsa")))
    from bibojsa import landscape as L
    from bibojsa.materials.bibo import BIBO
    X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]
    for seq in XZ_SEQS:
        rows = []
        for lp in np.arange(0.60, 1.301, 0.005):
            lam_d = 2 * lp
            if lam_d > L.LAM_MAX:
                break
            th, _ = L.solve_theta_grid(float(lp), float(lam_d) * 0.9999999,
                                       "xz", seq, X, Y, Z)
            if not np.isfinite(th):
                continue
            g = L.gvm_functions(lp, lam_d, lam_d, th, "xz", seq, X, Y, Z)
            # theta_PM is an orientation, defined mod 180 deg; wrap into
            # (-90, 90] so that the anticorrelated orientation is a single
            # branch at -45 deg rather than jumping between -45 and +135.
            t = (float(g["theta_pm"]) + 90.0) % 180.0 - 90.0
            if rows and t - rows[-1][1] < -120.0:
                t += 180.0          # unwrap: keep the branch continuous
            rows.append([round(float(lp), 6), round(t, 4),
                         round(float(np.degrees(th)), 4)])
        wr(f"tpmdeg_{seq}.csv", ["lam_p_um", "theta_pm_deg", "theta_deg"], rows)


JSA_COLS = 150

# Axis window of each panel of fig5_jsa_maps.tex (symmetric, rad/ps).  The
# survey grids are cropped to these windows BEFORE subsampling, so the
# side-lobes are drawn at the survey's own resolution wherever that fits in
# JSA_COLS cells.  Keep in step with the xmin/xmax of fig5_jsa_maps.tex.
JSA_WINDOW = {"A": 70, "B": 80, "C": 60, "D": 115, "E": 509, "R": 92}


def fig5_and_6():
    cols = {}
    for tag in FIGURE_TAGS:
        src = os.path.join(SURVEY, f"headline_jsa_{tag}.csv")
        if not os.path.exists(src):
            continue
        with open(src) as fh:
            rows = list(csv.reader(fh))
        head, body = rows[0], rows[1:]
        win = JSA_WINDOW.get(tag)
        if win is not None:
            # keep one grid point beyond the window on each side so the
            # cells at the axis edge are drawn in full
            xs_all = sorted({float(r[0]) for r in body})
            ys_all = sorted({float(r[1]) for r in body})
            dx = xs_all[1] - xs_all[0]
            dy = ys_all[1] - ys_all[0]
            body = [r for r in body
                    if abs(float(r[0])) <= win + dx and abs(float(r[1])) <= win + dy]
        xs = sorted({r[0] for r in body}, key=float)
        ys = sorted({r[1] for r in body}, key=float)
        # Subsample onto at most a JSA_COLS x JSA_COLS mesh: enough to
        # resolve the ridge, small enough that fig5 compiles with lualatex in a few minutes.
        keepx = {xs[i] for i in np.unique(
            np.linspace(0, len(xs) - 1, min(len(xs), JSA_COLS)).astype(int))}
        keepy = {ys[i] for i in np.unique(
            np.linspace(0, len(ys) - 1, min(len(ys), JSA_COLS)).astype(int))}
        out = [r for r in body if r[0] in keepx and r[1] in keepy]
        cols[tag] = (len(keepx), len(keepy))
        wr(f"jsa_{tag}.csv", head, out)
    # mesh dimensions for pgfplots' matrix plot (one macro per panel)
    with open(os.path.join(OUT, "jsa_mesh.tex"), "w") as fh:
        fh.write("% written by make_figure_data.py -- mesh/cols and\n")
        fh.write("% mesh/rows of data/jsa_<tag>.csv for fig5_jsa_maps.tex\n")
        for tag, (nx, ny) in cols.items():
            fh.write(f"\\expandafter\\def\\csname jsacols{tag}\\endcsname{{{nx}}}\n")
            fh.write(f"\\expandafter\\def\\csname jsarows{tag}\\endcsname{{{ny}}}\n")
    print("  data/jsa_mesh.tex:", cols)
    filt = rd("headline_filter.csv")
    by = {}
    for r in filt:
        by.setdefault(r["tag"], []).append(r)
    for tag, rows in by.items():
        if tag not in FIGURE_TAGS:
            continue
        # keep the transmission = 1 row (empty ratio): it is the unfiltered
        # purity that the caption of Fig. 6 refers to
        rows.sort(key=lambda r: float(r["transmission"]))
        wr(f"filter_{tag}.csv", ["transmission", "purity", "ratio",
                                 "filter_fwhm_nm"],
           [[r["transmission"], r["purity"], r["ratio"], r["filter_fwhm_nm"]]
            for r in rows])


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    print("fig1 / fig3 / fig4 data")
    fig1_and_3_and_4()
    print("fig2 data")
    fig2()
    print("fig5 / fig6 data")
    fig5_and_6()
    print("done")
