"""Append 1-nm-pump rows for G and H to survey/length_extension.csv.

20 Oct 2026 (physics check, issue m19): the manuscript quotes H's 30-mm
optimum, 0.968 with a 1-nm pump, which lay below the 2-nm edge of the
length_extension.csv grid.  This script computes the 1-nm rows for G and H
with the same functions, window set (6/12/24 sigma) and pixel as
window_tail.py's length_rows() and inserts them before each (tag, L) block,
leaving every existing row untouched.  (A full regeneration with
`python examples/window_tail.py length` -- whose grid now includes 1 nm --
recomputes the A, B, C rows at the 0.25-rad/ps pixel as well.)

Run from bibojsa/:  python examples/add_1nm_rows.py
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from window_tail import FEATURED, OUT, PIXEL_RAD_PS, axes_for, purity_on

from bibojsa.jsa import fwhm_omega_from_fwhm_wavelength, gaussian_pump_sigma
from bibojsa.sellmeier import lambda_um_to_omega

rows = []
for pt in [f for f in FEATURED if f[0] in ("G", "H")]:
    tag, plane, seq, _, pump_nm, lam_s, lam_i, theta_deg, *_ = pt
    lam_p = pump_nm / 1000.0
    wp0, ws0, wi0 = (lambda_um_to_omega(v) for v in (lam_p, lam_s, lam_i))
    for L_mm in (10.0, 15.0, 20.0, 30.0):
        bw = 1.0
        fw = fwhm_omega_from_fwhm_wavelength(bw * 1e-3, lam_p)
        sig = gaussian_pump_sigma(fw)
        Ps, Ws = [], []
        for w in (6.0, 12.0, 24.0):
            ws, wi, clipped = axes_for(ws0, wi0, w * sig, PIXEL_RAD_PS)
            if clipped:
                continue
            Ps.append(purity_on(pt, L_mm, bw, ws, wi).purity)
            Ws.append(w)
        A = np.vstack([np.ones(len(Ws)), 1.0 / np.array(Ws)]).T
        (p_inf, c), *_ = np.linalg.lstsq(A, np.array(Ps), rcond=None)
        rows.append({"tag": tag, "L_mm": L_mm, "pump_fwhm_nm": bw,
                     "P_6sig": round(Ps[0], 5), "P_inf": round(p_inf, 4),
                     "C_per_sigma": round(c, 4)})
        print(rows[-1], flush=True)

path = os.path.join(OUT, "length_extension.csv")
with open(path) as fh:
    old = list(csv.DictReader(fh))
out, used = [], []
for r in old:
    for n in rows:
        if (n["tag"] == r["tag"] and float(n["L_mm"]) == float(r["L_mm"])
                and n not in used):
            out.append({k: str(v) for k, v in n.items()})
            used.append(n)
    out.append(r)
with open(path, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(old[0].keys()))
    w.writeheader()
    w.writerows(out)
print(len(old), "->", len(out), "rows")
