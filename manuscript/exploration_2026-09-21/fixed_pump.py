"""0-deg (asymmetric GVM, k_p'=k_s') family of the (f,s,f) interaction at fixed pump
wavelengths, traced in phi: does off-plane propagation extend the 678-722 nm
in-plane pump range of the manuscript's 0-deg family?"""
import numpy as np, csv, sys
from scipy.optimize import brentq
import offplane as o

LI_MAX = 3.05


def signal_roots(lp, s, assign, n=80):
    mp, ms, mi = assign
    ls_min = 1.0 / (1.0 / lp - 1.0 / LI_MAX) * 1.0005
    lss = np.linspace(ls_min, 2 * lp, n)
    f = np.array([o.dk(lp, l, s, mp, ms, mi) for l in lss])
    return [brentq(lambda l: o.dk(lp, l, s, mp, ms, mi), lss[i], lss[i + 1], xtol=1e-9)
            for i in range(n - 1) if np.isfinite(f[i]) and np.isfinite(f[i + 1]) and f[i] * f[i + 1] < 0]


def trace(lp, phis, assign=('f', 's', 'f'), th_lo=15.0, th_hi=75.0, step=0.25):
    out = []
    thetas = np.arange(th_lo, th_hi + 1e-9, step)
    for phi in phis:
        prev = None
        for th in thetas:
            s = o.direction(th, phi)
            cur = [(ls, o.gvm(lp, ls, s, *assign)[0]) for ls in signal_roots(lp, s, assign)]
            if prev is not None and len(prev) == len(cur) and cur:
                for b, ((ls0, g0p), (ls1, g0c)) in enumerate(zip(prev, cur)):
                    if abs(ls1 - ls0) > 0.1 or not (g0p * g0c < 0):
                        continue

                    def G(t):
                        ss = o.direction(t, phi); rr = signal_roots(lp, ss, assign)
                        return np.nan if len(rr) != len(cur) else o.gvm(lp, rr[b], ss, *assign)[0]
                    try:
                        tz = brentq(G, th - step, th, xtol=1e-4)
                    except Exception:
                        continue
                    ss = o.direction(tz, phi); ls = signal_roots(lp, ss, assign)[b]; li = 1 / (1 / lp - 1 / ls)
                    g0, g45, g90, d31, d32 = o.gvm(lp, ls, ss, *assign)
                    de = o.deff(lp, ls, li, ss, *assign)
                    rmax = max(o.walkoff(lp, ss, assign[0]), o.walkoff(ls, ss, assign[1]), o.walkoff(li, ss, assign[2])) * 1e3
                    out.append(dict(pump_nm=lp * 1e3, phi_deg=phi, theta_deg=tz, lambda_s_um=ls, lambda_i_um=li,
                                    d31_fs_mm=d31, abs_deff=abs(de), rho_max_mrad=rmax, fom=abs(de) / rmax))
                    print(f"  lp={lp*1e3:.0f} phi={phi:5.1f} theta={tz:7.3f} ls={ls:.4f} li={li:.4f} |D31|={abs(d31):.1f} fs/mm "
                          f"|deff|={abs(de):.3f} rho={rmax:.1f} fom={abs(de)/rmax:.4f}", flush=True)
            prev = cur
    return out


rows = []
for lpnm in [float(x) for x in sys.argv[1:]] or (722.0, 740.0, 760.0, 780.0, 800.0):
    rows += trace(lpnm * 1e-3, np.arange(0.0, 90.0 + 1e-9, 10.0))
if rows:
    with open('offplane_zero_deg_fixed_pump.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("done")
