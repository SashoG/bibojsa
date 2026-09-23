"""Degenerate EPM (45-deg) locus off-plane: at degeneracy g45 = d(Delta k)/d(omega_p),
so the EPM point is the turning point of the degenerate phase-matching curve
theta_pm(lambda_p).  For each phi, minimise/maximise theta_pm over lambda_p."""
import numpy as np, csv
from scipy.optimize import brentq, minimize_scalar
import offplane as o

def theta_pm_deg(lp, phi, assign, lo=0.5, hi=89.5, n=180):
    mp, ms, mi = assign
    ths = np.linspace(lo, hi, n)
    f = np.array([o.dk(lp, 2*lp, o.direction(t, phi), mp, ms, mi) for t in ths])
    roots = [brentq(lambda t: o.dk(lp, 2*lp, o.direction(t, phi), mp, ms, mi), ths[i], ths[i+1], xtol=1e-6)
             for i in range(n-1) if f[i]*f[i+1] < 0]
    return roots

rows = []
for assign in (('f','s','f'),):
    for phi in np.arange(0.0, 90.0+1e-9, 5.0):
        # the x-z oeo degenerate PM angle has a minimum near 870 nm; track the branch nearest 45 deg
        def th_of_lp(lp):
            r = theta_pm_deg(lp, phi, assign)
            if not r: return np.nan
            return min(r, key=lambda t: abs(t-45))
        res = minimize_scalar(th_of_lp, bounds=(0.75, 1.0), method='bounded', options={'xatol':1e-6})
        lp, th = res.x, res.fun
        s = o.direction(th, phi)
        g0,g45,g90,d31,d32 = o.gvm(lp, 2*lp, s, *assign)
        de = o.deff(lp,2*lp,2*lp,s,*assign)
        rmax = max(o.walkoff(lp,s,assign[0]), o.walkoff(2*lp,s,assign[1]), o.walkoff(2*lp,s,assign[2]))*1e3
        # the mirror octant (180-theta)
        s2 = o.direction(180-th, phi); de2 = o.deff(lp,2*lp,2*lp,s2,*assign)
        rows.append(dict(assign='fsf', target='45', phi_deg=phi, theta_deg=th, pump_nm=lp*1e3, pair_um=2*lp,
                         g45_fs_mm=g45, d31_fs_mm=d31, abs_deff=abs(de), abs_deff_mirror=abs(de2), rho_max_mrad=rmax,
                         fom=abs(de)/rmax, fom_mirror=abs(de2)/rmax))
        print(f"  fsf 45deg phi={phi:5.1f} theta={th:7.3f} lp={lp*1e3:8.2f} nm g45={g45:+.3f} |D31|={abs(d31):.1f} "
              f"|deff|={abs(de):.3f} / mirror {abs(de2):.3f}  rho_max={rmax:.1f}  fom={abs(de)/rmax:.4f} / {abs(de2)/rmax:.4f}")
with open('offplane_epm_locus.csv','w',newline='') as fh:
    w=csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
