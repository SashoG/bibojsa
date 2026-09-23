"""Off-principal-plane exploration for the BIBO GVM manuscript (21 Sep 2026).

Independent numpy/scipy implementation (no bibojsa import) of:
  * biaxial eigenindices and D-eigenvectors for an arbitrary direction
    (index-ellipsoid section, 2x2 eigenproblem in the plane normal to k);
  * group index k'(omega) along a fixed direction (Richardson derivative);
  * walk-off angle rho = angle(D, E) with E_i = D_i / n_i^2;
  * effective nonlinearity from the full point-group-2 d tensor (2-fold
    axis = dielectric x, Kleinman-averaged LPR2010 coefficients) with the
    D-eigenvectors as polarization unit vectors (walk-off correction to
    d_eff neglected, as in LPR2010 and the manuscript).

Interactions are labelled by fast (f, smaller index) / slow (s) eigenmode,
in the LPR2010 order (idler, signal, pump).  In the x-z plane between the
optic axes the o wave (n_y) is the fast one, so the manuscript's oeo is
(f, s, f) here and eoe is (s, f, s).

Outputs (CSV + PNG) go to the same folder.
"""
import numpy as np
from scipy.optimize import brentq
import csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
C_UM_PS = 299.792458  # um/ps
SELL = {'x': (3.0759, 0.03169, 0.03323, 0.01402),
        'y': (3.1698, 0.03666, 0.03599, 0.01819),
        'z': (3.6546, 0.05116, 0.03713, 0.02299)}
LMIN, LMAX = 0.34, 3.05
# Kleinman-averaged d_il (pm/V), 2-fold axis along x: nonzero
# d11 d12 d13 d14 d25 d26 d35 d36
D_IL = np.zeros((3, 6))
D_IL[0, 0] = 2.53; D_IL[0, 1] = 3.2; D_IL[0, 2] = -1.76; D_IL[0, 3] = 1.66
D_IL[1, 4] = 1.66; D_IL[1, 5] = 3.2
D_IL[2, 4] = -1.76; D_IL[2, 5] = 1.66
VOIGT = {0: [(0, 0)], 1: [(1, 1)], 2: [(2, 2)], 3: [(1, 2), (2, 1)],
         4: [(0, 2), (2, 0)], 5: [(0, 1), (1, 0)]}
DTEN = np.zeros((3, 3, 3))
for i in range(3):
    for l, jks in VOIGT.items():
        for (j, k) in jks:
            DTEN[i, j, k] = D_IL[i, l]


def n_axes(lam):
    out = []
    for ax in 'xyz':
        A, B, C, D = SELL[ax]
        out.append(np.sqrt(A + B / (lam ** 2 - C) - D * lam ** 2))
    return np.array(out)


def direction(theta_deg, phi_deg):
    t, p = np.radians(theta_deg), np.radians(phi_deg)
    return np.array([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)])


def eigen(lam, s):
    """Return ((n_fast, D_fast), (n_slow, D_slow)) for direction s."""
    a = 1.0 / n_axes(lam) ** 2                     # eigenvalues of eps^-1
    # orthonormal basis of the plane normal to s
    ref = np.array([0.0, 0.0, 1.0]) if abs(s[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = np.cross(s, ref); u /= np.linalg.norm(u)
    v = np.cross(s, u)
    M = np.array([[np.sum(u * a * u), np.sum(u * a * v)],
                  [np.sum(v * a * u), np.sum(v * a * v)]])
    w, V = np.linalg.eigh(M)                       # ascending eigenvalues 1/n^2
    # larger 1/n^2 -> smaller n -> fast
    Df = V[0, 1] * u + V[1, 1] * v
    Ds = V[0, 0] * u + V[1, 0] * v
    return (1.0 / np.sqrt(w[1]), Df / np.linalg.norm(Df)), (1.0 / np.sqrt(w[0]), Ds / np.linalg.norm(Ds))


def n_mode(lam, s, mode):
    f, sl = eigen(lam, s)
    return f[0] if mode == 'f' else sl[0]


def kprime(lam, s, mode, h=2e-4):
    """dk/domega in ps/um, k = n omega / c, along fixed direction s.
    Richardson-extrapolated central difference in lambda."""
    def k(l):
        return n_mode(l, s, mode) * (2 * np.pi / l)   # rad/um
    def dk_dl(l, hh):
        return (k(l + hh) - k(l - hh)) / (2 * hh)
    d1 = dk_dl(lam, h); d2 = dk_dl(lam, h / 2)
    dkdl = (4 * d2 - d1) / 3
    # omega = 2 pi c / lam  ->  d omega / d lam = -2 pi c / lam^2
    return dkdl / (-2 * np.pi * C_UM_PS / lam ** 2)


def walkoff(lam, s, mode):
    f, sl = eigen(lam, s)
    D = f[1] if mode == 'f' else sl[1]
    a = 1.0 / n_axes(lam) ** 2
    E = a * D
    E /= np.linalg.norm(E)
    return np.arccos(np.clip(abs(np.dot(D, E)), -1, 1))   # rad


def deff(lp, ls, li, s, mp, ms, mi):
    """d_eff = sum e_p,i d_ijk e_s,j e_i,k with D-eigenvectors as polarization
    unit vectors (pm/V)."""
    ep = eigen(lp, s)[0 if mp == 'f' else 1][1]
    es = eigen(ls, s)[0 if ms == 'f' else 1][1]
    ei = eigen(li, s)[0 if mi == 'f' else 1][1]
    return float(np.einsum('i,ijk,j,k->', ep, DTEN, es, ei))


def dk(lp, ls, s, mp, ms, mi):
    """Collinear phase mismatch k_p - k_s - k_i in rad/um."""
    li = 1.0 / (1.0 / lp - 1.0 / ls)
    return (n_mode(lp, s, mp) / lp - n_mode(ls, s, ms) / ls - n_mode(li, s, mi) / li) * 2 * np.pi


def degenerate_pumps(s, mp, ms, mi, lo=LMIN, hi=1.52, n=120):
    """All lambda_p in [lo,hi] with degenerate phase matching along s."""
    lps = np.linspace(lo, hi, n)
    f = np.array([dk(lp, 2 * lp, s, mp, ms, mi) for lp in lps])
    roots = []
    for i in range(n - 1):
        if np.isfinite(f[i]) and np.isfinite(f[i + 1]) and f[i] * f[i + 1] < 0:
            roots.append(brentq(lambda lp: dk(lp, 2 * lp, s, mp, ms, mi), lps[i], lps[i + 1], xtol=1e-10))
    return roots


def gvm(lp, ls, s, mp, ms, mi):
    li = 1.0 / (1.0 / lp - 1.0 / ls)
    kp, ks, ki = kprime(lp, s, mp), kprime(ls, s, ms), kprime(li, s, mi)
    d32, d31 = (kp - ks) * 1e3, (kp - ki) * 1e3    # fs/mm
    g0, g45, g90 = d32, (d32 + d31) / 2, d31
    return g0, g45, g90, d31, d32


# ------------------------------------------------------------------ checks
def validate():
    print("== validation against the manuscript's x-z (phi=0) results ==")
    s = direction(45.968, 0.0)
    f, sl = eigen(0.6822, s)
    print(f"  pump 682.2 nm along theta=45.968: n_fast={f[0]:.5f} (n_y={n_axes(0.6822)[1]:.5f}), n_slow={sl[0]:.5f}; "
          f"D_fast={np.round(f[1],3)}, D_slow={np.round(sl[1],3)}")
    lp = degenerate_pumps(s, 'f', 's', 'f')
    print("  degenerate PM pumps for (f,s,f) at theta=45.968, phi=0:", np.round(lp, 5))
    for l in lp:
        g0, g45, g90, d31, d32 = gvm(l, 2 * l, s, 'f', 's', 'f')
        print(f"    lp={l*1e3:.2f} nm: D32={d32:.3f} D31={d31:.3f} fs/mm, deff={deff(l,2*l,2*l,s,'f','s','f'):.4f} pm/V, "
              f"rho_slow={walkoff(2*l,s,'s')*1e3:.2f} mrad rho_fast={walkoff(2*l,s,'f')*1e3:.2f} mrad")
    # d_eff planar formulas
    th = np.radians(45.968)
    print(f"  LPR2010 x-z oeo: d12 cos(theta) = {3.2*np.cos(th):.4f}")
    s2 = direction(11.036, 0.0)
    print(f"  x-z ooe at theta=11.036 (f,f,s) deff={deff(0.8187,1.6374,1.6374,s2,'f','f','s'):.4f}  vs -d14 sin2theta = {-1.66*np.sin(2*np.radians(11.036)):.4f}")
    # y-z plane eeo: theta measured from z, phi=90
    for thd in (11.1, 168.9):
        s3 = direction(thd, 90.0)
        t = np.radians(thd)
        val = deff(1.064, 2.128, 2.128, s3, 's', 's', 'f')
        print(f"  y-z eeo (s,s,f) theta={thd}: |deff|={abs(val):.4f} vs LPR2010 |-/+d12cos2 -/+d13 sin2 + d14 sin2t| = "
              f"{abs(-3.2*np.cos(t)**2 - (-1.76)*np.sin(t)**2 + 1.66*np.sin(2*t)):.4f} / {abs(3.2*np.cos(t)**2 + (-1.76)*np.sin(t)**2 + 1.66*np.sin(2*t)):.4f}")
    # x-y plane: theta=90, phi
    s4 = direction(90.0, 30.0)
    print(f"  x-y (f,f,s) at phi=30: |deff|={abs(deff(0.4,0.8,0.8,s4,'f','f','s')):.4f} vs |d13 sin phi|={abs(-1.76*np.sin(np.radians(30))):.4f};"
          f"  (s,f,s): {abs(deff(0.4,0.8,0.8,s4,'s','f','s')):.4f} vs |d14 sin2phi|={abs(1.66*np.sin(np.radians(60))):.4f}")
    # cross-check with bibojsa eigenindices_general if available
    try:
        sys.path.insert(0, os.path.join(HERE, '..', 'verify', 'bibojsa_copy'))
        from bibojsa.biaxial import eigenindices_general
        from bibojsa.materials.bibo import BIBO
        for (thd, phd) in ((45.968, 0.0), (50.0, 30.0), (120.0, 70.0)):
            r = eigenindices_general(1.364, np.radians(thd), np.radians(phd), BIBO['x'], BIBO['y'], BIBO['z'])
            f, sl = eigen(1.364, direction(thd, phd))
            print(f"  bibojsa eigenindices_general({thd},{phd}) -> {r} ; mine (fast, slow) = ({f[0]:.6f}, {sl[0]:.6f})")
    except Exception as e:
        print("  (bibojsa cross-check skipped:", e, ")")


# ------------------------------------------------ degenerate loci over (theta, phi)
def trace_degenerate_loci(assign, phis, theta_step=0.25):
    """For each phi, scan theta, find degenerate PM pumps, evaluate g0/g45/g90
    and bracket their zeros in theta (per PM branch)."""
    mp, ms, mi = assign
    rows = []
    thetas = np.arange(0.0, 180.0 + 1e-9, theta_step)
    for phi in phis:
        prev = None   # list of (lp, g0, g45, g90) per branch at previous theta
        for th in thetas:
            s = direction(th, phi)
            lps = degenerate_pumps(s, mp, ms, mi)
            cur = []
            for lp in lps:
                g0, g45, g90, d31, d32 = gvm(lp, 2 * lp, s, mp, ms, mi)
                cur.append((lp, g0, g45, g90))
            if prev is not None and len(prev) == len(cur) and len(cur) > 0:
                # pair branches by pump wavelength ordering
                for b, ((lp0, *g_prev), (lp1, *g_cur)) in enumerate(zip(prev, cur)):
                    if abs(lp1 - lp0) > 0.05:
                        continue
                    for tgt, name in ((0, '0'), (1, '45'), (2, '90')):
                        if g_prev[tgt] * g_cur[tgt] < 0:
                            # refine in theta by bisection on the composite function
                            def G(t):
                                ss = direction(t, phi)
                                r = degenerate_pumps(ss, mp, ms, mi)
                                if len(r) != len(cur):
                                    return np.nan
                                return gvm(r[b], 2 * r[b], ss, mp, ms, mi)[tgt]
                            try:
                                tz = brentq(G, th - theta_step, th, xtol=1e-4)
                            except Exception:
                                continue
                            ss = direction(tz, phi)
                            lp = degenerate_pumps(ss, mp, ms, mi)[b]
                            g0, g45, g90, d31, d32 = gvm(lp, 2 * lp, ss, mp, ms, mi)
                            de = deff(lp, 2 * lp, 2 * lp, ss, mp, ms, mi)
                            rp, rs, ri = walkoff(lp, ss, mp), walkoff(2 * lp, ss, ms), walkoff(2 * lp, ss, mi)
                            rmax = max(rp, rs, ri) * 1e3
                            rows.append(dict(assign=''.join(assign), target=name, phi_deg=phi, theta_deg=tz,
                                             pump_nm=lp * 1e3, pair_um=2 * lp, deff_pm_V=de, abs_deff=abs(de),
                                             d31_fs_mm=d31, d32_fs_mm=d32, scale_fs_mm=max(abs(d31), abs(d32)),
                                             rho_pump_mrad=rp * 1e3, rho_signal_mrad=rs * 1e3, rho_idler_mrad=ri * 1e3,
                                             rho_max_mrad=rmax, fom=abs(de) / rmax if rmax > 0 else np.inf))
            prev = cur
    return rows


# ------------------------------------------- nondegenerate 0-deg family at fixed pump
def zero_deg_family_at_pump(lp, phis, assign=('f', 's', 'f'), theta_step=0.25):
    """Points (theta, lambda_s) with dk=0 and g0=0 for a fixed pump, per phi."""
    mp, ms, mi = assign
    out = []
    thetas = np.arange(0.0, 180.0 + 1e-9, theta_step)
    for phi in phis:
        prev = None
        for th in thetas:
            s = direction(th, phi)
            # signal roots of dk=0 for ls in (lp, 2lp]  (signal = shorter daughter)
            lss = np.linspace(lp * 1.02, 2 * lp, 60)
            f = np.array([dk(lp, l, s, mp, ms, mi) for l in lss])
            roots = []
            for i in range(len(lss) - 1):
                if f[i] * f[i + 1] < 0:
                    roots.append(brentq(lambda l: dk(lp, l, s, mp, ms, mi), lss[i], lss[i + 1], xtol=1e-9))
            cur = [(ls, gvm(lp, ls, s, mp, ms, mi)[0]) for ls in roots]
            if prev is not None and len(prev) == len(cur) and cur:
                for b, ((ls0, g0p), (ls1, g0c)) in enumerate(zip(prev, cur)):
                    if abs(ls1 - ls0) > 0.1:
                        continue
                    if g0p * g0c < 0:
                        def G(t):
                            ss = direction(t, phi)
                            ll = np.linspace(lp * 1.02, 2 * lp, 60)
                            ff = np.array([dk(lp, l, ss, mp, ms, mi) for l in ll])
                            rr = [brentq(lambda l: dk(lp, l, ss, mp, ms, mi), ll[i], ll[i + 1], xtol=1e-9)
                                  for i in range(len(ll) - 1) if ff[i] * ff[i + 1] < 0]
                            if len(rr) != len(cur):
                                return np.nan
                            return gvm(lp, rr[b], ss, mp, ms, mi)[0]
                        try:
                            tz = brentq(G, th - theta_step, th, xtol=1e-4)
                        except Exception:
                            continue
                        ss = direction(tz, phi)
                        ll = np.linspace(lp * 1.02, 2 * lp, 60)
                        ff = np.array([dk(lp, l, ss, mp, ms, mi) for l in ll])
                        rr = [brentq(lambda l: dk(lp, l, ss, mp, ms, mi), ll[i], ll[i + 1], xtol=1e-9)
                              for i in range(len(ll) - 1) if ff[i] * ff[i + 1] < 0]
                        ls = rr[b]; li = 1 / (1 / lp - 1 / ls)
                        g0, g45, g90, d31, d32 = gvm(lp, ls, ss, mp, ms, mi)
                        de = deff(lp, ls, li, ss, mp, ms, mi)
                        rmax = max(walkoff(lp, ss, mp), walkoff(ls, ss, ms), walkoff(li, ss, mi)) * 1e3
                        out.append(dict(pump_nm=lp * 1e3, phi_deg=phi, theta_deg=tz, lambda_s_um=ls, lambda_i_um=li,
                                        d31_fs_mm=d31, abs_deff=abs(de), rho_max_mrad=rmax, fom=abs(de) / rmax))
            prev = cur
    return out


# --------------------------------------------------------------- JSA purity
def purity(lp, ls, s, assign, L_mm, pump_fwhm_nm, W_sigma=6.0, pix=0.5):
    mp, ms, mi = assign
    li = 1 / (1 / lp - 1 / ls)
    wp0 = 2 * np.pi * C_UM_PS / lp; ws0 = 2 * np.pi * C_UM_PS / ls; wi0 = 2 * np.pi * C_UM_PS / li   # rad/ps
    fwhm_w = 2 * np.pi * C_UM_PS * (pump_fwhm_nm * 1e-3) / lp ** 2
    sig = fwhm_w / (2 * np.sqrt(np.log(2)))
    half = W_sigma * sig; n = int(2 * half / pix) + 1
    ws = np.linspace(ws0 - half, ws0 + half, n); wi = np.linspace(wi0 - half, wi0 + half, n)
    # exact dk on the grid
    def kk(w, mode):
        lam = 2 * np.pi * C_UM_PS / w
        return np.array([n_mode(l, s, mode) for l in lam]) * w / C_UM_PS   # rad/um
    ks = kk(ws, ms); ki = kk(wi, mi)
    WS, WI = np.meshgrid(ws, wi, indexing='ij')
    WP = WS + WI
    lamp = 2 * np.pi * C_UM_PS / WP.ravel()
    npv = np.array([n_mode(l, s, mp) for l in np.unique(np.round(lamp, 9))])
    # pump index varies little; interpolate over the unique set
    ul = np.unique(np.round(lamp, 9))
    kp = np.interp(lamp, ul, npv * (2 * np.pi / ul)).reshape(WP.shape)
    DK = kp - ks[:, None] - ki[None, :]          # rad/um
    L = L_mm * 1e3
    x = DK * L / 2
    f = np.exp(-((WP - wp0) ** 2) / (2 * sig ** 2)) * np.sinc(x / np.pi) * np.exp(1j * x)
    sv = np.linalg.svd(f, compute_uv=False)
    lam_n = sv ** 2 / np.sum(sv ** 2)
    return float(np.sum(lam_n ** 2))


if __name__ == '__main__':
    validate()
    phis = np.arange(0.0, 90.0 + 1e-9, 5.0)
    allrows = []
    for assign in (('f', 's', 'f'), ('s', 'f', 's')):
        print(f"\n== tracing degenerate GVM loci for {''.join(assign)} over phi = 0..90 ==")
        rows = trace_degenerate_loci(assign, phis)
        allrows += rows
        for r in rows:
            print(f"  {r['assign']} {r['target']:>2}deg phi={r['phi_deg']:5.1f} theta={r['theta_deg']:7.3f} lp={r['pump_nm']:8.2f} nm "
                  f"|deff|={r['abs_deff']:.3f} rho_max={r['rho_max_mrad']:.1f} scale={r['scale_fs_mm']:.1f} fom={r['fom']:.4f}")
    with open(os.path.join(HERE, 'offplane_degenerate_loci.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(allrows[0].keys())); w.writeheader(); w.writerows(allrows)

    print("\n== 0-deg family of (f,s,f) at fixed pump vs phi ==")
    fam = []
    for lpnm in (700.0, 722.0, 750.0, 780.0, 800.0):
        r = zero_deg_family_at_pump(lpnm * 1e-3, np.arange(0.0, 90.0 + 1e-9, 5.0))
        fam += r
        for q in r:
            print(f"  lp={q['pump_nm']:.0f} phi={q['phi_deg']:5.1f} theta={q['theta_deg']:7.3f} ls={q['lambda_s_um']:.4f} li={q['lambda_i_um']:.4f} "
                  f"D31={q['d31_fs_mm']:.1f} |deff|={q['abs_deff']:.3f} rho={q['rho_max_mrad']:.1f} fom={q['fom']:.4f}")
    if fam:
        with open(os.path.join(HERE, 'offplane_zero_deg_fixed_pump.csv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(fam[0].keys())); w.writeheader(); w.writerows(fam)
