"""Independent re-verification of BIBO GVM / JSA manuscript numbers.

Written from scratch with numpy/scipy only (no bibojsa import).

Units: wavelength lambda in um, angular frequency omega in rad/ps,
c = 299.792458 um/ps, k = n*omega/c in rad/um, k' = dk/domega in ps/um.
Group-delay differences are reported in fs/mm: ps/um * 1e6 = fs/mm.

Conventions chosen (stated explicitly, see report):
  * x-z principal plane, theta measured from z.  'o' wave is polarized
    along y (n_o = n_y, theta-independent).  'e' wave is polarized in the
    x-z plane with 1/n_e^2 = cos^2(theta)/n_x^2 + sin^2(theta)/n_z^2
    (theta = 0 -> k || z -> polarization along x -> n_e = n_x).
  * LPR2010 sequence (lambda1, lambda2, lambda3) = (idler, signal, pump),
    so 'oeo' = idler o, signal e, pump o; 'ooe' = daughters o, pump e;
    'eeo' = daughters e, pump o.  (Same as bibojsa examples: cfg_for maps
    seq[0] -> idler, seq[1] -> signal, seq[2] -> pump.)
  * Delta_32 = k_p' - k_s', Delta_31 = k_p' - k_i'.
    theta_PM = 0  <=> k_p' = k_s' ; 90 <=> k_p' = k_i' ; 45 <=> k_p' = mean.
  * Pump amplitude alpha = exp(-(w_s+w_i-w_p0)^2 / (2 sigma^2)), so
    |alpha|^2 = exp(-x^2/sigma^2): sigma is the 1/e half-width of the
    pump INTENSITY; intensity FWHM = 2 sqrt(ln 2) sigma.
  * JSA f = alpha * sinc(dk L/2) * exp(i dk L/2)  (also tested without
    phase, and with |f|).  dk = k_p(w_s+w_i) - k_s(w_s) - k_i(w_i) exact.
  * Square window: w_s in [w_s0 - W sigma, w_s0 + W sigma], same for w_i,
    uniform grid, pixel <= 0.5 rad/ps (also 1.0 for comparison).
"""
import numpy as np
from scipy.optimize import brentq, fsolve

C = 299.792458  # um/ps

# ---- Sellmeier (Miyata-Umemura-Kato as in LPR2010 Table 2) --------------
SELL = {
    "x": (3.0759, 0.03169, 0.03323, 0.01402),
    "y": (3.1698, 0.03666, 0.03599, 0.01819),
    "z": (3.6546, 0.05116, 0.03713, 0.02299),
}
D12, D13, D14 = 3.2, -1.76, 1.66  # pm/V


def n2(ax, lam):
    A, B, Cc, D = SELL[ax]
    return A + B / (lam**2 - Cc) - D * lam**2


def dn2_dlam(ax, lam):
    A, B, Cc, D = SELL[ax]
    return -2.0 * B * lam / (lam**2 - Cc) ** 2 - 2.0 * D * lam


def n_ax(ax, lam):
    return np.sqrt(n2(ax, lam))


def dn_ax(ax, lam):
    return dn2_dlam(ax, lam) / (2.0 * n_ax(ax, lam))


def n_e(lam, th):
    """e-wave index in the x-z plane, theta from z."""
    u = np.cos(th) ** 2 / n2("x", lam) + np.sin(th) ** 2 / n2("z", lam)
    return 1.0 / np.sqrt(u)


def dn_e_dlam(lam, th):
    u = np.cos(th) ** 2 / n2("x", lam) + np.sin(th) ** 2 / n2("z", lam)
    du = (-np.cos(th) ** 2 * dn2_dlam("x", lam) / n2("x", lam) ** 2
          - np.sin(th) ** 2 * dn2_dlam("z", lam) / n2("z", lam) ** 2)
    return -0.5 * u ** (-1.5) * du


def dn_e_dth(lam, th):
    u = np.cos(th) ** 2 / n2("x", lam) + np.sin(th) ** 2 / n2("z", lam)
    du = 2 * np.sin(th) * np.cos(th) * (1.0 / n2("z", lam) - 1.0 / n2("x", lam))
    return -0.5 * u ** (-1.5) * du


def index(pol, lam, th):
    return n_ax("y", lam) if pol == "o" else n_e(lam, th)


def dindex(pol, lam, th):
    return dn_ax("y", lam) if pol == "o" else dn_e_dlam(lam, th)


def k(pol, lam, th):
    """k = n*omega/c = 2 pi n / lam  [rad/um]."""
    return 2 * np.pi * index(pol, lam, th) / lam


def kprime(pol, lam, th):
    """dk/domega = n_g / c, n_g = n - lam dn/dlam  [ps/um]."""
    ng = index(pol, lam, th) - lam * dindex(pol, lam, th)
    return ng / C


def kprime_numeric(pol, lam, th, h=1e-3):
    """Richardson-extrapolated central difference on k(omega) as a check."""
    w = 2 * np.pi * C / lam

    def kk(ww):
        return k(pol, 2 * np.pi * C / ww, th)
    d1 = (kk(w + h) - kk(w - h)) / (2 * h)
    d2 = (kk(w + h / 2) - kk(w - h / 2)) / h
    return (4 * d2 - d1) / 3


FS_MM = 1e6  # ps/um -> fs/mm


def to_um(lam):
    return lam


def w_of(lam):
    return 2 * np.pi * C / lam


def lam_of(w):
    return 2 * np.pi * C / w


# ---- phase matching -------------------------------------------------------
def dk_collinear(th, lp, ls, li, seq):
    pi_, ps_, pp_ = seq[0], seq[1], seq[2]
    return k(pp_, lp, th) - k(ps_, ls, th) - k(pi_, li, th)


def pm_angle(lp, ls, li, seq):
    f = lambda th: dk_collinear(th, lp, ls, li, seq)
    ths = np.linspace(1e-4, np.pi / 2 - 1e-4, 400)
    v = np.array([f(t) for t in ths])
    idx = np.where(np.diff(np.sign(v)) != 0)[0]
    if len(idx) == 0:
        return np.nan
    return brentq(f, ths[idx[0]], ths[idx[0] + 1], xtol=1e-14)


def gvm_triplet(lp, ls, li, th, seq):
    pi_, ps_, pp_ = seq[0], seq[1], seq[2]
    kp = kprime(pp_, lp, th)
    ks = kprime(ps_, ls, th)
    ki = kprime(pi_, li, th)
    return kp, ks, ki


def d32_d31(lp, ls, li, th, seq):
    kp, ks, ki = gvm_triplet(lp, ls, li, th, seq)
    return (kp - ks) * FS_MM, (kp - ki) * FS_MM


def deg_condition(lp, seq, which):
    """Degenerate GVM mismatch (fs/mm) as a function of pump wavelength."""
    ld = 2 * lp
    th = pm_angle(lp, ld, ld, seq)
    if np.isnan(th):
        return np.nan, np.nan
    kp, ks, ki = gvm_triplet(lp, ld, ld, th, seq)
    if which == "0":
        m = kp - ks
    elif which == "90":
        m = kp - ki
    else:
        m = kp - 0.5 * (ks + ki)
    return m * FS_MM, th


def solve_deg(seq, which, lo, hi, n=200):
    lps = np.linspace(lo, hi, n)
    vals = np.array([deg_condition(l, seq, which)[0] for l in lps])
    ok = np.isfinite(vals)
    roots = []
    for j in range(n - 1):
        if ok[j] and ok[j + 1] and np.sign(vals[j]) != np.sign(vals[j + 1]):
            r = brentq(lambda l: deg_condition(l, seq, which)[0], lps[j], lps[j + 1], xtol=1e-12)
            roots.append((r, deg_condition(r, seq, which)[1]))
    return roots


# ---- JSA ------------------------------------------------------------------
def sigma_from_fwhm_nm(fwhm_nm, lp_um):
    dw = 2 * np.pi * C / lp_um**2 * (fwhm_nm * 1e-3)  # intensity FWHM rad/ps
    return dw / (2 * np.sqrt(np.log(2))), dw


def jsa_grid(lp, ls0, li0, th, seq, L_mm, sigma, W, pixel, phase=True):
    L = L_mm * 1e3  # um
    ws0, wi0, wp0 = w_of(ls0), w_of(li0), w_of(lp)
    half = W * sigma
    n = int(round(2 * half / pixel)) + 1
    ws = np.linspace(ws0 - half, ws0 + half, n)
    wi = np.linspace(wi0 - half, wi0 + half, n)
    WS, WI = np.meshgrid(ws, wi, indexing="ij")
    pi_, ps_, pp_ = seq[0], seq[1], seq[2]
    dk = (k(pp_, lam_of(WS + WI), th) - k(ps_, lam_of(WS), th) - k(pi_, lam_of(WI), th))
    x = dk * L / 2
    pmf = np.sinc(x / np.pi)
    if phase:
        pmf = pmf * np.exp(1j * x)
    alpha = np.exp(-(WS + WI - wp0) ** 2 / (2 * sigma**2))
    return ws, wi, alpha * pmf


def purity(f):
    s = np.linalg.svd(f, compute_uv=False)
    lam = s**2 / np.sum(s**2)
    return float(np.sum(lam**2))


def filtered_purity(f, wi, wi0, T_target, sigma_pump):
    def trans(fw):
        sg = fw / (2 * np.sqrt(np.log(2)))
        filt = np.exp(-(wi - wi0) ** 2 / (2 * sg**2))
        ff = f * filt[None, :]
        return np.sum(np.abs(ff) ** 2) / np.sum(np.abs(f) ** 2), ff
    fw = brentq(lambda x: trans(x)[0] - T_target, 1e-3 * sigma_pump, 100 * sigma_pump, xtol=1e-10)
    T, ff = trans(fw)
    return fw, purity(ff), T


def fit_tail(Ws, Ps):
    A = np.vstack([np.ones(len(Ws)), 1.0 / np.asarray(Ws)]).T
    coef, *_ = np.linalg.lstsq(A, np.asarray(Ps), rcond=None)
    resid = A @ coef - np.asarray(Ps)
    return coef[0], coef[1], resid


def main():
    out = []
    P = out.append

    # ---- sanity: analytic vs numeric k' ------------------------------
    for pol, lam, th in [("o", 0.6822, 0.8), ("e", 1.3645, 0.8)]:
        a, b = kprime(pol, lam, th), kprime_numeric(pol, lam, th)
        P(f"[check] k' analytic {a:.12e} vs numeric {b:.12e} rel {abs(a-b)/a:.1e} ({pol},{lam})")
    P(f"[check] n at 632.8 nm: nx={n_ax('x',0.6328):.5f} ny={n_ax('y',0.6328):.5f} nz={n_ax('z',0.6328):.5f}"
      " (LPR2010 Table 3: 1.77668, 1.80641, 1.94582)")

    # ---- (a) degenerate oeo loci ------------------------------------
    P("\n=== (a) oeo (idler o, signal e, pump o), x-z plane, degenerate ===")
    res_a = {}
    for which, claim in [("0", (682.2, 45.97)), ("45", (869.6, 43.77)), ("90", (1116.8, 46.11))]:
        roots = solve_deg("oeo", which, 0.45, 1.5)
        for r, th in roots:
            d32, d31 = d32_d31(r, 2 * r, 2 * r, th, "oeo")
            P(f"  theta_PM={which:>2}: pump = {r*1e3:.2f} nm, theta = {np.degrees(th):.3f} deg, "
              f"D32={d32:.2f} D31={d31:.2f} fs/mm   [claim {claim[0]} nm, {claim[1]} deg]")
            res_a[which] = (r, th)
    # also state signal/idler pols check: at theta_PM=0 pump(o) matches signal(e)
    lpA, thA = res_a["0"]
    kp, ks, ki = gvm_triplet(lpA, 2 * lpA, 2 * lpA, thA, "oeo")
    P(f"  at A: k_p'(o)={kp*FS_MM:.3f}, k_s'(e)={ks*FS_MM:.3f}, k_i'(o)={ki*FS_MM:.3f} fs/mm*1e-? (ps/um*1e6)")

    # ---- (b) nondegenerate 0-degree member at 722 nm ----------------
    P("\n=== (b) oeo 0-degree family member at 722.0 nm ===")
    lp = 0.722

    def eqs(v):
        ls, th = v
        li = 1.0 / (1.0 / lp - 1.0 / ls)
        return [dk_collinear(th, lp, ls, li, "oeo") * 1e3,
                (kprime("o", lp, th) - kprime("e", ls, th)) * FS_MM]
    sol, info, ier, msg = fsolve(eqs, [1.16, np.radians(41.0)], full_output=True, xtol=1e-13)
    ls, th = sol
    li = 1.0 / (1.0 / lp - 1.0 / ls)
    d32, d31 = d32_d31(lp, ls, li, th, "oeo")
    P(f"  signal = {ls:.4f} um, idler = {li:.4f} um, theta = {np.degrees(th):.3f} deg, "
      f"D32 = {d32:.3f}, D31 = {d31:.2f} fs/mm  (residuals {eqs(sol)})   "
      f"[claim 1.156, 1.924, 41.36 deg, |D31| = 137.2]")
    res_b = (lp, ls, li, th)

    # ---- (c) deff and walk-off ---------------------------------------
    P("\n=== (c) d_eff = d12 cos(theta) (oeo, x-z, between optic axes) and walk-off ===")
    pts = {"A": (res_a["0"][0], 2 * res_a["0"][0], 2 * res_a["0"][0], res_a["0"][1], 2.22, 85.4),
           "B": (res_b[0], res_b[1], res_b[2], res_b[3], 2.40, 84.4),
           "D": (res_a["45"][0], 2 * res_a["45"][0], 2 * res_a["45"][0], res_a["45"][1], 2.31, 83.4),
           "C": (res_a["90"][0], 2 * res_a["90"][0], 2 * res_a["90"][0], res_a["90"][1], 2.22, 82.0)}
    for tag, (lp_, ls_, li_, th_, dclaim, rclaim) in pts.items():
        deff = D12 * np.cos(th_)
        rho_s = np.arctan(dn_e_dth(ls_, th_) / n_e(ls_, th_))
        rho_p = np.arctan(dn_e_dth(lp_, th_) / n_e(lp_, th_))
        P(f"  {tag}: deff = {deff:.3f} pm/V [claim {dclaim}];  walk-off of e-signal at {ls_:.4f} um: "
          f"{abs(rho_s)*1e3:.2f} mrad [claim {rclaim}];  (e-index at pump wavelength would give {abs(rho_p)*1e3:.2f} mrad)")
    # optic axis half-angle at pump A to confirm 'between the optic axes'
    lamp = res_a["0"][0]
    nx, ny, nz = n_ax("x", lamp), n_ax("y", lamp), n_ax("z", lamp)
    Om = np.degrees(np.arcsin((nz / ny) * np.sqrt((ny**2 - nx**2) / (nz**2 - nx**2))))
    P(f"  optic-axis half angle Omega at {lamp*1e3:.1f} nm = {Om:.2f} deg (theta_A = {np.degrees(res_a['0'][1]):.2f} is between Omega and 180-Omega)")

    # ---- (d) triple GVM points, type I --------------------------------
    P("\n=== (d) type-I triple-GVM degenerate points ===")
    res_d = {}
    for seq, claim in [("ooe", (818.7, 11.04)), ("eeo", (896.8, 34.91))]:
        roots = solve_deg(seq, "0", 0.5, 1.4)
        for r, th in roots:
            kp, ks, ki = gvm_triplet(r, 2 * r, 2 * r, th, seq)
            P(f"  {seq}: pump = {r*1e3:.2f} nm, theta = {np.degrees(th):.3f} deg; k_p'-k_s' = {(kp-ks)*FS_MM:.2e}, "
              f"k_p'-k_i' = {(kp-ki)*FS_MM:.2e}, k_s'-k_i' = {(ks-ki)*FS_MM:.2e} fs/mm; "
              f"deff = {abs(D12*np.cos(th)) if seq=='ooe' else abs(D14*np.sin(2*th)):.3f} pm/V   [claim {claim}]")
            res_d[seq] = (r, th)

    # ---- (e) type-I bound on ooe 0-degree family -----------------------
    P("\n=== (e) ooe 0-degree family (k_p' = k_s'), nondegenerate members ===")
    P("   lp(nm)   ls(um)   li(um)  theta(deg)  |kp'-ks'|  |kp'-ki'|  |ks'-ki'|  (fs/mm)")
    lpE = res_d["ooe"][0]
    for lp in np.array([0.70, 0.74, 0.78, 0.80, 0.81, 0.83, 0.85, 0.90, 0.95, 1.0]):
        def eqsI(v, lp=lp):
            ls, th = v
            li = 1.0 / (1.0 / lp - 1.0 / ls)
            return [dk_collinear(th, lp, ls, li, "ooe") * 1e3,
                    (kprime("e", lp, th) - kprime("o", ls, th)) * FS_MM]
        # initial guess: signal slightly below degeneracy
        best = None
        for g in np.linspace(1.02 * lp, 1.98 * lp, 25):
            sol, info, ier, msg = fsolve(eqsI, [g, np.radians(12.0)], full_output=True, xtol=1e-12)
            if ier == 1 and np.max(np.abs(eqsI(sol))) < 1e-6:
                ls, th = sol
                li = 1.0 / (1.0 / lp - 1.0 / ls)
                if lp < ls <= li and 0 < th < np.pi / 2:
                    best = (ls, li, th)
                    break
        if best is None:
            P(f"   {lp*1e3:.1f}  no solution found")
            continue
        ls, li, th = best
        kp, ks, ki = gvm_triplet(lp, ls, li, th, "ooe")
        a, b, c_ = abs(kp - ks) * FS_MM, abs(kp - ki) * FS_MM, abs(ks - ki) * FS_MM
        P(f"   {lp*1e3:6.1f}  {ls:.4f}  {li:.4f}  {np.degrees(th):8.3f}   {a:8.3f}   {b:8.3f}   {c_:8.3f}   bound ok={max(a,b) <= c_+1e-9}")

    # ---- (f) JSA purity vs window at A and D ---------------------------
    P("\n=== (f) Schmidt purity on square windows (pixel 0.5 rad/ps) ===")
    Wlist = [6, 8, 12, 16, 24]
    cases = {
        "A(mine)": (res_a["0"][0], 2 * res_a["0"][0], 2 * res_a["0"][0], res_a["0"][1], 10.0, 8.0,
                    {6: 0.9488, 12: 0.9455, 24: 0.9439}, (0.942, 0.039)),
        "A(author params 682.241nm/45.96827deg)": (0.682241, 1.364482, 1.364482, np.radians(45.96827), 10.0, 8.0,
                    {6: 0.9488, 12: 0.9455, 24: 0.9439}, (0.942, 0.039)),
        "D(mine)": (res_a["45"][0], 2 * res_a["45"][0], 2 * res_a["45"][0], res_a["45"][1], 2.0, 12.0,
                    {6: 0.833, 24: 0.806}, (0.798, 0.21)),
    }
    fA = None
    for tag, (lp_, ls_, li_, th_, L_mm, bw, claims, fitclaim) in cases.items():
        sigma, fwhm_w = sigma_from_fwhm_nm(bw, lp_)
        P(f"  {tag}: L={L_mm} mm, pump FWHM {bw} nm -> intensity FWHM {fwhm_w:.2f} rad/ps, sigma = {sigma:.3f} rad/ps")
        Ps, Ps_abs, Ps_nophase = [], [], []
        for W in Wlist:
            ws, wi, f = jsa_grid(lp_, ls_, li_, th_, "oeo", L_mm, sigma, W, 0.5, phase=True)
            p = purity(f)
            pa = purity(np.abs(f))
            _, _, f0 = jsa_grid(lp_, ls_, li_, th_, "oeo", L_mm, sigma, W, 0.5, phase=False)
            p0 = purity(f0)
            Ps.append(p); Ps_abs.append(pa); Ps_nophase.append(p0)
            cl = claims.get(W, "")
            P(f"    W={W:2d}sigma  n={len(ws):5d}  P(complex, phase)={p:.5f}  P(|f|)={pa:.5f}  P(real sinc, no phase)={p0:.5f}   claim {cl}")
            if tag.startswith("A") and W == 6:
                fA = (ws, wi, f, sigma, w_of(li_), fwhm_w)
        pinf, Cc, resid = fit_tail(Wlist, Ps)
        P(f"    fit P = P_inf + C/W (complex f): P_inf = {pinf:.4f}, C = {Cc:.4f}, max|resid| = {np.max(np.abs(resid)):.1e}   claim {fitclaim}")
        pinf2, Cc2, resid2 = fit_tail(Wlist, Ps_abs)
        P(f"    fit on |f|:                      P_inf = {pinf2:.4f}, C = {Cc2:.4f}, max|resid| = {np.max(np.abs(resid2)):.1e}")
        # pixel check at W=6 and 24 with 1.0 rad/ps
        for W in (6, 24):
            _, _, f1 = jsa_grid(lp_, ls_, li_, th_, "oeo", L_mm, sigma, W, 1.0, phase=True)
            P(f"    pixel 1.0 rad/ps, W={W}: P = {purity(f1):.5f}")
        if tag == "A(mine)":
            _, _, f1 = jsa_grid(lp_, ls_, li_, th_, "oeo", L_mm, sigma, 6, 0.25, phase=True)
            P(f"    pixel 0.25 rad/ps, W=6: P = {purity(f1):.5f}")

    # ---- (g) heralding filter at A -------------------------------------
    P("\n=== (g) Gaussian idler filter at A (6-sigma window, pixel 0.5) ===")
    ws, wi, f, sigma, wi0, fwhm_w = fA
    liA = lam_of(wi0)
    for T in (0.9, 0.7):
        fw, p, Tach = filtered_purity(f, wi, wi0, T, sigma)
        fw_nm = liA**2 / (2 * np.pi * C) * fw * 1e3
        P(f"  T={T}: filter intensity FWHM = {fw:.3f} rad/ps = {fw_nm:.2f} nm, purity = {p:.4f} (T achieved {Tach:.4f})"
          f"   [claim {'10.7 nm, 0.988' if T==0.9 else '3.6 nm, 0.994'}]")
    # window independence check at 16 sigma
    lp_, ls_, li_, th_ = res_a["0"][0], 2 * res_a["0"][0], 2 * res_a["0"][0], res_a["0"][1]
    ws16, wi16, f16 = jsa_grid(lp_, ls_, li_, th_, "oeo", 10.0, sigma, 16, 0.5)
    for T in (0.9, 0.7):
        fw, p, _ = filtered_purity(f16, wi16, w_of(li_), T, sigma)
        P(f"  (16-sigma window) T={T}: FWHM = {liA**2/(2*np.pi*C)*fw*1e3:.2f} nm, purity = {p:.4f}")

    # ---- (h) first sinc zeros ------------------------------------------
    P("\n=== (h) first sinc zero ===")
    lpA, thA = res_a["0"]
    lsA = 2 * lpA
    d32, d31 = d32_d31(lpA, lsA, lsA, thA, "oeo")
    L = 10.0
    lin = 2 * np.pi / (abs(d31) * 1e-15 * L) * 1e-12
    P(f"  A: |D31| = {abs(d31):.2f} fs/mm -> linear estimate 2pi/(D31 L) = {lin:.3f} rad/ps [claim 3.7]")
    # exact: scan idler detuning at fixed signal
    wsA, wiA = w_of(lsA), w_of(lsA)
    dws = np.linspace(-40, 40, 400001)
    dk_i = (k("o", lam_of(wsA + wiA + dws), thA) - k("e", lam_of(wsA), thA) - k("o", lam_of(wiA + dws), thA)) * L * 1e3 / 2
    zc = dws[np.where(np.diff(np.sign(dk_i - np.pi)) != 0)[0]]
    zc2 = dws[np.where(np.diff(np.sign(dk_i + np.pi)) != 0)[0]]
    P(f"  A exact (Sellmeier) idler-axis zeros of sinc: dk L/2 = +pi at d_omega_i = {zc}, -pi at {zc2} rad/ps")
    lpD, thD = res_a["45"]
    lsD = 2 * lpD
    d32D, d31D = d32_d31(lpD, lsD, lsD, thD, "oeo")
    LD = 2.0
    linD = 2 * np.pi / (abs(d31D - d32D) * 1e-15 * LD) * 1e-12
    P(f"  D: D32 = {d32D:.2f}, D31 = {d31D:.2f} fs/mm -> linear antidiagonal zero 2pi/((D31-D32) L) = {linD:.3f} rad/ps [claim 21.5]")
    wsD = w_of(lsD)
    dk_d = (k("o", lam_of(2 * wsD), thD) - k("e", lam_of(wsD + dws), thD) - k("o", lam_of(wsD - dws), thD)) * LD * 1e3 / 2
    zcD = dws[np.where(np.diff(np.sign(np.abs(dk_d) - np.pi)) != 0)[0]]
    P(f"  D exact antidiagonal zeros (|dk L/2| = pi) at d_omega = {zcD} rad/ps")
    sigD, fwD = sigma_from_fwhm_nm(12.0, lpD)
    P(f"  D pump intensity FWHM = {fwD:.2f} rad/ps [claim 29.9]")

    txt = "\n".join(out)
    print(txt)
    with open(__file__.replace("independent.py", "independent_output.txt"), "w") as fh:
        fh.write(txt)


if __name__ == "__main__":
    main()
