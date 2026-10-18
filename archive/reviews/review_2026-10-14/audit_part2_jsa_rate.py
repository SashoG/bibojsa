"""Audit part 2: JSA/purity conventions, filter arm, pair rate, off-plane reduction,
and the manuscript's derived arithmetic."""
import sys, os, csv
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from scipy.special import erf

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "bibojsa"))
sys.path.insert(0, os.path.join(HERE, "..", "bibojsa", "examples"))
sys.path.insert(0, os.path.join(HERE, "..", "manuscript", "exploration_2026-09-21"))
from bibojsa.materials.bibo import BIBO
from bibojsa import jsa as J
from bibojsa.phasematch import ThreeWaveConfig
from bibojsa.sellmeier import lambda_um_to_omega, C_LIGHT
from landscape_points import adaptive_axes, cfg_for, run_case

X, Y, Z = BIBO["x"], BIBO["y"], BIBO["z"]
C_UM_PS = 299.792458
out = []
P = lambda *a: (print(*a), out.append(" ".join(str(x) for x in a)))

# ------------------------------------------------------------ 1. pump sigma
P("=== 1. Pump-envelope convention ===")
fw = J.fwhm_omega_from_fwhm_wavelength(0.008, 0.682241)
sig = J.gaussian_pump_sigma(fw)
P(f"  8-nm pump at 682.241 nm: intensity FWHM = {fw*1e-12:.3f} rad/ps; sigma = {sig*1e-12:.3f} rad/ps (paper: 19.4); FWHM/sigma = {fw/sig:.4f} = 2 sqrt(ln2) = {2*np.sqrt(np.log(2)):.4f}")
# intensity at Omega = sigma is 1/e?
Om = np.array([sig])
a = J.pump_envelope_gaussian(Om, np.zeros(1), 0.0, sig)
P(f"  |alpha(Omega=sigma)|^2 = {float(a[0]**2):.5f} (1/e = {np.exp(-1):.5f}) -> sigma is the 1/e half-width of the INTENSITY, as the paper says")
a2 = J.pump_envelope_gaussian(np.array([fw / 2]), np.zeros(1), 0.0, sig)
P(f"  |alpha(FWHM/2)|^2 = {float(a2[0]**2):.5f} (0.5 expected)")
# exact vs linear Jacobian for 8 nm at 682 nm
w1 = lambda_um_to_omega(0.682241 - 0.004); w2 = lambda_um_to_omega(0.682241 + 0.004)
P(f"  exact omega-width of an 8-nm band: {(w1-w2)*1e-12:.4f} rad/ps vs linear Jacobian {fw*1e-12:.4f} (rel {abs((w1-w2)/fw-1):.1e})")
# sech2 FWHM constant
P(f"  sech^2 half-max argument arccosh(sqrt2) = {np.arccosh(np.sqrt(2)):.6f} (code 0.881374; 2x = {2*np.arccosh(np.sqrt(2)):.6f} vs 1.762747)")
# filter
sigf = J.gaussian_pump_sigma(1.0)
filt = np.exp(-(0.5)**2 / (2 * sigf**2))
P(f"  heralding filter with intensity FWHM 1: amplitude transmission at +FWHM/2 = {filt:.5f}, intensity {filt**2:.5f} (0.5 expected) -> filter FWHM is an INTENSITY FWHM; applied as amplitude; 'fraction of pairs' = sum|f filt|^2 / sum|f|^2 (checked in landscape_points.run_case)")
# sinc
x = np.array([0.3, 1.0, np.pi, 4.0])
P(f"  np.sinc(x/pi) vs sin(x)/x max diff: {np.max(np.abs(np.sinc(x/np.pi)-np.sin(x)/x)):.1e}")

# -------------------------------------------- 2. Schmidt purity: grid independence
P("\n=== 2. Schmidt decomposition on the discrete grid ===")
rng = np.random.default_rng(0)
# analytic test: separable-ish kernel with known purity: f = exp(-(a x^2 + b y^2 + 2 c x y)); purity = sqrt(1 - c^2/(ab))... use two-mode Schmidt test instead
xs = np.linspace(-6, 6, 301); ys = np.linspace(-6, 6, 451)
Xg, Yg = np.meshgrid(xs, ys, indexing="ij")
for c in (0.0, 0.3, 0.6):
    f = np.exp(-(Xg**2 + Yg**2 + 2 * c * Xg * Yg) / 2)
    # exact: for a correlated Gaussian amplitude exp(-(x^2+y^2+2cxy)/2), the purity is sqrt(1-c^2)  (Schmidt number K = 1/sqrt(1-c^2))
    _, p, K = J.schmidt_decompose(f)
    # same with different (coarser, anisotropic) grid
    xs2 = np.linspace(-6, 6, 121); ys2 = np.linspace(-6, 6, 701)
    X2, Y2 = np.meshgrid(xs2, ys2, indexing="ij")
    _, p2, _ = J.schmidt_decompose(np.exp(-(X2**2 + Y2**2 + 2 * c * X2 * Y2) / 2))
    P(f"  Gaussian test c={c}: purity {p:.6f} (grid 301x451), {p2:.6f} (121x701); exact sqrt(1-c^2) = {np.sqrt(1-c**2):.6f}")
P("  -> unweighted SVD of samples on a uniform (even anisotropic) grid gives the continuum purity: grid spacings enter as a global factor that the normalisation removes.")

# -------------------------------------- 3. featured points: own JSA vs package; filter arm at C
P("\n=== 3. Point A/C: package purity on the adaptive window; filter on idler vs signal arm ===")
FEAT = {"A": ("xz", "oeo", 682.241, 1.364482, 1.364482, 45.96827, 10.0, 8.0),
        "B": ("xz", "oeo", 722.000, 1.155982, 1.923690, 41.35980, 10.0, 12.0),
        "C": ("xz", "oeo", 1116.841, 2.233682, 2.233682, 46.10995, 10.0, 16.0),
        "D": ("xz", "oeo", 869.583, 1.739166, 1.739166, 43.77035, 2.0, 12.0)}


def filtered(res, axis_vals, w0, fw_pump, herald_axis, target):
    tot = float(np.sum(np.abs(res.f)**2))

    def trans(k):
        ff = J.apply_heralding_filter_gaussian(res.f, axis_vals, w0, k * fw_pump, herald_axis)
        return float(np.sum(np.abs(ff)**2)) / tot
    k = brentq(lambda k: trans(k) - target, 0.01, 60.0)
    ff = J.apply_heralding_filter_gaussian(res.f, axis_vals, w0, k * fw_pump, herald_axis)
    return J.schmidt_decompose(ff)[1], k * fw_pump


for tag in ("A", "B", "C", "D"):
    plane, seq, pnm, ls, li, th, Lmm, bw = FEAT[tag]
    pt = {"plane": plane, "interaction": seq, "pump_nm": pnm, "lambda_s_um": ls, "lambda_i_um": li, "angle_deg": th}
    res, edge, filt, wsa, wia, wp0, ws0, wi0 = run_case(pt, Lmm, bw, n=200)
    fwp = J.fwhm_omega_from_fwhm_wavelength(bw * 1e-3, pnm / 1000)
    # marginal widths (intensity) along each axis
    I = np.abs(res.f)**2
    ms = I.sum(axis=1); mi = I.sum(axis=0)

    def fwhm(axis, m):
        half = m.max() / 2
        idx = np.where(m >= half)[0]
        return (axis[idx[-1]] - axis[idx[0]]) * 1e-12
    P(f"  {tag}: P = {res.purity:.5f}, edge = {edge:.4f}; marginal intensity FWHM signal {fwhm(wsa, ms):.2f} rad/ps, idler {fwhm(wia, mi):.2f} rad/ps; "
      f"window half-widths {(wsa[-1]-wsa[0])/2*1e-12:.1f}/{(wia[-1]-wia[0])/2*1e-12:.1f} rad/ps = {(wsa[-1]-wsa[0])/2*(ls*1e-6)**2/(2*np.pi*C_LIGHT)*1e9:.1f}/{(wia[-1]-wia[0])/2*(li*1e-6)**2/(2*np.pi*C_LIGHT)*1e9:.1f} nm")
    for T in (0.9, 0.7):
        pI, fI = filtered(res, wia, wi0, fwp, 1, T)
        pS, fS = filtered(res, wsa, ws0, fwp, 0, T)
        nmI = fI * (li * 1e-6)**2 / (2 * np.pi * C_LIGHT) * 1e9
        nmS = fS * (ls * 1e-6)**2 / (2 * np.pi * C_LIGHT) * 1e9
        P(f"      T={T:.1f}: filter on IDLER: P={pI:.4f} (FWHM {nmI:.1f} nm);  filter on SIGNAL: P={pS:.4f} (FWHM {nmS:.1f} nm)")

# ------------------------------------------------------- 4. independent JSA at A (own code)
P("\n=== 4. Independent JSA at A (analytic Sellmeier, own grid) ===")
SELL = {"x": (3.0759, 0.03169, 0.03323, 0.01402), "y": (3.1698, 0.03666, 0.03599, 0.01819), "z": (3.6546, 0.05116, 0.03713, 0.02299)}


def n2(ax, lam):
    A, Bc, Cc, D = SELL[ax]
    return A + Bc / (lam**2 - Cc) - D * lam**2


def nxz(pol, lam, th):
    return np.sqrt(n2("y", lam)) if pol == "o" else 1 / np.sqrt(np.cos(th)**2 / n2("x", lam) + np.sin(th)**2 / n2("z", lam))


def own_jsa(lp, ls, li, thd, Lmm, bw_nm, Wsig, pix=0.5):
    th = np.radians(thd)
    wp0, ws0, wi0 = (2 * np.pi * C_UM_PS / v for v in (lp, ls, li))
    fw = 2 * np.pi * C_UM_PS * bw_nm * 1e-3 / lp**2
    sg = fw / (2 * np.sqrt(np.log(2)))
    half = Wsig * sg; n = int(2 * half / pix) + 1
    ws = np.linspace(ws0 - half, ws0 + half, n); wi = np.linspace(wi0 - half, wi0 + half, n)
    WS, WI = np.meshgrid(ws, wi, indexing="ij")
    lam_s, lam_i, lam_p = 2 * np.pi * C_UM_PS / WS, 2 * np.pi * C_UM_PS / WI, 2 * np.pi * C_UM_PS / (WS + WI)
    dk = (nxz("o", lam_p, th) * (WS + WI) - nxz("e", lam_s, th) * WS - nxz("o", lam_i, th) * WI) / C_UM_PS   # rad/um
    x = dk * Lmm * 1e3 / 2
    f = np.exp(-(WS + WI - wp0)**2 / (2 * sg**2)) * np.sinc(x / np.pi) * np.exp(1j * x)
    s = np.linalg.svd(f, compute_uv=False); lam = s**2 / np.sum(s**2)
    return float(np.sum(lam**2)), ws, wi, f, sg


for W in (6, 12, 24):
    p, *_ = own_jsa(0.682241, 1.364482, 1.364482, 45.96827, 10, 8, W)
    P(f"  A, W={W} sigma: own P = {p:.5f} (window_tail.csv: 0.94877/0.94546/0.94386)")
p, *_ = own_jsa(0.869583, 1.739166, 1.739166, 43.77035, 2, 12, 6)
P(f"  D, W=6 sigma: own P = {p:.5f} (window_tail.csv 0.83282)")

# ------------------------------------------------------------- 5. pair rate
P("\n=== 5. Pair rate at A (Ling 2008 Eq. 36 structure, re-derived) ===")
eps0 = 8.8541878128e-12; c = 299792458.0
lp, ls = 0.682241, 1.364482; th = np.radians(45.96827)


def dn_dlam(pol, lam, th, h=1e-5):
    return (nxz(pol, lam + h, th) - nxz(pol, lam - h, th)) / (2 * h)


n_p, n_s, n_i = nxz("o", lp, th), nxz("e", ls, th), nxz("o", ls, th)
ng_s = n_s - ls * dn_dlam("e", ls, th); ng_i = n_i - ls * dn_dlam("o", ls, th)
dng = abs(ng_s - ng_i)
d = 2.224e-12; wp = 2 * np.pi * c / (lp * 1e-6)
P(f"  n_p={n_p:.5f} n_s={n_s:.5f} n_i={n_i:.5f}; group indices n_g,s={ng_s:.5f} n_g,i={ng_i:.5f}; Delta n_g = {dng:.5f} (paper 0.051; script 0.05080); phase-index diff {n_s-n_i:.5f}")


def R_thin(L, W, Pw=1e-3):
    return 4 * d**2 * Pw * L * wp**2 / (9 * np.pi * n_p * n_s * n_i * eps0 * c**2 * W**2 * dng)


P("  derivation check: H_I = 4 eps0 d E_p E_s E_i (consistent with P(w3) = 2 eps0 d A1 A2), Gaussian modes with equal waists,")
P("  transverse overlap (2/(pi W^2)) * pi W^2/3 = 2/3 -> 4/9; P = n_p eps0 c E0^2 pi W^2/4; spectral integral 2 pi c/(L Delta n_g):")
P("  R = 16 d^2 P L w_s w_i /(9 pi eps0 c^2 n_p n_s n_i W^2 Dn_g) = 4 d^2 P L w_p^2/(9 pi ...) at degeneracy  [matches pair_rate.py]")
rho = 85.38e-3
for name, L, W in (("10 mm, 482 um", 10e-3, 482e-6), ("10 mm, 480 um", 10e-3, 480e-6), ("1 mm, 60 um", 1e-3, 60e-6), ("1 mm, 30 um", 1e-3, 30e-6), ("10 mm, 60 um", 10e-3, 60e-6)):
    R = R_thin(L, W)
    # script's walk-off factor: (L_eff/L)^2 with 2-Gaussian overlap exp(-(rho z)^2/(2W^2))
    Leff = quad(lambda z: np.exp(-(rho * z)**2 / (2 * W**2)), 0, L)[0]
    f_script = (Leff / L)**2
    # Parseval (spectrally integrated) factor with the same 2-Gaussian overlap: (1/L) int h^2 dz
    f_pars2 = quad(lambda z: np.exp(-(rho * z)**2 / W**2), 0, L)[0] / L
    # three-beam overlap: pump & o-idler fixed, e-signal displaced: h = exp(-2 (rho z)^2/(3 W^2)); Parseval
    f_pars3 = quad(lambda z: np.exp(-4 * (rho * z)**2 / (3 * W**2)), 0, L)[0] / L
    P(f"  {name:14s}: R_thin = {R:7.0f} /s/mW; walk-off factor: script (L_eff/L)^2 = {f_script:.3f}; Parseval 2-beam = {f_pars2:.3f}; Parseval 3-beam overlap = {f_pars3:.3f}"
      f" -> R = {R*f_script:.0f} / {R*f_pars2:.0f} / {R*f_pars3:.0f}")
P(f"  z_R check (10 mm / 482 um): pi W^2 n/lambda = {np.pi*(482e-6)**2*n_s/(ls*1e-6)*1e3:.0f} mm >> L (thin-crystal limit ok)")
P(f"  KTP comparison: 25.1 kHz/mW detected at 21% coincidence-to-singles -> generated 25.1e3/0.21^2 = {25.1e3/0.21**2*1e-6:.2f} MHz/mW (paper ~0.6)")
P(f"  ratios vs BIBO at A: 0.6e6/25071 = {0.6e6/25071:.0f}, 0.6e6/6268 = {0.6e6/6268:.0f}, 0.6e6/971 = {0.6e6/971:.0f}, 0.6e6/415 = {0.6e6/415:.0f}  (paper '20-1400 times')")

# ---------------------------------------------------------- 6. off-plane
P("\n=== 6. Off-plane script: reduction at phi=0 and the 11-12 % figure ===")
import offplane as o
s0 = o.direction(45.96827, 0.0)
(nf, Df), (ns_, Ds) = o.eigen(1.364482, s0)
P(f"  phi=0, theta=45.968, 1.3645 um: offplane n_fast={nf:.6f} (n_y={np.sqrt(n2('y',1.364482)):.6f}), n_slow={ns_:.6f} (n_e(theta)={nxz('e',1.364482,th):.6f})")
P(f"  offplane deff (f,s,f) at A = {o.deff(lp, ls, ls, s0, 'f','s','f'):.4f} pm/V vs d12 cos(theta) = {3.2*np.cos(th):.4f}")
P(f"  offplane kprime units: kprime(1.3645,'s') = {o.kprime(ls, s0, 's'):.6e} ps/um -> *1e6 = {o.kprime(ls, s0, 's')*1e6:.1f} fs/mm; gvm() multiplies by 1e3 -> ps/mm (CSV label slip confirmed)")
g0, g45, g90, d31, d32 = o.gvm(lp, ls, s0, 'f', 's', 'f')
P(f"  offplane gvm at A: d31 = {d31:.5f} ('fs/mm' column; actually ps/mm = {d31*1e3:.1f} fs/mm), d32 = {d32:.2e}")
rows = list(csv.DictReader(open(os.path.join(HERE, "..", "manuscript", "exploration_2026-09-21", "offplane_degenerate_loci.csv"))))
for tgt, tag in (("0", "A"), ("45", "D"), ("90", "C")):
    rs = [r for r in rows if r["assign"] == "fsf" and r["target"] == tgt]
    inpl = [r for r in rs if float(r["phi_deg"]) == 0 and float(r["theta_deg"]) < 90]
    if not inpl:
        P(f"  {tag}: no phi=0 row in CSV (EPM locus is in offplane_epm_locus.csv)")
        continue
    f0 = float(inpl[0]["fom"])
    oct1 = [r for r in rs if float(r["theta_deg"]) < 90]
    oct2 = [r for r in rs if float(r["theta_deg"]) > 90]
    b1 = max(oct1, key=lambda r: float(r["fom"])); b2 = max(oct2, key=lambda r: float(r["fom"]))
    P(f"  {tag} locus: in-plane fom {f0:.4f}; best theta<90: {float(b1['fom']):.4f} at phi={b1['phi_deg']} (+{100*(float(b1['fom'])/f0-1):.1f} %); best mirror octant {float(b2['fom']):.4f} at phi={b2['phi_deg']} (+{100*(float(b2['fom'])/f0-1):.1f} %);"
      f" scale phi=0 {float(inpl[0]['scale_fs_mm'])*1e3:.1f} -> phi=90 {float([r for r in oct1 if float(r['phi_deg'])==90][0]['scale_fs_mm'])*1e3:.1f} fs/mm")
    r5 = [r for r in rs if float(r["phi_deg"]) == 5]
    P(f"      |deff| at phi=5: {', '.join(f'{float(r[chr(97)+chr(98)+chr(115)+chr(95)+chr(100)+chr(101)+chr(102)+chr(102)]):.3f}' for r in r5)}")
epm = list(csv.DictReader(open(os.path.join(HERE, "..", "manuscript", "exploration_2026-09-21", "offplane_epm_locus.csv"))))
P(f"  EPM locus columns: {list(epm[0].keys())}")
for r in epm[:1] + epm[-1:]:
    P("   ", {k: (round(float(v), 4) if v not in ('', None) and v.replace('.', '', 1).replace('-', '', 1).isdigit() else v) for k, v in r.items()})
try:
    fom0 = float([r for r in epm if abs(float(r["phi_deg"])) < 1e-9][0]["fom"])
    best = max(epm, key=lambda r: float(r["fom"]))
    P(f"  D locus (EPM file): in-plane fom {fom0:.4f}, best {float(best['fom']):.4f} at phi={best['phi_deg']} (+{100*(float(best['fom'])/fom0-1):.1f} %)")
except Exception as e:
    P("  (EPM fom check skipped:", e, ")")

# -------------------------------------------------------- 7. manuscript arithmetic
P("\n=== 7. Manuscript arithmetic re-computed ===")
P(f"  x-z birefringence at 1.364 um: (nz-nx)/nx = {100*(np.sqrt(n2('z',1.3645))/np.sqrt(n2('x',1.3645))-1):.2f} % (paper 'close to 9 %')")
P(f"  aperture length sqrt(pi) w0/rho at 60 um, 85.4 mrad: {np.sqrt(np.pi)*0.060/0.0854:.3f} mm (paper 1.25)")
for Lmm in (0.5, 2, 10, 30):
    P(f"  matched waist rho L/sqrt(pi) at L={Lmm} mm: {0.0854*Lmm/np.sqrt(np.pi)*1e3:.0f} um (paper 24/96/482/1450)")
P(f"  L^-0.6 scaling of optimum pump bandwidth from 8 nm at 10 mm: 15 mm -> {8*(10/15)**0.6:.1f}, 20 -> {8*(10/20)**0.6:.1f}, 30 -> {8*(10/30)**0.6:.1f} nm (paper 6, 6, 4);"
  f" 3 nm optimum at L = {10*(8/3)**(1/0.6):.0f} mm (paper ~50)")
# fit the exponent from the grid optima (8,6,6,4 at 10,15,20,30)
Ls = np.array([10, 15, 20, 30.]); bws = np.array([8, 6, 6, 4.])
slope = np.polyfit(np.log(Ls), np.log(bws), 1)[0]
P(f"  least-squares exponent of (8,6,6,4) vs (10,15,20,30): {slope:.2f} (coarse grid 2/3/4/6/8/12/16 nm)")
P(f"  d_eff/rho table: " + ", ".join(f"{t}:{d/r:.4f}" for t, d, r in (("A", 2.22, 85.4), ("B", 2.40, 84.4), ("C", 2.22, 82.0), ("D", 2.31, 83.4), ("E", 3.14, 30.7), ("F", 1.56, 76.3), ("R", 2.47, 83.1))))
P(f"  E over others: {0.1023/0.0297:.2f} (R) .. {0.1023/0.0260:.2f} (A) -> 'at least 3.4', '3.4-3.9'; squared {0.1023**2/0.0297**2:.1f}-{0.1023**2/0.0260**2:.1f} (paper 12-15)")
P(f"  C ratio symmetric/asymmetric: 0.21/0.039 = {0.21/0.039:.1f}, 0.21/0.038 = {0.21/0.038:.1f}, 0.21/0.056 = {0.21/0.056:.1f} (paper 'four to six')")
P(f"  Jin deff ratios: 2.50/2.22 = {2.50/2.22-1:.3f}, 2.48/2.31 = {2.48/2.31-1:.3f}, 2.30/2.22 = {2.30/2.22-1:.3f} (paper 13, 7, 4 %)")
P(f"  Sellmeier spread: 688.91/682.24-1 = {688.91/682.24-1:.4f}, 879.73/869.58-1 = {879.73/869.58-1:.4f}, 1131.58/1116.84-1 = {1131.58/1116.84-1:.4f} (paper 1.0-1.3 %)")
P(f"  angle spread: max pairwise |dtheta|: 0deg {max(45.968,45.795,45.961)-min(45.968,45.795,45.961):.3f}, 45deg {max(43.770,43.601,43.794)-min(43.770,43.601,43.794):.3f}, 90deg {max(46.110,45.942,46.023)-min(46.110,45.942,46.023):.3f} (paper 'up to 0.19')")
P(f"  KDP: K = 1/0.9199 = {1/0.9199:.4f} (paper 1.087)")
P(f"  Fig.5 panel half-widths E/D = 509/115 = {509/115:.1f}, E/A = 509/70 = {509/70:.1f} (paper 'four times', 'seven times')")
P(f"  10 mm at 1 cm^-1: transmission e^-1 = {np.exp(-1):.3f}")
P(f"  D: antidiagonal first sinc zero 2pi/((D31-D32) L) = {2*np.pi/(145.72e-15/1e-3*2e-3)*1e-12:.2f} rad/ps (paper 21.6); pump FWHM 12 nm at 869.583: {2*np.pi*C_UM_PS*0.012/0.869583**2:.2f} rad/ps (paper 29.9)")
P(f"  A: idler first sinc zero 2pi/(D31 L) = {2*np.pi/(169.44e-15/1e-3*10e-3)*1e-12:.2f} rad/ps (paper 3.7); pump FWHM {2*np.pi*C_UM_PS*0.008/0.682241**2:.2f} (paper 32.4)")
P(f"  filter FWHM at A: 10.7 nm at 1.3645 um = {10.7e-3*2*np.pi*C_UM_PS/1.364482**2:.2f} rad/ps = {10.7e-3*2*np.pi*C_UM_PS/1.364482**2/32.38:.3f} pump FWHM")
P(f"  Kleinman-average check: (2.93+3.48)/2={3.205:.3f}->3.2; (1.93+1.58)/2={1.755:.3f}->1.76; (1.63+1.67+1.67)/3={(1.63+1.67+1.67)/3:.3f}->1.66")
P(f"  (2/pi) d24(KTP~3.6-3.9 pm/V) = {2/np.pi*3.6:.2f}-{2/np.pi*3.9:.2f} pm/V (paper 'comparable' to 2.2-2.4)")
P(f"  Eq. tail sign: tail weight excluded ~ 1/W, P ~ 1 - 2*(excluded weight) -> C>0: consistent with fits C=0.039..0.21")

with open(os.path.join(HERE, "audit_part2_output.txt"), "w") as fh:
    fh.write("\n".join(out))
