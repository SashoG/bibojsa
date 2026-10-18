"""Audit part 1: Sellmeier, derivatives, biaxial optics, GVM, d_eff.

Independent numpy reference (analytic derivatives, general Fresnel equation,
full Voigt-tensor contraction) compared against the bibojsa package.
Run from BIBO_JSA/bibojsa so that `import bibojsa` resolves to the package.
"""
import sys, os
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bibojsa"))
from bibojsa.materials.bibo import BIBO, D_COEFFS_PM_PER_V, BIBO_HELLWIG, BIBO_UMEMURA
from bibojsa import sellmeier as S
from bibojsa import biaxial as B
from bibojsa import gvm as G
from bibojsa import deff as DF
from bibojsa.phasematch import ThreeWaveConfig, find_phase_matching_angle
from bibojsa import landscape as L

C_UM_PS = 299.792458  # um/ps
C_SI = 299792458.0

SELL = {"x": (3.0759, 0.03169, 0.03323, 0.01402),
        "y": (3.1698, 0.03666, 0.03599, 0.01819),
        "z": (3.6546, 0.05116, 0.03713, 0.02299)}

out = []
P = lambda *a: (print(*a), out.append(" ".join(str(x) for x in a)))


# ---------------------------------------------------------------- 1. Sellmeier
def n2(ax, lam):
    A, Bc, Cc, D = SELL[ax]
    return A + Bc / (lam**2 - Cc) - D * lam**2


def dn2(ax, lam):
    A, Bc, Cc, D = SELL[ax]
    return -2 * Bc * lam / (lam**2 - Cc)**2 - 2 * D * lam


def d2n2(ax, lam):
    A, Bc, Cc, D = SELL[ax]
    return (-2 * Bc / (lam**2 - Cc)**2 + 8 * Bc * lam**2 / (lam**2 - Cc)**3 - 2 * D)


def n(ax, lam):
    return np.sqrt(n2(ax, lam))


def dn(ax, lam):
    return dn2(ax, lam) / (2 * n(ax, lam))


def d2n(ax, lam):
    # n = sqrt(n2): n'' = n2''/(2n) - (n2')^2/(4 n^3)
    return d2n2(ax, lam) / (2 * n(ax, lam)) - dn2(ax, lam)**2 / (4 * n(ax, lam)**3)


def ng(ax, lam):
    return n(ax, lam) - lam * dn(ax, lam)


def gvd_ps2_per_um(ax, lam):
    # k'' = lam^3/(2 pi c^2) d2n/dlam2   [ps^2/um] with c in um/ps
    return lam**3 / (2 * np.pi * C_UM_PS**2) * d2n(ax, lam)


P("=== 1. Sellmeier coefficients and derivatives ===")
P("code coefficients vs LPR2010 Table 2 (MUK fit) as transcribed in my reference:")
for ax in "xyz":
    c = BIBO[ax]
    P(f"  {ax}: code (A,B,C,D)=({c.A},{c.B},{c.C},{c.D}) ref={SELL[ax]} match={np.allclose([c.A,c.B,c.C,c.D],SELL[ax])}")
P("LPR2010 Table 3 reference indices:")
ref = {0.6328: (1.77668, 1.80641, 1.94582), 1.064: (1.75752, 1.78400, 1.91711)}
for lam, (rx, ry, rz) in ref.items():
    vals = [n(a, lam) for a in "xyz"]
    P(f"  {lam} um: nx,ny,nz = {vals[0]:.6f} {vals[1]:.6f} {vals[2]:.6f}; |diff| = "
      f"{abs(vals[0]-rx):.1e} {abs(vals[1]-ry):.1e} {abs(vals[2]-rz):.1e}")
# optic axis angle
lam = 0.6328
nx, ny, nz = (n(a, lam) for a in "xyz")
Om = np.degrees(np.arcsin(nz / ny * np.sqrt((ny**2 - nx**2) / (nz**2 - nx**2))))
P(f"  optic axis 2*Omega(632.8 nm) = {2*Om:.3f} deg (LPR2010: 52.6); code: "
  f"{2*S.optic_axis_angle_deg(lam, BIBO['x'], BIBO['y'], BIBO['z']):.3f}")
lam = 1.064
nx, ny, nz = (n(a, lam) for a in "xyz")
Om = np.degrees(np.arcsin(nz / ny * np.sqrt((ny**2 - nx**2) / (nz**2 - nx**2))))
P(f"  optic axis 2*Omega(1064 nm) = {2*Om:.3f} deg")
# Vz formula alternative: tan^2 Vz = (1/nx^2 - 1/ny^2)/(1/ny^2 - 1/nz^2)
Vz = np.degrees(np.arctan(np.sqrt((1/nx**2 - 1/ny**2) / (1/ny**2 - 1/nz**2))))
P(f"  cross-check with tan^2 Vz formula at 1064: 2Vz = {2*Vz:.3f}")

P("\nk' and k'' : code finite differences vs analytic (rel. error)")
for ax in "xyz":
    for lam in (0.4, 0.6822, 1.0, 1.3645, 2.2337, 3.0):
        w = 2 * np.pi * C_SI / (lam * 1e-6)  # rad/s
        kp_code = S.kprime(w, BIBO[ax])             # s/m
        kpp_code = S.kdoubleprime(w, BIBO[ax])      # s^2/m
        kp_ref = ng(ax, lam) / C_SI
        kpp_ref = gvd_ps2_per_um(ax, lam) * 1e-24 / 1e-6   # ps^2/um -> s^2/m
        ngc = S.group_index(w, BIBO[ax])
        P(f"  {ax} {lam:6.4f} um: n_g code {ngc:.9f} ref {ng(ax,lam):.9f} rel {abs(ngc/ng(ax,lam)-1):.1e};"
          f"  k'' code {kpp_code:.6e} ref {kpp_ref:.6e} rel {abs(kpp_code/kpp_ref-1):.1e}")

P("\nunit conventions: omega_to_lambda_um / lambda_um_to_omega round trip and fs/mm bookkeeping")
w = S.lambda_um_to_omega(0.6822)
P(f"  omega(682.2 nm) = {w:.6e} rad/s = {w*1e-12:.4f} rad/ps;  back: {S.omega_to_lambda_um(w):.6f} um")
P(f"  k'(s/m) * 1e12 = fs/mm?  1 s/m = 1e15 fs / 1e3 mm = 1e12 fs/mm -> factor 1e12: {1e12:.0e} (landscape uses gd_ps_per_mm = k'*1e9: 1 s/m = 1e12 ps/1e3 mm = 1e9 ps/mm ok)")
P(f"  GD at A pump (o, 682.2): n_g/c = {ng('y',0.6822)/C_UM_PS*1e3:.6f} ps/mm (CSV GD_p 6.1973)")

# ---------------------------------------------------- 2. biaxial principal planes
P("\n=== 2. Biaxial optics: principal-plane eigenindices vs general Fresnel equation ===")


def fresnel_indices(lam, theta, phi):
    nx, ny, nz = (n(a, lam) for a in "xyz")
    s = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
    a = np.array([1 / nx**2, 1 / ny**2, 1 / nz**2])
    # eigenproblem of eps^-1 projected on plane normal to s
    ref = np.array([0, 0, 1.0]) if abs(s[2]) < 0.9 else np.array([1.0, 0, 0])
    u = np.cross(s, ref); u /= np.linalg.norm(u); v = np.cross(s, u)
    M = np.array([[u @ (a * u), u @ (a * v)], [v @ (a * u), v @ (a * v)]])
    wv, V = np.linalg.eigh(M)
    Dfast = V[0, 1] * u + V[1, 1] * v
    Dslow = V[0, 0] * u + V[1, 0] * v
    return (1 / np.sqrt(wv[1]), Dfast), (1 / np.sqrt(wv[0]), Dslow)   # (n_fast, D), (n_slow, D)


lam = 1.3645
for plane, phi, thetas in (("xz", 0.0, (10, 26.2, 30, 45.968, 70, 90)), ("yz", 90.0, (11.1, 45, 168.9)), ("xy", None, (10, 30, 60))):
    for thd in thetas:
        th = np.radians(thd)
        if plane == "xy":
            nf, ns_ = fresnel_indices(lam, np.pi / 2, th)
        else:
            nf, ns_ = fresnel_indices(lam, th, np.radians(phi))
        n_o, n_e = B.principal_plane_index(lam, th, plane, BIBO["x"], BIBO["y"], BIBO["z"])
        # which code branch is fast?
        P(f"  {plane} angle={thd:7.3f}: code (o,e)=({n_o:.6f},{n_e:.6f}); Fresnel (fast,slow)=({nf[0]:.6f},{ns_[0]:.6f}); "
          f"D_fast={np.round(nf[1],3)} D_slow={np.round(ns_[1],3)}; sets agree={np.allclose(sorted([n_o,n_e]),sorted([nf[0],ns_[0]]),atol=1e-12)}")

P("  -> x-z plane: o-wave (code 'o') is y-polarised, n_o = n_y; e-wave is x-z polarised; e is SLOW (n_e > n_y) between the optic axes.")
Om = S.optic_axis_angle_deg(0.6822, BIBO["x"], BIBO["y"], BIBO["z"])
P(f"  optic-axis half angle Omega(682.2 nm) = {Om:.3f} deg; theta_A = 45.97 lies in (Omega, 180-Omega): {Om < 45.97 < 180-Om}")
# e index equals n_y at theta = Omega
n_o, n_e = B.principal_plane_index(0.6822, np.radians(Om), "xz", BIBO["x"], BIBO["y"], BIBO["z"])
P(f"  n_e(theta=Omega) - n_y = {n_e-n_o:.2e} (should be 0: theta is measured from z, optic axes at +-Omega from z)")

# walk-off check: arctan((1/n) dn/dtheta) vs angle between D and E
P("\nwalk-off: code arctan((1/n)dn/dtheta) vs angle(D,E) from the Fresnel eigenvector, and closed form")
for lam_, thd in ((1.364471, 45.96847), (2.233682, 46.11023), (1.739152, 43.77045), (1.155982, 41.3598), (0.818696, 11.03554)):
    th = np.radians(thd)
    rho_code = B.walkoff_angle(lam_, th, "xz", BIBO["x"], BIBO["y"], BIBO["z"])
    (nf, Df), (ns_, Ds) = fresnel_indices(lam_, th, 0.0)
    a = np.array([1 / n("x", lam_)**2, 1 / n("y", lam_)**2, 1 / n("z", lam_)**2])
    # the e-wave is the x-z polarised one: pick whichever D has zero y component
    D = Ds if abs(Ds[1]) < 1e-9 else Df
    E = a * D; E /= np.linalg.norm(E)
    rho_DE = np.arccos(abs(D @ E))
    nx_, nz_ = n("x", lam_), n("z", lam_)
    ne = 1 / np.sqrt(np.cos(th)**2 / nx_**2 + np.sin(th)**2 / nz_**2)
    rho_cf = np.arctan(ne**2 / 2 * (1 / nx_**2 - 1 / nz_**2) * np.sin(2 * th))
    P(f"  lam={lam_:.4f} theta={thd:.3f}: code {rho_code*1e3:.3f} mrad, angle(D,E) {rho_DE*1e3:.3f}, closed-form {abs(rho_cf)*1e3:.3f}")

# ----------------------------------------------------------- 3. GVM points
P("\n=== 3. Collinear phase matching and the three degenerate oeo GVM points (independent) ===")


def n_xz(pol, lam, th):
    if pol == "o":
        return n("y", lam)
    return 1 / np.sqrt(np.cos(th)**2 / n2("x", lam) + np.sin(th)**2 / n2("z", lam))


def dn_xz(pol, lam, th):
    if pol == "o":
        return dn("y", lam)
    u = np.cos(th)**2 / n2("x", lam) + np.sin(th)**2 / n2("z", lam)
    du = -np.cos(th)**2 * dn2("x", lam) / n2("x", lam)**2 - np.sin(th)**2 * dn2("z", lam) / n2("z", lam)**2
    return -0.5 * u**-1.5 * du


def k_xz(pol, lam, th):
    return 2 * np.pi * n_xz(pol, lam, th) / lam          # rad/um


def kp_xz(pol, lam, th):
    return (n_xz(pol, lam, th) - lam * dn_xz(pol, lam, th)) / C_UM_PS   # ps/um


def dk(th, lp, ls, li, seq):   # seq = idler,signal,pump
    return k_xz(seq[2], lp, th) - k_xz(seq[1], ls, th) - k_xz(seq[0], li, th)


def pm_angle(lp, ls, li, seq):
    ths = np.linspace(1e-3, np.pi / 2 - 1e-3, 600)
    v = np.array([dk(t, lp, ls, li, seq) for t in ths])
    idx = np.where(np.diff(np.sign(v)) != 0)[0]
    return brentq(lambda t: dk(t, lp, ls, li, seq), ths[idx[0]], ths[idx[0] + 1], xtol=1e-14) if len(idx) else np.nan


def deg_mismatch(lp, seq, which):
    ld = 2 * lp
    th = pm_angle(lp, ld, ld, seq)
    if np.isnan(th):
        return np.nan, th
    kp, ks, ki = kp_xz(seq[2], lp, th), kp_xz(seq[1], ld, th), kp_xz(seq[0], ld, th)
    m = {"0": kp - ks, "90": kp - ki, "45": kp - 0.5 * (ks + ki)}[which]
    return m * 1e6, th    # fs/mm


for which, claim in (("0", (682.2, 45.97)), ("45", (869.6, 43.77)), ("90", (1116.8, 46.11))):
    lps = np.linspace(0.5, 1.5, 300)
    vals = np.array([deg_mismatch(l, "oeo", which)[0] for l in lps])
    roots = []
    for j in range(len(lps) - 1):
        if np.isfinite(vals[j]) and np.isfinite(vals[j + 1]) and np.sign(vals[j]) != np.sign(vals[j + 1]):
            roots.append(brentq(lambda l: deg_mismatch(l, "oeo", which)[0], lps[j], lps[j + 1], xtol=1e-12))
    for r in roots:
        m, th = deg_mismatch(r, "oeo", which)
        kp, ks, ki = kp_xz("o", r, th), kp_xz("e", 2 * r, th), kp_xz("o", 2 * r, th)
        tpm = np.degrees(np.arctan2(-(kp - ks), kp - ki))
        P(f"  theta_PM={which:>2}: pump {r*1e3:.2f} nm, theta {np.degrees(th):.3f} deg, D32={(kp-ks)*1e6:.2f} D31={(kp-ki)*1e6:.2f} fs/mm, "
          f"theta_PM(eq.4)={tpm:.2f} deg   [paper {claim}]")
        # compare with package root at same pump
        cfg = ThreeWaveConfig(plane="xz", pol_pump="o", pol_signal="e", pol_idler="o")
        th_code = find_phase_matching_angle(r, 2 * r, 2 * r, cfg, BIBO["x"], BIBO["y"], BIBO["z"])
        P(f"      package theta at that pump: {np.degrees(th_code):.6f} vs mine {np.degrees(th):.6f} (diff {abs(np.degrees(th_code-th)):.1e} deg)")

# theta_PM sign convention check with the triple point
P("\n  type-I degenerate: tan(theta_PM) = -1 identically (k_s'=k_i')?")
for seq in ("ooe", "eeo"):
    for lp in (0.70, 0.8187, 0.95):
        th = pm_angle(lp, 2 * lp, 2 * lp, seq)
        kp, ks, ki = kp_xz(seq[2], lp, th), kp_xz(seq[1], 2 * lp, th), kp_xz(seq[0], 2 * lp, th)
        P(f"    {seq} lp={lp}: theta={np.degrees(th):.3f}, ks'-ki'={(ks-ki)*1e6:.2e} fs/mm, theta_PM={np.degrees(np.arctan2(-(kp-ks),kp-ki)):.2f}")

# landscape vectorised vs scalar at point A
gf = L.gvm_functions(0.682235, 1.364471, 1.364471, np.radians(45.96847), "xz", "oeo", BIBO["x"], BIBO["y"], BIBO["z"])
P(f"\n  landscape.gvm_functions at A: d31 = {gf['g90']*1e12:.3f} fs/mm, d32 = {gf['g0']*1e12:.3f}, theta_pm = {gf['theta_pm']:.3f}")

# ----------------------------------------------------------------- 4. d_eff
P("\n=== 4. d_eff: full point-group-2 Voigt contraction vs LPR2010 Table-1 expressions in the code ===")
d = D_COEFFS_PM_PER_V
DIL = np.zeros((3, 6))
DIL[0, 0], DIL[0, 1], DIL[0, 2], DIL[0, 3] = d["d11"], d["d12"], d["d13"], d["d14"]
DIL[1, 4], DIL[1, 5] = d["d25"], d["d26"]
DIL[2, 4], DIL[2, 5] = d["d35"], d["d36"]
VOIGT = {0: [(0, 0)], 1: [(1, 1)], 2: [(2, 2)], 3: [(1, 2), (2, 1)], 4: [(0, 2), (2, 0)], 5: [(0, 1), (1, 0)]}
DT = np.zeros((3, 3, 3))
for i in range(3):
    for l, jk in VOIGT.items():
        for (j, k) in jk:
            DT[i, j, k] = DIL[i, l]
# sanity: which indices odd in x survive? point group 2 with 2-fold axis x -> d_ijk nonzero iff number of x indices odd
ok = True
for i in range(3):
    for j in range(3):
        for k in range(3):
            nxidx = (i == 0) + (j == 0) + (k == 0)
            if (DT[i, j, k] != 0) and (nxidx % 2 == 0):
                ok = False
P(f"  all nonzero tensor elements have an odd number of x indices (2-fold axis || x): {ok}")


def unit_vectors(plane, angle):
    """propagation s and the o (normal to plane) and e (in-plane, perp. to s) unit vectors."""
    if plane == "xz":
        s = np.array([np.sin(angle), 0, np.cos(angle)]); o = np.array([0, 1.0, 0])
    elif plane == "yz":
        s = np.array([0, np.sin(angle), np.cos(angle)]); o = np.array([1.0, 0, 0])
    else:
        s = np.array([np.cos(angle), np.sin(angle), 0]); o = np.array([0, 0, 1.0])
    e = np.cross(o, s); e /= np.linalg.norm(e)
    return s, o, e


def deff_tensor(plane, angle, seq):   # seq idler,signal,pump
    s, o, e = unit_vectors(plane, angle)
    vec = {"o": o, "e": e}
    ei, es, ep = vec[seq[0]], vec[seq[1]], vec[seq[2]]
    return float(np.einsum("i,ijk,j,k->", ep, DT, es, ei))


cases = [("xz", "oeo"), ("xz", "eoo"), ("xz", "eeo"), ("xz", "ooe"), ("xz", "eoe"), ("xz", "oee"),
         ("yz", "eeo"), ("yz", "oeo"), ("yz", "eoo"), ("xy", "ooe"), ("xy", "eoe"), ("xy", "oee")]
angles = np.radians(np.array([5, 11.1, 30, 45.968, 60, 85, 95, 120, 150, 168.9, 175]))
for plane, seq in cases:
    maxabs, maxsgn = 0, 0
    for a in angles:
        if plane == "xy" and a > np.pi / 2:
            continue
        ref_val = deff_tensor(plane, a, seq)
        code_val = DF.deff_bibo_pm_per_v(plane, seq[2], seq[1], seq[0], a, enforce_range=False)
        maxabs = max(maxabs, abs(abs(ref_val) - abs(code_val)))
        maxsgn = max(maxsgn, abs(ref_val - code_val))
    P(f"  {plane} {seq}: max | |deff_ref| - |deff_code| | = {maxabs:.2e} pm/V;  max |deff_ref - deff_code| (signed) = {maxsgn:.2e}")

# xy eeo / oeo not implemented in the code: tensor value and phase-matchability remark
for seq in ("eeo", "oeo"):
    vals = [deff_tensor("xy", a, seq) for a in np.radians([10, 30, 45, 60, 80])]
    P(f"  xy {seq} (not in LPR2010 Table 1; pump o = z-polarised, largest index, not phase-matchable): tensor deff = {np.round(vals,3)}")

# phase-matchability of xy eeo/oeo: need n_p(o=z) * 1/lp <= ...: check dk sign cannot vanish
P("  check: in x-y, pump o has n_z (largest) -> k_p > k_s + k_i always? sample dk at 0.4->0.8 um, phi=45:")
for phi in (0.2, 0.8, 1.4):
    nz_p = n("z", 0.4); ne_s = 1 / np.sqrt(np.cos(phi)**2 / n2("y", 0.8) + np.sin(phi)**2 / n2("x", 0.8))
    P(f"    phi={np.degrees(phi):.0f}: n_p(z)={nz_p:.4f}  n_e(daughters)={ne_s:.4f}  -> Delta k = 2pi(n_p - n_e)/0.4 = {2*np.pi*(nz_p-ne_s)/0.4:+.3f} rad/um (>0, no PM)")

# octant inequivalence yz
for thd in (11.1, 168.9):
    a = np.radians(thd)
    P(f"  yz eeo theta={thd}: tensor {deff_tensor('yz', a, 'eeo'):+.4f}, code phi=+90 {DF.deff_bibo_pm_per_v('yz','o','e','e',a):+.4f}, code phi=-90 {DF.deff_bibo_pm_per_v('yz','o','e','e',a,phi_sign=-1):+.4f}")
# yz eeo: max |deff| by octant (over all theta, no PM restriction)
th = np.radians(np.linspace(0, 180, 18001))
f = -3.2 * np.cos(th)**2 + 1.76 * np.sin(th)**2 + 1.66 * np.sin(2 * th)
i1 = th < np.pi / 2
P(f"  yz eeo |deff| max over theta<90: {np.abs(f[i1]).max():.3f} at {np.degrees(th[i1][np.argmax(np.abs(f[i1]))]):.1f} deg; "
  f"theta>90: {np.abs(f[~i1]).max():.3f} at {np.degrees(th[~i1][np.argmax(np.abs(f[~i1]))]):.1f} deg (paper: 3.70 / 3.19 along the PM curve)")

# featured points
P("\n  |deff| at the featured points (tensor contraction):")
for tag, plane, seq, lp, thd, claim in (("A", "xz", "oeo", 0.682235, 45.96847, 2.22), ("B", "xz", "oeo", 0.722, 41.3598, 2.40),
                                        ("C", "xz", "oeo", 1.116841, 46.11023, 2.22), ("D", "xz", "oeo", 0.869576, 43.77045, 2.31),
                                        ("E", "xz", "ooe", 0.818696, 11.03554, 3.14), ("F", "xz", "eeo", 0.896757, 34.9082, 1.56),
                                        ("R", "xz", "oeo", 0.810, 39.42211, 2.47)):
    v = deff_tensor(plane, np.radians(thd), seq)
    P(f"    {tag}: {abs(v):.4f} pm/V (paper {claim}); code {abs(DF.deff_bibo_pm_per_v(plane, seq[2], seq[1], seq[0], np.radians(thd), lambda_pump_um=lp)):.4f}")

# Kleinman averages
P("\n  Kleinman averages of the Ghotbi/Hellwig set: d12/d26 (2.93,3.48) -> %.3f; d13/d35 (-1.93,-1.58) -> %.3f; d14,d25,d36 (1.63,1.67,1.67) -> %.3f" %
  ((2.93 + 3.48) / 2, (-1.93 - 1.58) / 2, (1.63 + 1.67 + 1.67) / 3))
P("  'Kleinman violations': (3.48-2.93)/3.205 = %.1f %%, (1.93-1.58)/1.755 = %.1f %%" % (100 * 0.55 / 3.205, 100 * 0.35 / 1.755))

# --------------------------------------------- 5. x-y plane and y-z eeo negative results (independent)
P("\n=== 5. Negative results: x-y plane and y-z eeo (independent re-trace on a coarse grid) ===")


def n_plane(plane, pol, lam, th):
    if plane == "xz":
        return n_xz(pol, lam, th)
    if plane == "yz":
        return n("x", lam) if pol == "o" else 1 / np.sqrt(np.cos(th)**2 / n2("y", lam) + np.sin(th)**2 / n2("z", lam))
    return n("z", lam) if pol == "o" else 1 / np.sqrt(np.cos(th)**2 / n2("y", lam) + np.sin(th)**2 / n2("x", lam))


def kprime_plane(plane, pol, lam, th, h=1e-4):
    # numerical d k / d omega via lambda: k = 2 pi n / lam ; omega = 2 pi c / lam
    def k(l):
        return 2 * np.pi * n_plane(plane, pol, l, th) / l
    dkdl = (k(lam + h) - k(lam - h)) / (2 * h)
    dkdl2 = (k(lam + h / 2) - k(lam - h / 2)) / h
    dkdl = (4 * dkdl2 - dkdl) / 3
    return dkdl / (-2 * np.pi * C_UM_PS / lam**2)


def scan_plane(plane, seq, lp_list, n_ls=60, th_lo=0.001, th_hi=None):
    th_hi = th_hi or (np.pi / 2 if plane == "xy" else np.pi)
    ths = np.linspace(th_lo, th_hi - 0.001, 721)
    found = []
    for lp in lp_list:
        for ls in np.linspace(lp * 1.001, 2 * lp, n_ls):
            li = 1 / (1 / lp - 1 / ls)
            if li > 3.05 or ls < 0.33:
                continue
            v = np.array([2 * np.pi * (n_plane(plane, seq[2], lp, t) / lp - n_plane(plane, seq[1], ls, t) / ls - n_plane(plane, seq[0], li, t) / li) for t in ths])
            idx = np.where(np.diff(np.sign(v)) != 0)[0]
            for i in idx:
                t = brentq(lambda t: 2 * np.pi * (n_plane(plane, seq[2], lp, t) / lp - n_plane(plane, seq[1], ls, t) / ls - n_plane(plane, seq[0], li, t) / li), ths[i], ths[i + 1])
                kp, ks, ki = kprime_plane(plane, seq[2], lp, t), kprime_plane(plane, seq[1], ls, t), kprime_plane(plane, seq[0], li, t)
                tpm = np.degrees(np.arctan2(-(kp - ks), kp - ki))
                found.append((lp, ls, li, np.degrees(t), tpm))
    return found


res = scan_plane("xy", "ooe", np.arange(0.34, 1.52, 0.04))
P(f"  xy ooe: phase-matching roots found: {len(res)} (paper: none)")
res = scan_plane("xy", "eoe", np.arange(0.33, 0.60, 0.005))
if res:
    arr = np.array(res)
    P(f"  xy eoe: roots for pump {arr[:,0].min()*1e3:.0f}-{arr[:,0].max()*1e3:.0f} nm; theta_PM range {arr[:,4].min():.1f} .. {arr[:,4].max():.1f} deg (paper: <396 nm, -21.4..-19.4)")
res = scan_plane("yz", "eeo", np.arange(0.34, 1.52, 0.06), n_ls=40)
arr = np.array(res)
P(f"  yz eeo: {len(res)} roots; theta_PM range {arr[:,4].min():.1f} .. {arr[:,4].max():.1f} deg (paper: [-45, -6.9]); theta range {arr[:,3].min():.1f}-{arr[:,3].max():.1f}")
de = np.array([deff_tensor("yz", np.radians(t), "eeo") for t in arr[:, 3]])
P(f"  yz eeo |deff| along PM roots: max theta<90: {np.abs(de[arr[:,3]<90]).max():.3f}, theta>90: {np.abs(de[arr[:,3]>90]).max():.3f} (paper 3.70 / 3.19)")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "audit_part1_output.txt"), "w") as fh:
    fh.write("\n".join(out))
