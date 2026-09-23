"""Absolute pair-rate estimate at point A (BIBO, x-z oeo, 682.2 nm -> 2 x 1.364 um).

Formula: Ling, Lamas-Linares & Kurtsiefer, PRA 77, 043834 (2008), Eq. (36):
collinear degenerate SPDC into matched single Gaussian modes
(W_p = W_s = W_i = W), thin-crystal (loose-focusing) limit, no walk-off:

    R = 4 d^2 P L omega_p^2 / (9 pi n_p n_s n_i eps0 c^2 W^2 Delta_n)

Ling et al. write Delta_n = n_i - n_s (phase indices, dispersion
neglected in d(Delta k)/d omega_s); the spectral integral of sinc^2 is
really 2 pi / (L |k_s' - k_i'|), so Delta_n is the GROUP-index
difference c |k_s' - k_i'| here (0.0508 at A).

Walk-off: the e-polarized daughter beam is displaced by rho z; with
amplitude modes of waist W the z-dependent overlap is exp(-(rho z)^2 / (2 W^2))
(product of two Gaussians of equal waist displaced by rho z), so the
coherent effective length is L_eff = int_0^L exp(-(rho z)^2/(2W^2)) dz and the
rate scales as (L_eff/L)^2 relative to the no-walk-off value.  This is a
transverse-overlap estimate only (the longitudinal phase is unaffected
in the plane-wave limit); it is indicative, like Eq. (8) of the paper.
"""
import numpy as np
from scipy.integrate import quad
import offplane as o

eps0 = 8.8541878128e-12
c = 299792458.0
lp, ls = 0.6822, 1.3645          # um
s = o.direction(45.968, 0.0)
n_p = o.n_mode(lp, s, 'f'); n_s = o.n_mode(ls, s, 's'); n_i = o.n_mode(ls, s, 'f')
ks, ki = o.kprime(ls, s, 's'), o.kprime(ls, s, 'f')      # ps/um
dn_g = abs(ks - ki) * o.C_UM_PS                            # dimensionless
rho = o.walkoff(ls, s, 's')                                 # rad
d = 2.224e-12                                               # m/V
wp = 2 * np.pi * c / (lp * 1e-6)
print(f"n_p={n_p:.5f} n_s={n_s:.5f} n_i={n_i:.5f}  Delta n_g={dn_g:.5f}  (phase-index diff n_s-n_i={n_s-n_i:.5f})  rho={rho*1e3:.2f} mrad")


def rate_thin(L_m, W_m, P_W=1e-3):
    return 4 * d ** 2 * P_W * L_m * wp ** 2 / (9 * np.pi * n_p * n_s * n_i * eps0 * c ** 2 * W_m ** 2 * dn_g)


def walkoff_factor(L_m, W_m):
    Leff = quad(lambda z: np.exp(-(rho * z) ** 2 / (2 * W_m ** 2)), 0, L_m)[0]
    return (Leff / L_m) ** 2


def bandwidth_nm(L_m):
    # sinc^2 FWHM of the narrow (idler) marginal: dk L/2 = +-1.39 -> d omega = 2*1.39*2/(L |ks'-ki'|)
    dw = 2 * 1.3916 * 2 / (L_m * 1e3 * abs(ks - ki) * 1e3)   # rad/ps  (L in mm, k' in fs/mm -> careful units)
    # k' difference in ps/um: dn_g / C_UM_PS ; L in um
    dw = 4 * 1.3916 / (L_m * 1e6 * abs(ks - ki))             # rad/ps
    return dw * (ls ** 2) / (2 * np.pi * o.C_UM_PS) * 1e3   # nm


print("\n geometry                      L(mm)  W(um)  R_thin(/s/mW)  walk-off factor  R_est(/s/mW)  idler FWHM(nm)")
for name, L, W in (("walk-off-matched, 30 mm", 30e-3, 1450e-6),
                   ("walk-off-matched, 10 mm", 10e-3, 482e-6),
                   ("walk-off-matched, 2 mm", 2e-3, 96e-6),
                   ("walk-off-matched, 0.5 mm", 0.5e-3, 24e-6),
                   ("10 mm, 60-um waist", 10e-3, 60e-6),
                   ("1 mm, 60-um waist", 1e-3, 60e-6),
                   ("1 mm, 30-um waist", 1e-3, 30e-6)):
    Rt = rate_thin(L, W); wf = walkoff_factor(L, W)
    zR = np.pi * W ** 2 * n_s / (ls * 1e-6)
    print(f" {name:28s} {L*1e3:5.1f} {W*1e6:6.0f} {Rt:12.0f} {wf:15.2f} {Rt*wf:12.0f} {bandwidth_nm(L):12.2f}   (z_R = {zR*1e3:.0f} mm)")

# sanity: Ling's BBO example, their own inputs (phase-index difference as they used it)
# BBO type-II at 702 nm (351.1 nm pump), theta_pm ~ 49.2 deg: n_o(702)=1.6642, n_e(702)=1.5487
no, ne = 1.6642, 1.5487
th = np.radians(49.2); ne_th = 1 / np.sqrt(np.cos(th) ** 2 / no ** 2 + np.sin(th) ** 2 / ne ** 2)
n_pB = 1.6640  # e-pump at 351 nm ~ n_e(theta) at 351: n_o(351)=1.7062, n_e(351)=1.5772 -> n_e(theta)~1.641; use 1.641
n_pB = 1 / np.sqrt(np.cos(th) ** 2 / 1.7062 ** 2 + np.sin(th) ** 2 / 1.5772 ** 2)
dB = 0.9e-12; wpB = 2 * np.pi * c / 351.1e-9
RB = 4 * dB ** 2 * 1e-3 * 2e-3 * wpB ** 2 / (9 * np.pi * n_pB * no * ne_th * eps0 * c ** 2 * (82e-6) ** 2 * abs(no - ne_th))
print(f"\n Ling 2008 BBO check with their inputs and phase-index difference {abs(no-ne_th):.4f}: R_T = {RB:.0f} /s/mW "
      f"(they report 2 x 0.23 x R_T = 1100 -> R_T ~ 2390; observed ~800)")
