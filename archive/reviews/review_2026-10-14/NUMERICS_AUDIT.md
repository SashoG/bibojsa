# Code-level numerics audit of the BIBO GVM manuscript (v0.8) and `bibojsa`

Date: 14 October 2026. Scope: the eight items of the audit brief (Sellmeier and
derivatives; biaxial optics and GVM conventions; d_eff; JSA/purity/filter
conventions; GVM classification; Ling pair rate; off-plane script; derived
arithmetic in the text). Tables 1-2, family ranges and point counts were NOT
re-verified wholesale (done in September); where a CSV number was needed it was
read, not recomputed, unless stated.

Method: two scripts in this folder, written against numpy/scipy only and
compared with the package where relevant:

* `audit_part1_optics.py` / `audit_part1_output.txt` - analytic Sellmeier
  derivatives, general Fresnel eigenproblem with D-eigenvectors, independent
  phase matching and GVM root finding, full Voigt-tensor d_eff contraction,
  independent re-trace of the x-y and y-z negative results.
* `audit_part2_jsa_rate.py` / `audit_part2_output.txt` - pump/filter
  conventions, Schmidt-SVD grid-independence test with an analytic kernel,
  own JSA at A and D, filter on signal vs idler arm at A-D, Ling (2008)
  rate re-derived and recomputed, walk-off factor alternatives, off-plane
  reductions and the 11-12 % figure, all derived arithmetic of the text.

Run both from `BIBO_JSA/bibojsa/` (so that `import bibojsa` resolves).
Nothing in the package, manuscript or CSVs was modified.

---------------------------------------------------------------------------

## Findings ranked by severity

### F1 (medium, physics/presentation): at point C the heralding filter is applied to the BROAD arm

* Where: `bibojsa/examples/landscape_points.py:101`, `landscape_headline.py:88`,
  `window_tail.py:75,79` (filter always on `herald_axis=1`, the idler);
  manuscript `manuscript_oe.tex:425-428` ("Gaussian amplitude filter on the
  idler"), Table 2 rows `P_90`/`P_70` for C (`tex:828-829`), Fig. 6 curve C,
  Sec. 4.2 `tex:877-891`.
* What the code does: for every featured point the Gaussian filter multiplies
  the idler axis. At A and B (theta_PM = 0, k_p' = k_s') the phase-matching
  ridge is horizontal and the idler is the narrow arm, so this is the right
  arm. At C (theta_PM = 90, k_p' = k_i') the ridge is vertical: the SIGNAL is
  the narrow arm (marginal intensity FWHM 2.45 rad/ps) and the idler is the
  broad one (23.65 rad/ps). Filtering the idler at C therefore cuts into the
  pump-limited broad marginal.
* What I computed (package JSA on the paper's adaptive window at C, 10 mm,
  16-nm pump, brentq on the transmitted fraction):

  | C, filter arm | T = 90 % | T = 70 % |
  |---|---|---|
  | idler (as in Table 2) | P = 0.9628, FWHM 143.4 nm | P = 0.9825, 66.4 nm |
  | signal (narrow arm)   | P = 0.9841, FWHM 28.8 nm  | P = 0.9928, 9.9 nm   |

  For comparison, A: idler 0.9876 (10.6 nm) vs signal 0.9683 (70 nm); B:
  idler 0.9847 (26.9 nm) vs signal 0.9638 (67 nm); D is symmetric (0.9184
  either way).
* Impact: Table 2 row C (`P_90`, `P_70`, `Delta lambda_90/70`), Fig. 6 (C
  curve), and the reading of Sec. 4.2 ("the asymmetric points are cheap to
  improve" - C currently appears to need a 144-nm filter while a 29-nm filter
  on the other arm does better). No other number changes. The current numbers
  are not wrong for what they are (an idler filter), but C is the one point
  where "idler" is not the heralding/narrow arm, so the row is not comparable
  with A and B. Suggested fix: filter the signal arm at C (or, generally, the
  arm whose marginal is narrower) and say so in the Table 2 caption and
  `tex:425`; or keep the idler and state explicitly that at C it is the broad
  arm. Also note Sec. 5.3 (ii) `tex:1128` "Point C heralds at 2.234 um" - both
  arms are at 2.234 um; the narrow (filterable) one is the e-polarised signal.

### F2 (low, modelling choice in an "indicative" estimate): walk-off factor in `pair_rate.py`

* Where: `manuscript/exploration_2026-09-21/pair_rate.py:43-45` and docstring
  lines 14-20; manuscript `tex:1038-1041` ("roughly 0.4e3 with a
  Gaussian-overlap estimate", "3-5e3 with walk-off").
* What the code does: `walkoff_factor = (L_eff/L)^2` with
  `L_eff = int_0^L exp(-(rho z)^2/(2W^2)) dz`.
* Two issues, which partly cancel:
  1. `(L_eff/L)^2` is the on-axis (Delta k = 0) amplitude squared of ONE
     spectral component. The spectrally integrated rate is, by Parseval,
     proportional to `int_0^L h(z)^2 dz`, i.e. the factor is
     `(1/L) int h^2 dz`, not `((1/L) int h dz)^2`.
  2. The overlap `exp(-(rho z)^2/(2W^2))` is for two equal Gaussians; with
     the pump and the o-daughter fixed and only the e-daughter displaced the
     three-beam overlap is `exp(-2 (rho z)^2/(3 W^2))`.
* Numbers (my recomputation, rho = 85.38 mrad):

  | geometry | script factor | Parseval, 2-beam h | Parseval, 3-beam h | R_thin (/s/mW) |
  |---|---|---|---|---|
  | 10 mm, 482 um | 0.427 | 0.494 | 0.432 | 971 |
  | 1 mm, 60 um   | 0.554 | 0.595 | 0.528 | 6267 |
  | 1 mm, 30 um   | 0.192 | 0.311 | 0.270 | 25067 |
  | 10 mm, 60 um  | 0.008 | 0.062 | 0.054 | 62668 |

* Impact: the quoted "~0.4e3" at 10 mm/0.48 mm is unchanged (0.43 either
  way); "3-5e3 with walk-off" at 1 mm becomes about 3.3-6.8e3 with the
  three-beam/Parseval version; the "20-1400 times" range is unaffected at its
  ends. The manuscript already labels this indicative. Optional wording:
  "3-7e3".

### F3 (low, misleading self-check, no effect on the paper): Ling BBO "sanity check" in `pair_rate.py:68-77`

The script compares the collinear Eq. (36) with Ling et al.'s BBO
experiment and concludes the implementation is "right to the factor-of-two
level" (EXPLORATION md, Sec. 2). Ling's BBO measurement was NON-collinear
(external emission angle ~3.1 deg, Xi = 0.933 between thin and thick limits),
so the 1.7x gap (4121 vs ~2390) is not a test of the collinear formula. The
formula itself is correct (see V6 below); the check should just not be cited
as validation.

### F4 (info): unit-label slip in `offplane_degenerate_loci.csv` confirmed; manuscript used it correctly

`offplane.py:130` multiplies ps/um by 1e3, giving ps/mm in columns named
`*_fs_mm` (A row: `d31_fs_mm = 0.16945` = 169.4 fs/mm). The manuscript's
169 -> 174 and 166 -> 189 fs/mm (`tex:1120-1121`) are the correct x1000
values (`landscape_degenerate.csv` y-z oeo rows give 173.88 and 189.33
independently). `offplane_epm_locus.csv` is in ps/mm too (`d31_fs_mm =
0.0729` at phi = 0), despite the EXPLORATION note saying it is "correct as
labelled". Nothing in the paper is affected; fix the labels before the Zenodo
release.

### F5 (info): sign conventions of d_eff and theta_PM

* Signed d_eff: the code reproduces LPR2010 Table 1 verbatim. A direct
  contraction with the unit vectors `o` = plane normal, `e = o x s` gives the
  same |d_eff| to 1e-15 pm/V for all 12 label variants (incl. `oee`, `eoo`),
  but the opposite sign for x-z `ooe`, y-z `eeo` (at phi = +90) and x-y `eoe`.
  This is only the e-vector sign convention; the paper uses |d_eff| only.
  Octant assignment in y-z agrees: 2.389 pm/V at theta = 11.1 deg, 3.643 at
  168.9 deg for k = (0, sin theta, cos theta), both in the tensor and in the
  code with `phi_sign=+1`.
* theta_PM: `gvm.theta_pm_deg` docstring says "(-90, 90]" but arctan2 returns
  (-180, 180]; `landscape_degenerate.csv` lists 135.0 for the type-I
  degenerate rows while the text (`tex:648`) says -45 deg. Equivalent modulo
  180 deg; consider normalising the CSV for Data File 1.

### F6 (info, cannot verify from here): crystallographic orientation statement

`tex:321` "crystallographic b axis (parallel to dielectric x)". The eight
surviving coefficients (d11 d12 d13 d14 d25 d26 d35 d36) are exactly those
with an odd number of x indices, so the tensor, the LPR2010 Table-1
expressions, Eq. (2) and the "y-z plane is the plane normal to the two-fold
axis" statement are mutually consistent with a two-fold axis along x. Whether
Hellwig's dielectric frame has b || x was not checked against the primary
source in this audit.

---------------------------------------------------------------------------

## Verified as correct (brief)

V1 Sellmeier (`materials/bibo.py`, `sellmeier.py`). MUK coefficients match
the LPR2010 Table-2 values; LPR2010 Table-3 indices reproduced to <= 3.8e-6
(both wavelengths, all axes); 2*Omega(632.8 nm) = 52.599 deg (cross-checked
with the tan^2 Vz form). The step-halved 4th-order central differences on
k(omega) agree with analytic n_g = n - lambda dn/dlambda to 1e-13 relative
and with analytic k'' = lambda^3/(2 pi c^2) d2n/dlambda2 to <= 1e-7 over
0.4-3.0 um, all axes. Units: omega in rad/s internally, lambda in um at the
API; k' [s/m] x 1e12 = fs/mm and x 1e9 = ps/mm as used; GD_p at A = 6.197329
ps/mm matches the CSV. Linear Jacobian for 8 nm at 682 nm vs exact: 3.4e-5.

V2 Biaxial optics (`biaxial.py`). `principal_plane_index` agrees with the
general Fresnel eigenproblem in all three planes to 1e-12; in x-z the "o"
wave is y-polarised (n_y) and theta is measured from z (n_e(theta=Omega) -
n_y = 2e-16, optic axes at +-Omega from z); between the optic axes the e wave
is slow, so "+" interactions (pump o) apply there - consistent with
Omega(682 nm) = 26.15 deg < theta_A = 45.97 deg. Walk-off
arctan((1/n)dn/dtheta) = angle(D, E) from the eigenvector = closed form
(n^2/2)(1/nx^2 - 1/nz^2) sin 2theta at all featured points (85.380, 81.977,
83.399, 84.434, 30.665 mrad). `eigenindices_general` agrees with
`offplane.eigen` (log).

V3 GVM points (independent root finding, analytic k'): oeo degenerate
682.24 nm / 45.968 deg (D31 = 169.45 fs/mm), 869.58 / 43.770 (|D| = 72.86),
1116.84 / 46.110 (D32 = -166.38); package angles agree to 1e-13 deg. Type-I
degenerate k_s' - k_i' = 0 exactly, theta_PM = -45 (= 135). Eq. (4) sign
convention checked numerically: 0 at A, +45 at D, 90 at C; theta_PM = 0 gives
a narrow idler marginal (3.19 vs 31.9 rad/ps at A). Eq. (6) identity holds at
all 1160 CSV points (0 violations at 1e-3 relative).

V4 d_eff: full point-group-2 contraction reproduces every LPR2010 Table-1
expression in |d_eff| to 1e-15 (12 label variants x 11 angles); y-z oeo/eoo
= 0 exactly; x-y eeo/oeo correctly absent (pump on n_z, Delta k > 0
everywhere). Featured points: 2.2242, 2.4018, 2.2185, 2.3108, 3.1408, 1.5581,
2.4720 pm/V = Table 2. Kleinman averages 3.205/1.755/1.657 -> 3.2/1.76/1.66;
"violations" 17.2 % and 19.9 %. Eq. (2) values 2.389/3.643 at 11.1/168.9 deg.
Along the y-z eeo phase-matching curve |d_eff| peaks at 3.70 (theta > 90) and
3.18-3.20 (theta -> 0), text 3.70/3.19.

V5 Negative results re-traced independently (coarse grid): x-y ooe no root;
x-y eoe roots only for 330-395 nm pump with theta_PM in [-21.4, -19.4]; y-z
eeo theta_PM in [-45.0, -7.7] (text -6.9 on the finer survey grid).

V6 JSA/purity (`jsa.py`): sigma = FWHM/(2 sqrt ln2) is the 1/e half-width of
|alpha|^2 (checked numerically), FWHM_intensity = 2 sqrt(ln 2) sigma, 8 nm at
682.241 nm -> 32.375 rad/ps FWHM, sigma = 19.443 (text 19.4). PMF =
sinc(dk L/2) exp(i dk L/2) with np.sinc(x/pi) = sin x / x. Delta k exact at
fixed cut. Schmidt: unweighted SVD on a uniform grid gives the continuum
purity (analytic correlated-Gaussian test: 0.953939 and 0.800000 on two
unrelated anisotropic grids, exact sqrt(1-c^2)). Filter: Gaussian amplitude
whose intensity transmission has the quoted FWHM (0.5 at +-FWHM/2);
"fraction of pairs" = sum|f filt|^2 / sum|f|^2. nm widths use the local
Jacobian lambda^2/(2 pi c) at each arm's centre (B's 102/281 nm = equal
rad/ps windows). Edge metric = max boundary |f| / max |f|. My own
analytic-index JSA at A gives 0.94878/0.94546/0.94386 at W = 6/12/24 sigma and
0.83290 at D (CSV 0.94877/0.94546/0.94386, 0.83282).

V7 Pair rate: Ling et al. Eq. (36) read from arXiv:0801.2220 is
R = 4 d^2 P l omega_p^2 / (9 n_s n_i n_p eps0 pi W_p^2 (n_i - n_s) c^2);
`pair_rate.py:40` implements it verbatim with the group-index difference
in place of (n_i - n_s). Re-derived from H_I = 4 eps0 d E_p E_s E_i
(consistent with P(omega_3) = 2 eps0 d A1 A2), equal-waist Gaussian overlap
2/3, P = n_p eps0 c E0^2 pi W^2/4 and the exact sinc^2 integral
2 pi/(L |k_s' - k_i'|): same coefficient. The group-index replacement is the
correct one (the sinc^2 integral involves d(Delta k)/d omega_s, i.e. k'),
Delta n_g = 0.05080 (phase-index diff 0.05068). SI units consistent. Values:
971 /s/mW (10 mm, 482 um; 979 at exactly 480 um), 6267 (1 mm, 60 um), 25067
(1 mm, 30 um); z_R = 977 mm >> L. KTP comparison 25.1 kHz/0.21^2 = 0.57 MHz
("~0.6"); ratios 24-1446 ("20-1400").

V8 Off-plane (`offplane.py`): at phi = 0 the eigenproblem returns n_fast =
n_y = 1.776517 and n_slow = n_e(theta) = 1.827201 at 1.3645 um; tensor d_eff
2.2242 = d12 cos theta; FOM gains from the CSVs: A +12.4 % (phi = 10), D
+11.0 % (phi = 10, EPM file), C +12.7 % (phi = 10); mirror octant +9.3/+6.0 %;
|d_eff| 2.450 vs 1.815 at phi = 5 (text 2.45/1.82); pump drift 682 -> 645.5,
870 -> 843.7, 1117 -> 1148.6 nm.

V9 Derived arithmetic in the text: x-z birefringence 8.94 % ("close to 9");
L_a = 1.245 mm at 60 um; matched waists 24/96/482/1445 um; L^-0.6 gives
6.3/5.3/4.1 nm at 15/20/30 mm vs grid optima 6/6/4 (fitted exponent -0.59),
3-nm optimum at 51 mm; d_eff/rho row of Table 2 all correct to the digit
given; E/others 3.44-3.93, squared 11.9-15.5 ("12-15"); C ratio 3.8-5.5
("four to six"); Jin d_eff differences 12.6/7.4/3.6 % ("13, 7, 4"); Sellmeier
spread 0.98-1.32 % and up to 0.193 deg (Hellwig-Umemura at 45 deg); K =
1.0871; Fig. 5 panel ratios 4.4 and 7.3; e^-1 at 1 cm^-1 x 10 mm; D
antidiagonal sinc zero 21.56 rad/ps, pump FWHM 29.89; A idler sinc zero 3.71,
pump FWHM 32.38; 10.7 nm = 10.83 rad/ps; grid convergence 120 -> 340 nodes
changes P by <= 3.5e-4 (R) and <= 1.3e-4 elsewhere ("< 5e-4"); filtered
purities at A/B/C move <= 4.3e-4 between 6 and 16 sigma; length_extension
numbers 0.967/0.963 (A, 30 mm, 4 nm), 0.962 (3 nm), 0.965 (B), 0.960 (C)
match the text; P_inf + C/W at W = 6 reproduces 0.9485 (A) and 0.833 (D).

## Not verified / out of scope

* Literature values themselves (LPR2010 Table 2 transcription beyond the
  Table-3 reproduction; Hellwig/Umemura alternative fits; b || x orientation;
  Ghotbi/Hellwig d coefficients; Jin 2020 numbers; kim2025 ratio).
* `points2/` per-point grids (not in this copy): the A-point length series
  0.795 (0.5 mm) / 0.899 (2 mm) in `tex:1022-1024` could not be re-read from
  a CSV.
* Family purity optimisation (Fig. 4) beyond reading `family_purity.csv`
  ranges (0.935-0.954 for oeo 0 deg; 0.770-0.939 for 90 deg; 0.791-0.855 for
  45 deg - consistent with the text).
