# bibojsa — Progress Report

Computational core for "Group-velocity-matched photon-pair generation in the
monoclinic crystal BiB₃O₆: joint spectral engineering and heralded-photon
purity," built per `claude_code_prompt_bibojsa.md`, milestone by milestone,
with each validation stopped and reviewed before proceeding.

## Package layout (implemented so far)

```
bibojsa/
  sellmeier.py      Sellmeier dispersion (3 functional forms) + converged
                     finite-difference k', k'' (dispersion derivatives)
  biaxial.py         Biaxial index ellipsoid: general Fresnel eigenindices,
                     principal-plane closed forms, walk-off angle
  phasematch.py      Delta k(theta), phase-matching angle solver, QPM
                     grating support, QPM-period solver
  gvm.py             Group index/delay, GVM ridge-angle classification,
                     symmetric-GVM pump-wavelength solver, pairwise
                     group-delay helper
  deff.py            Effective nonlinearity, LPR2010 Table 1 verbatim
                     (monoclinic point group 2, Kleinman symmetry,
                     both inequivalent octants, optic-axis range checks)
  jsa.py             Joint spectral amplitude, SVD/Schmidt decomposition,
                     purity, Schmidt number, heralding-arm filter
  materials/
    bibo.py           BIBO Sellmeier (Miyata-Umemura-Kato / LPR2010 Table 2)
                       + Kleinman d_il coefficients (not yet used in code)
    kdp.py            KDP Sellmeier (Zernike 1964)
    ktp.py             KTP Sellmeier: n_x Kato-Takaoka (2002), n_y König
                       & Wong (2004, sourced directly, see Milestone 2
                       update), n_z Fradkin (1999)
tests/
  test_milestone1_kdp.py
  test_milestone2_ktp.py
  test_milestone3_bibo.py
  test_deff_bibo.py
  test_milestone5_bibo_epm_jsa.py
examples/
  locate_epm.py       solve the theta_PM = 45 deg point (x-z ooe, 800 nm)
  epm_jsa_scan.py     JSA/Schmidt/purity scans at that point (+ figures)
```

62/62 tests passing (40 original + 19 deff + 3 milestone-5). Milestone 2's
previously-failing "known gap" test was rewritten into a passing
assertion once König & Wong's n_y was sourced directly (see below).

## Milestone 0 — BIBO Sellmeier vs. LPR2010 Table 3

**Status: done.** Hard-coded the Miyata-Umemura-Kato Sellmeier fit
(LPR2010 Table 2) and validated computed n_x, n_y, n_z and the optic-axis
half-angle Ω against LPR2010's own Table 3 reference values (632.8 nm and
1064 nm), matching to the precision literature gives.

## Milestone 1 — KDP asymmetric GVM (Mosley et al., PRL 2008)

**Status: done, fully validated.** Type-II KDP configuration reproduces
Mosley 2008's asymmetric GVM point (415/830 nm): the phase-matching-ridge
angle θ_PM lands at 0°/90° as expected, confirming the finite-difference
group-index machinery and the phase-matching-angle solver.

## Milestone 2 — KTP symmetric GVM / EPM (König & Wong 2004 / Laudenbach 2017)

**Status: fully validated** (updated — you located König & Wong's paper
itself in `papers/09.`, which resolved the previous gap).

- n_y is now König & Wong's own "modified Fan" fit (their Eq. 3:
  A1=2.09930, A2=0.922683, A3=0.0467695, A4=0.0138408), read directly from
  the paper instead of approximated by a Kato-Takaoka proxy. n_z needed
  no change — the paper confirms they used Fradkin et al.'s n_z
  unmodified, which is what this package already had.
- All five headline numbers now match Laudenbach et al. (2017) to well
  under 1%:

  | quantity | computed | literature | rel. error |
  |---|---|---|---|
  | GD_p (791 nm, o) | 6.0309 ps/mm | 6.027 | 0.06% |
  | GD_s (1582 nm, o) | 5.8833 ps/mm | 5.880 | 0.06% |
  | GD_i (1582 nm, e) | 6.1783 ps/mm | 6.175 | 0.05% |
  | D (dispersion parameter) | 1.0007 | 0.993 | 0.77% |
  | QPM period | 46.146 µm | 46.1 | 0.10% |

- θ_PM = −135.0° (≡ 45.0° mod 180°, up to the ridge-orientation's
  mod-180 ambiguity), correctly landing on the symmetric/EPM locus —
  previously this was misclassified as asymmetric (D≈0) with the
  Kato-Takaoka n_y proxy.
- Previous gap (documented for the record): the earlier n_y proxy
  reproduced n_y's *absolute* value well but not its dispersion *slope*
  near 791 nm, which was enough to shift the classified GVM regime away
  from the literature's symmetric/EPM point. That proxy was adopted only
  because König & Wong's paper was believed paywalled/unsourceable at
  the time; finding it directly removed the need for a proxy entirely.

**Does this affect Milestone 3 or deff? No.** Checked directly in code:
`bibojsa/deff.py` and `tests/test_milestone3_bibo.py` import only from
`bibojsa.materials.bibo` and `bibojsa.sellmeier`/`phasematch`/`gvm` —
nothing from `bibojsa.materials.ktp`. BIBO's dispersion data (Miyata-
Umemura-Kato, LPR2010 Table 2) and KTP's are entirely independent
datasets; KTP only ever served as an external validation crystal for the
GVM/phase-matching machinery itself (milestones 1 and 2), not as an input
to any BIBO computation. So Milestone 3's comparison table and the deff
implementation/validation (Milestone 4 below) stand exactly as before —
this update only strengthens confidence that the shared machinery
(finite-difference group index, phase-matching solver, GVM classification)
is correct, since it now reproduces a second independent literature
benchmark to <1% instead of ~3%.

## Milestone 3 — BIBO sanity vs. LPR2010 Fig. 10b

**Status: done, presented and approved.** Configuration: BIBO, x–z
principal plane, negative type-I (ooe: idler=o, signal=o, pump=e),
λ_pump = 800 nm fixed, phase-matching angle θ solved from Δk=0 at each
scanned signal wavelength.

Computed comparison table:

| λ_s (µm) | λ_i (µm) | θ (deg) | GD_p | GD_s | GD_i | Δ31 | Δ21 | Δ32 |
|---|---|---|---|---|---|---|---|---|
| 1.10 | 2.933 | 10.02 | 6.0281 | 6.0488 | 6.1032 | −75.1 | −54.4 | −20.7 |
| 1.20 | 2.400 | 10.72 | 6.0304 | 6.0379 | 6.0576 | −27.2 | −19.6 | −7.5 |
| 1.30 | 2.080 | 10.92 | 6.0311 | 6.0309 | 6.0381 | −7.0 | −7.2 | +0.2 |
| **1.34** | **1.985** | **10.96** | **6.0312** | **6.0289** | **6.0337** | **−2.5** | **−4.8** | **+2.3** |
| 1.369 | 1.926 | 10.97 | 6.0313 | 6.0277 | 6.0313 | 0.0 | −3.6 | +3.6 |
| 1.45 | 1.785 | 11.00 | 6.0314 | 6.0253 | 6.0268 | +4.6 | −1.5 | +6.0 |
| 1.59 | 1.610 | 11.01 | 6.0314 | 6.0241 | 6.0241 | +7.3 | −0.1 | +7.3 |

(GD in ps/mm; Δ_ij = pairwise group-delay difference in fs/mm.)

Reproduces LPR2010 Fig. 10b qualitatively:

- Cut angle stays in the ~10–11° band described in the paper's text
  ("maximum polar angle slightly above 11°").
- Δ21 shrinks monotonically toward degeneracy.
- Central claim confirmed: Δ31 and Δ32 both cross zero in a band
  bracketing λ_s ≈ 1.34–1.37 µm, matching the paper's stated λ2 ≈ 1340 nm
  point where "signal and idler travel with group velocities very close
  to that of the pump" — all three group delays agree to <0.15% there,
  versus spreads of tens of fs/mm at λ_s = 1.10 µm. This is the project's
  first candidate extended-phase-matching (EPM) point.

4 tests added (`tests/test_milestone3_bibo.py`): phase-matching angle
stays in the expected band, Δ21 shrinks toward degeneracy, all three
group delays converge near 1340 nm, and Δ31/Δ32 both change sign in
[1.25, 1.40] µm.

No deff (nonlinear-coefficient) calculation was implemented yet — not
required for this GVM-only sanity check, and left for explicit go-ahead
per the brief's "ask before adding any physics not specified" instruction.

## Milestone 4 — deff (LPR2010 Table 1), go-ahead given

**Status: done.** Table 1 was extracted directly from the LPR2010 PDF in
`papers/01.` (not from memory) and implemented verbatim in `deff.py`:

| plane | interaction (λ1λ2λ3 = i,s,p) | deff |
|---|---|---|
| x–y (θ=90°) | ooe (I⁻) | −d13 sin φ |
| x–y (θ=90°) | eoe≡oee (II⁻≡III⁻) | d14 sin 2φ |
| y–z (φ=±90°) | eeo (I⁺) | ∓d12 cos²θ ∓ d13 sin²θ + d14 sin 2θ |
| y–z (φ=±90°) | oeo≡eoo (II⁺≡III⁺) | 0 (exactly) |
| x–z (φ=0°, Ω<θ<180°−Ω) | eeo (I⁺) | −d14 sin 2θ |
| x–z (φ=0°, Ω<θ<180°−Ω) | oeo≡eoo (II⁺≡III⁺) | d12 cos θ |
| x–z (φ=0°, θ<Ω or >180°−Ω) | ooe (I⁻) | −d12 cos θ |
| x–z (φ=0°, θ<Ω or >180°−Ω) | eoe≡oee (II⁻≡III⁻) | −d14 sin 2θ |

Octant convention per LPR2010 Sec. 3: first octant plus second octant
(0≤φ≤90°, 90°≤θ≤180°); the paper's φ=180°/270° rows are 2/m duplicates
and are not implemented. θ runs over [0°, 180°] in the y–z and x–z
planes, covering both *inequivalent* monoclinic octants; no
even/odd-in-θ symmetrization anywhere (two-fold axis ∥ dielectric x,
which is why exactly d11, d12, d13, d14, d25, d26, d35, d36 survive).

Validation (19 tests, three independent lines):

- Analytic spot checks of every Table 1 row, equivalences (oee≡eoe,
  eoo≡oeo), invalid-combination rejection, and the x–z optic-axis
  range enforcement (Ω(800 nm) = 25.9° computed and sanity-checked).
- **First-principles cross-check**: |deff| from the Table 1 expressions
  matches a direct Voigt contraction of the point-group-2 d-tensor with
  explicit o/e unit polarization vectors, on a 73-point angle grid per
  row, across all 8 rows, to 1e-10 pm/V. The transcription from the PDF
  is therefore internally consistent with the crystal symmetry, not
  merely copied.
- Literature benchmarks: y–z eeo (the flagship type-I SHG config,
  1064→532 nm) gives |deff| = 3.644 pm/V at θ = 168.9° (second octant)
  vs 2.389 pm/V at the mirror angle 11.1° (first octant) — the ~1.5×
  monoclinic octant inequivalence, matching the known ~3.6–3.7 pm/V
  optimum. The milestone-3 SPDC design point (x–z, ooe, θ ≈ 10.96°,
  800 nm pump) gives |deff| = d12 cos θ = 3.14 pm/V — BIBO's headline
  advantage over BBO (~2 pm/V) for bright pair generation.

API: `deff.deff_bibo_pm_per_v(plane, pol_pump, pol_signal, pol_idler,
angle_rad, lambda_pump_um, phi_sign, enforce_range)` returns the signed
deff in pm/V (only |deff| is physical for our deff²-based brightness
proxy); `deff_bibo_from_config(...)` accepts a `ThreeWaveConfig`
directly. x–z calls check the requested θ against the (wavelength-
dependent) optic-axis regions of Table 1 unless `enforce_range=False`
(useful for plotting deff across the whole plane).

## Milestone 5 — JSA / heralded purity at the 1.34 µm point (x–z ooe, 800 nm)

**Status: done.** First application of `jsa.py` to BIBO, at the exact
θ_PM = 45° point of the milestone-3 candidate region
(`examples/locate_epm.py`, `examples/epm_jsa_scan.py`,
`tests/test_milestone5_bibo_epm_jsa.py`).

**Design point** (this Sellmeier fit): λ_s = 1.3409 µm, λ_i = 1.9831 µm,
θ_cut = 10.9577°, θ_PM = 45.0° (mod 180°), GD = 6.0312/6.0288/6.0336
ps/mm (p/s/i), |deff| = 3.14 pm/V, pump walk-off ≈ 30.5 mrad.

**Central finding (negative, and important for the paper's framing):**
this is a *triple*-GVM point — |Δ31| = |Δ32| ≈ 2.4 fs/mm, ~60× smaller
than the canonical KTP EPM design (147 fs/mm). θ_PM = 45° holds, but the
sinc phase-matching function is orders of magnitude too broadband for
its width to match any realistic pump bandwidth at 0.1–10 mm crystal
lengths; its true extent is set by second-order dispersion (and runs
into the Sellmeier validity window, λ ≤ 3.09 µm, before decaying at
short L — truncation tracked per point via an `edge` metric). The JSA is
a long frequency-anticorrelated stripe, not the near-separable round JSA
of a true EPM source:

| L (mm) | best pump FWHM scanned (nm) | unfiltered purity P | K |
|---|---|---|---|
| 0.1 | 40 | 0.149 | 6.7 |
| 0.5 | 40 | 0.157 | 6.4 |
| 2.0 | 40 | **0.175** | 5.7 |
| 5.0 | 40 | 0.146 | 6.8 |
| 10.0 | 24 | 0.125 | 8.0 |

Best point grid-converged: P = 0.17450/0.17451/0.17453/0.17453 at
N = 100/140/200/260 (per the brief's convergence-check requirement).
Full records in `examples/epm_scan_results.json`; figures (|JSA| maps,
purity scans, Schmidt spectrum) in `figures/` (currently delivered as
session outputs — PNGs can't be written to this folder by the file
tools; drop them in manually or re-run `epm_jsa_scan.py`).

**Interpretation:** the LPR2010 λ₂ ≈ 1340 nm point is a
huge-spectral-acceptance / ultra-broadband-biphoton operating point
(potentially interesting in its own right, e.g. for few-cycle biphotons
or broadband difference-frequency schemes), NOT a high-purity heralded
single-photon design. High-purity BIBO operating points must instead be
sought at asymmetric-GVM loci (θ_PM = 0°/90°) and symmetric loci with
*sizable* |Δ| — precisely what the pending 6-pump exploration survey
scans for. The milestone-5 tests pin all three facts (45° crossing
location, fs/mm triple-GVM scale, low-purity character) so a future
dispersion change that silently alters this conclusion fails CI.

**Engineering note:** `jsa.phase_matching_function` is now fully
vectorized (`delta_k` and the Sellmeier chain are numpy-array-safe; the
old per-point `np.vectorize` loop was ~100× slower). No numerical
change — the full pre-existing suite passes unmodified.

**Milestone-2 note:** the König & Wong n_y update (see Milestone 2
above) does not touch this milestone — BIBO computations import nothing
from `materials/ktp.py`. The full suite (62 tests, including the updated
milestone-2 assertions) passes against the updated package.

## Milestone 6 — Full exploration survey (data sets only, per Alex: no PNGs)

**Status: done.** `examples/bibo_survey.py` (loci + crossings) and
`examples/survey_points.py` (per-point JSA studies); all deliverables
plain CSV for OriginLab — see `survey/README.txt` for the column-by-
column format guide. Full data set delivered as a zip (survey_loci.csv:
1944 phase-matched points, survey_crossings.csv: 54 crossings,
survey_summary.csv: 25 ranked unique operating points, plus 5 CSVs per
point: purity grid, Schmidt spectrum, |JSA| matrix, walk-off/brightness).
`survey_summary.csv` + README are also in `bibojsa/survey/` here.

**Survey scope:** pumps 405/515/775/780/800/810 nm × all LPR2010 Table 1
principal-plane interactions × both octants; per phase-matched point:
θ_PM, signed deff, walk-offs, group delays.

**Headline findings:**

1. **No θ_PM = 0°/45°/90° crossings exist for 405 or 515 nm pumps**, and
   the x–y plane phase-matches nowhere in the scanned ranges. All 54
   crossings live in the x–z plane with near-IR pumps (775–810 nm),
   plus deff=0 y–z oeo crossings (excluded).
2. **Best defensible BIBO source: x–z oeo (II⁺≡III⁺) symmetric-GVM
   point, 810 nm pump** — λs = 1.139 µm / λi = 2.805 µm, θ = 39.4°,
   |deff| = 2.47 pm/V: unfiltered P = 0.81 (L = 5 mm, 8 nm pump FWHM),
   0.87 with a 3× idler filter. The 800 nm twin (λi = 2.87 µm) gives
   0.79/0.86.
3. **x–z eeo (I⁺) type-45 family (775–810 nm)**: P = 0.77–0.79
   unfiltered / 0.84–0.86 filtered at L = 10 mm, |deff| = 1.56 pm/V,
   λs ≈ 1.07–1.17 µm, λi ≈ 2.6–2.8 µm, θ ≈ 35.3° — a robust,
   pump-tunable EPM family (|Δ32| = |Δ31| ≈ 14 fs/mm — small but
   sufficient at 10 mm).
4. The milestone-5 conclusion generalizes: the entire ooe (I⁻) triple-GVM
   family (all crossings at θ ≈ 11°) stays below P ≈ 0.27 despite its
   |deff| = 3.14 pm/V — brightness champion, purity dud.
5. **Truncation caveat (documented in README + per-point `edge`
   metric):** oeo/eoe type-45 points at 775–800 nm have idlers at
   2.87–3.03 µm; their JSAs are clipped by the Sellmeier validity limit
   (3.09 µm) and their apparent P ≈ 0.89–0.91 values are truncation
   artifacts. `survey_summary.csv` therefore carries separate
   `best_lowtrunc_*` columns (edge < 0.15 only) — use those. Extending
   BIBO's dispersion data into the mid-IR (or an IR-validated Sellmeier
   fit) would be needed to firm up the 775–800 nm oeo points; the 810 nm
   oeo point already passes at edge < 0.15.

**Orientation figures (added on request):** `examples/survey_figures.py`
renders 6 PNGs *from the CSV data sets* (they are a viewing layer, not a
data source): fig1 phase-matching loci, fig2 θ_PM along the x–z loci
with crossing markers, fig3 |deff| along the loci, fig4 purity ranking
bar chart, fig5 |JSA| maps of the top-6 points (the Sellmeier-clip
truncation of the oeo points is directly visible at the 3.05 µm idler
edge), fig6 purity vs L and vs bandwidth for the top 3. Delivered as
session outputs (PNGs can't be written into this folder by the file
tools); regenerate any time with the script.

**Numerical hygiene note:** mid-survey, per-JSA-case cost dropped 100×
after discovering Python was executing stale bytecode (`__pycache__`
mtime confusion on the sandbox mount had silently reverted jsa.py to the
pre-vectorization version, and masked two truncated source files). All
data in the deliverable was produced/reproduced after the fix, and the
full test suite (62/62) passes on the exact code that generated it.

## Milestone 7 — Convergence tables for the top-3 operating points

**Status: done.** `examples/survey_convergence.py`; data in
`survey/convergence_all.csv` (per-point files in the survey zip's
`survey/convergence/`). 2-D convergence per the brief: grid size
N ∈ {80…300} × window scale s ∈ {0.7…1.3} (s multiplies the adaptive
window before the Sellmeier-validity clip), at each point's best
low-truncation (L, pump FWHM):

| point | config | L, bw | unfiltered P (s: 0.7→1.3, N=300) | filtered P (3× idler filter) |
|---|---|---|---|---|
| #45 | 810 nm xz/oeo t45 | 5 mm, 8 nm | 0.8194 → 0.8082 | **0.8696** (flat to <5·10⁻⁵) |
| #31 | 800 nm xz/oeo t45 | 10 mm, 4 nm | 0.8040 → 0.7885 | **0.8622** (flat) |
| #42 | 810 nm xz/eeo t45 | 10 mm, 16 nm | 0.8065 → 0.7701 (still falling) | **0.8540** (flat) |

Conclusions for the paper:

- **Grid resolution is a non-issue**: |ΔP| < 8×10⁻⁴ from N = 80 to 300
  at every window, every point.
- **Unfiltered purity is not window-convergent** on these hard-sinc
  JSAs: the slowly decaying (~1/x) sinc tails keep adding Schmidt
  weight, giving a monotone downward drift with window size (worst for
  #42, −0.036 across the s range and not yet plateaued). Unfiltered
  numbers must be quoted with an explicitly defined spectral window, or
  presented as window-dependent upper bounds. For #45/#31 the idler
  window is additionally pinned at the 3.05 µm validity clip at all s
  (edge ≈ 0.06–0.11), folding the known truncation systematic into the
  same caveat.
- **The filtered purity is the robust, quotable number**: with a
  Gaussian idler-arm heralding filter (intensity FWHM = 3× pump's, in
  angular frequency) purity converges in *both* N and s — drift
  < 2×10⁻³ over the entire window range and < 5×10⁻⁵ for s ≥ 1.15:
  **P_filt = 0.870 (#45), 0.862 (#31), 0.854 (#42)** — because the
  filter suppresses both the sinc tails and the region near the
  validity clip. This also cleanly resolves the milestone-6 truncation
  caveat for the two oeo points at the filtered level.

## Milestone 8 — Notebooks (closing the brief's deliverables list)

**Status: done.** Both notebooks built, executed end-to-end (0 errors)
against the current package (62/62 tests passing):

- `01_validation_kdp_ktp` — Milestones 0–2: BIBO Sellmeier vs LPR2010
  Table 3 (diffs < 4×10⁻⁶); KDP/Mosley asymmetric GVM (θ = 67.8°,
  θ_PM = −0.03°, unfiltered P ≈ 0.94 JSA figure); KTP/König & Wong EPM
  (all five numbers vs Laudenbach < 1%, table + history note on the n_y
  proxy episode).
- `02_bibo_survey` — Milestones 3–7: Fig. 10b group-delay convergence
  plot; deff vs angle across the planes (y–z octant inequivalence and
  x–z curves with optic-axis regions); survey summary table, #45 JSA and
  purity-vs-L curves read from the CSV data sets; convergence plots
  (window-dependence of unfiltered vs window-converged filtered purity)
  and the quotable numbers.

Format note: the repo carries jupytext-style percent-format sources
(`notebooks/01_validation_kdp_ktp.py`, `notebooks/02_bibo_survey.py`,
regenerable/executable with any percent-format tool); the executed
`.ipynb` twins (with embedded outputs/figures) are delivered as session
outputs — drop them into `notebooks/` alongside the sources.

## Milestone 9 — Manuscript draft

**Status: draft v0.2 — retargeted to Optics Express (canonical:
`manuscript/manuscript_oe.tex`).** Optica Publishing Group universal
manuscript template structure (class `opticajnl`, style
`opticajnl.bst`; Prism compiles via Overleaf and adds the class files
automatically). Requirements checked against Optica's author pages
(July 2026): single .tex, root-level zip, lowercase figure callouts,
mandatory backmatter (Funding / Acknowledgments / Disclosures / Data
availability — all present), line numbers available via `lineno`.
Abstract tightened to ~180 words for the Optica format. A stub-class
skeleton with the identical equations, table, figure and backmatter
plumbing plus the full .bib compiles with zero errors and all
citations resolved. Step-by-step submission tasks:
`manuscript/SUBMISSION_CHECKLIST_OE.md`. The REVTeX version
(`manuscript.tex`, v0.1) is retained as a reference copy only.

Structure: Intro (GVM purity engineering, BIBO's appeal, monoclinic
complications) → Theory/methods (Sellmeier, monoclinic deff with the
tensor-contraction verification, θ_PM, exact-Δk JSA, edge metric,
KDP/PPKTP validation incl. the n_y-slope cautionary note) → GVM
landscape (three structuring negative results; all useful crossings in
x–z at 775–810 nm) → Operating points (triple-GVM reinterpretation of
the 1.34 µm point; oeo/eeo EPM families; Table of #45/#31/#42/#35 with
Δ, purities, walk-off, L_eff; convergence subsection with the
filtered-purity argument; walk-off/brightness trade-offs) → Discussion
(vs PPKTP; mid-IR herald pros/cons; four limitations incl. IR Sellmeier
validity) → Conclusions; data & code availability.

Numbers in the draft trace to the survey/convergence CSVs; the
equation/table/citation plumbing was compile-tested (a stub-class
skeleton with the same equations, table, figures and full .bib compiles
error-free with all references resolved; REVTeX itself wasn't available
in the sandbox).

TODOs left in the source for Alex: author/affiliation/funding lines,
bib entries marked "TODO verify" (grice1997/law2000 pages,
hellwig2000, heberle2013 title, graffitti2018optica/pra author lists),
repository/DOI for the data availability statement. (The former
"figure remakes in Origin" TODO is superseded by Milestone 10 below.)

## Milestone 10 — Manuscript figures as data-connected pgfplots (no Origin)

**Status: done.** All six survey figures are now LaTeX/pgfplots source
in `manuscript/figures/`, drawn at compile time from CSV slices that
are derived deterministically from the survey data sets — no bitmap or
Origin figures anywhere. Chain:
`bibojsa/survey/*.csv → make_figure_data.py → figures/data/*.csv →
fig<N>.tex → fig<N>.pdf` (see `manuscript/figures/README.md`).

- fig1 pm loci, fig2 θ_PM, fig3 |deff| (6-panel group plots, both
  octant branches, crossing markers), fig4 purity ranking (20
  low-truncation points + PPKTP benchmark line), fig5 |JSA| maps of the
  top-6 points (140×140 meshes, viridis), fig6 purity vs L / vs pump
  FWHM for the top 3 with filtered twins.
- `manuscript_oe.tex` gained a `\iffigsfromdata` toggle:
  false (default) includes the pre-built vector PDFs
  (submission mode, any engine); true redraws every figure from
  `figures/data/*.csv` at manuscript compile time (lualatex, 2 passes).
- All figures now have callouts in the text: fig:loci, fig:thetapm,
  fig:deff, fig:ranking, fig:jsa, fig:purityscans (figs 3, 4, 6 and
  their callouts are new).
- Build: `figures/Makefile` (`make data` + `make`); fig5 requires
  lualatex (the 6×19600-cell meshes exceed pdflatex's memory, ~2.5 min
  compile is normal). Verified: all six PDFs build; the manuscript
  compiles clean in BOTH toggle modes against a local opticajnl stub
  (11 pp, all refs/citations resolved).

## Known limitations (carried forward, documented in code)

- BIBO xyz dielectric frame's small wavelength-dependent rotation
  (≤ ±1°) about the crystallographic b-axis is neglected, per the brief.
- Walk-off corrections to deff are neglected (per the brief; noted in
  deff.py's docstring). Walk-off will enter only via L_eff in the
  exploration phase.
- BBO comparison (used in LPR2010 Fig. 10b/10a for contrast) is out of
  scope / not implemented — only BIBO's own internal GVM trend was
  checked in Milestone 3.

## Pending

- Firm up the 775–800 nm x–z oeo type-45 points (idlers 2.87–3.03 µm):
  needs an IR-extended/validated BIBO Sellmeier fit so the JSA window
  isn't clipped at 3.09 µm (see Milestone 6 truncation caveat).
- (Nothing else from the original brief — all deliverables complete.
  Optional follow-ups: IR-extended BIBO Sellmeier fit for the oeo
  points, and drafting the paper itself from PROGRESS.md + survey data.)
