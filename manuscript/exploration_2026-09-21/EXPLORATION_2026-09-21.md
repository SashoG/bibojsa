# Exploration of the two optional additions (21 September 2026)

The 18 Sep review left two additions to Alex's call: an off-principal-plane scan of the GVM loci, and an absolute pair-rate estimate at point A. Both were done here with independent code (no `bibojsa` import), both results went into manuscript v0.7 (Sec. 5.3 (i) and Sec. 4.5), and the page count is still 15. Files in this folder: `offplane.py` (eigenindices, D-eigenvectors, group index, walk-off, full point-group-2 d_eff, JSA purity along an arbitrary direction), `epm_locus.py`, `fixed_pump.py`, `pair_rate.py`, `plot_offplane.py`, the CSVs they wrote, the logs, and `offplane_loci.png`.

## 1. Off-principal-plane scan

### Method

For a direction **s**(θ, φ) in the dielectric frame (θ from z, φ from the x–z plane towards y) the two eigenindices and D-eigenvectors are the eigenpairs of the inverse-permittivity tensor projected on the plane normal to **s** (index-ellipsoid section). Group index k′ = ∂k/∂ω along the fixed direction by Richardson-extrapolated finite differences; walk-off ρ = angle(D, E) with E_i = D_i/n_i²; d_eff = Σ e_p,i d_ijk e_s,j e_i,k with the Kleinman-averaged LPR2010 tensor (2-fold axis along x: d11, d12, d13, d14, d25, d26, d35, d36) and the D-eigenvectors as polarization unit vectors (walk-off correction neglected, as in the paper). Interactions are labelled fast (f) / slow (s) in the (idler, signal, pump) order; in the x–z plane between the optic axes the o wave (n_y) is fast, so the manuscript's oeo is (f, s, f).

Validation (all in `offplane_log.txt`): point A reproduced along (45.968°, 0°) to all digits (682.25 nm, Δ₃₁ = 169.4 fs/mm, d_eff = 2.2242 pm/V, ρ = 85.38 mrad); d_eff reproduces every LPR2010 principal-plane expression (x–z oeo d12 cos θ; x–z ooe −d14 sin 2θ; y–z eeo 2.3889/3.6434 pm/V at 11.1°/168.9°; x–y d13 sin φ and d14 sin 2φ); eigenindices agree with `bibojsa.biaxial.eigenindices_general` at three off-plane directions to 1e-6.

The 2-fold axis makes (θ, φ) ≡ (θ, 180° − φ), and k ≡ −k, so φ ∈ [0°, 90°], θ ∈ [0°, 180°] covers every direction; the θ ↔ 180° − θ mirror is a symmetry only at φ = 0 (the manuscript's "θ ≤ 90° retained in x–z" statement), which is why the octants split off-plane.

### Result 1 — the degenerate points persist as loci at almost fixed pump

At degeneracy the phase-matching condition fixes λ_p for each (θ, φ), so each GVM condition is a curve in (θ, φ). The EPM point is the turning point of the degenerate phase-matching curve (g₄₅ = dΔk/dω_p at degeneracy), which sign-bracketing cannot see; `epm_locus.py` finds it as the extremum of θ_pm(λ_p). The (s, f, s) assignment (the manuscript's eoe) has no degenerate GVM locus anywhere, as in-plane.

| locus | φ = 0° (in plane) | φ = 90° (y–z) | pump drift | Δ₃₁ | purity (10 mm, 8 nm, 6σ) |
|---|---|---|---|---|---|
| A, θ_PM = 0° | θ 45.97°, 682.2 nm, d_eff 2.22, ρ 85.4 | θ 24.2°, 645.5 nm, d_eff 0, ρ 50.9 | −37 nm | 169 → 174 fs/mm | 0.9488 → 0.9511 |
| D, θ_PM = 45° | θ 43.77°, 869.6 nm, d_eff 2.31, ρ 83.4 | θ 19.7°, 843.7 nm, d_eff 0, ρ 42.5 | −26 nm | 72.9 (½·mismatch) | — |
| C, θ_PM = 90° | θ 46.11°, 1116.8 nm, d_eff 2.22, ρ 82.0 | θ 25.5°, 1148.6 nm, d_eff 0, ρ 51.9 | +32 nm | 166 → 189 (corrected 25 Sep 2026; was misprinted as 174) | — |

So the off-plane freedom does not tune the degenerate points: the GVM wavelengths are set by the dispersion of the two branches, which off-plane are mixtures of the same three principal indices. The group-delay scale and the purity barely change along the A locus (purity checked at φ = 10° and 45°, both octants); along the C locus the scale grows by 14 % (166 → 189 fs/mm, `offplane_degenerate_loci.csv` × 1000). Corrected 25 Sep 2026 (REVIEW_2026-09-25.md item 6).

### Result 2 — octant splitting and a weak optimum near φ ≈ 10°

|d_eff| splits immediately: on the A locus 2.45 vs 1.82 pm/V at φ = 5° (θ < 90° vs the mirror octant), 2.47 vs 1.26 at φ = 10°. The θ < 90° octant peaks at φ ≈ 10° (A 2.47, D 2.53, C 2.48 pm/V, +9–11 % over in-plane), passes through zero at φ ≈ 35–40° and recovers to ~1 pm/V; the mirror octant passes through zero at φ ≈ 20° and peaks at 1.8–1.9 pm/V near φ = 40–50°; both vanish in the y–z plane, where the manuscript already has d_eff(oeo) = 0. The walk-off falls monotonically (85 → 51 mrad on A) but not fast enough to compensate:

| locus | in-plane d_eff/ρ (pm/V per mrad) | best off-plane, θ < 90° | best off-plane, mirror octant |
|---|---|---|---|
| A | 0.0261 | 0.0293 at φ = 10° (+12 %) | 0.0285 at φ = 45° (+9 %) |
| D | 0.0277 | 0.0307 at φ = 10° (+11 %) | 0.0296 at φ = 45° (+7 %) |
| C | 0.0271 | 0.0305 at φ = 10° (+12 %) | 0.0287 at φ = 50° (+6 %) |

The principal plane is within about 10 % of the off-plane optimum in the brightness figure of merit. A referee who asks "why only principal planes?" now gets a number.

### Result 3 — the 0° family does not gain pump range off-plane

`fixed_pump.py` traces the nondegenerate k_p′ = k_s′ condition of (f, s, f) at fixed pump against φ (θ ∈ 15–75°). At 722 nm — the in-plane limit of the 0° family, where the two arms meet (point B at θ 41.36°, 1.156/1.924 µm, and its partner at θ 40.56°, 1.095/2.121 µm) — solutions exist only at φ = 0. At 700 nm they survive to φ = 30°, at 710 nm to 20°; at 740–800 nm there are none at any azimuth. The family's pump range contracts off-plane (consistent with A drifting to shorter pump), so a Ti:sapphire centred at 780–800 nm still cannot reach asymmetric GVM in BIBO.

### What went into the manuscript

Sec. 5.3 (i), replacing the speculative "would turn each degenerate member into a two-parameter family" sentence, ~10 lines; one clause in the conclusions ("which propagation off the principal plane improves by at most 10 %"). No figure was added (page budget); `offplane_loci.png` and the CSVs are here if you want a supplementary figure — a pgfplots version from `offplane_degenerate_loci.csv` + `offplane_epm_locus.csv` would take an hour.

Note on units (corrected 14 Oct 2026): until 14 Oct the `*_fs_mm` columns of `offplane_degenerate_loci.csv` and `offplane_epm_locus.csv` held ps/mm (the fixed-pump CSVs were already fs/mm) (`offplane.gvm()` scaled k' by 1e3 instead of 1e6). The CSVs were converted in place (x1000) and `offplane.gvm()` now returns fs/mm; `fixed_pump.py` no longer applies its own x1000. All columns are now fs/mm as labelled.

## 2. Absolute pair rate at point A

Formula: Ling, Lamas-Linares & Kurtsiefer, PRA 77, 043834 (2008), Eq. (36) — collinear degenerate SPDC into matched single Gaussian modes (W_p = W_s = W_i = W), thin-crystal limit, no walk-off:

R = 4 d² P L ω_p² / (9π n_p n_s n_i ε₀ c² W² Δn).

Ling et al. write Δn = n_i − n_s with phase indices (they neglect dispersion in dΔk/dω_s); the sinc² integral is 2π/(L|k_s′ − k_i′|), so Δn is the group-index difference c|k_s′ − k_i′| = 0.0508 at A (the phase-index difference happens to be 0.0507 here, so the distinction is cosmetic at A but not in general). Inputs: d = 2.224 pm/V, n_p = 1.8019, n_s = 1.8272, n_i = 1.7765, ρ = 85.4 mrad. Walk-off: the e-daughter is displaced by ρz; with equal-waist Gaussian amplitude modes the overlap at depth z is exp(−(ρz)²/2W²), so L_eff = ∫₀ᴸ exp(−(ρz)²/2W²) dz and the rate scales as (L_eff/L)² — a transverse-overlap estimate only.

| geometry | L (mm) | W (µm) | R, no walk-off (s⁻¹ mW⁻¹) | overlap factor | R estimate | idler sinc² FWHM |
|---|---|---|---|---|---|---|
| walk-off-matched, 30 mm | 30 | 1450 | 322 | 0.43 | 138 | 1.1 nm |
| walk-off-matched, 10 mm | 10 | 482 | 971 | 0.43 | 415 | 3.3 nm |
| walk-off-matched, 2 mm | 2 | 96 | 4 897 | 0.42 | 2 079 | 16 nm |
| 10 mm, 60-µm waist | 10 | 60 | 62 677 | 0.008 | 486 | 3.3 nm |
| 1 mm, 60-µm waist | 1 | 60 | 6 268 | 0.55 | 3 474 | 32 nm |
| 1 mm, 30-µm waist | 1 | 30 | 25 071 | 0.19 | 4 819 | 32 nm |

Two checks: (i) the walk-off-matched line gives R ∝ 1/L (322 → 971 → 4897 from 30 to 10 to 2 mm), exactly the manuscript's "total pair rate falls as 1/L along the same line"; and the 10-mm/60-µm case lands at the same ~450 s⁻¹ mW⁻¹ as the matched 10-mm case, which is the aperture-length argument of Sec. 4.5 in numbers. (ii) Ling's own BBO example (351 → 702 nm, 2 mm, 82 µm, d = 0.9 pm/V) with approximate BBO indices gives R_T ≈ 4100 s⁻¹ mW⁻¹ from the same code against their ≈ 2400 (they observed ~800 after their Ξ = 0.93 walk-off factor and detection), so the implementation is right to the factor-of-two level that the input indices allow.

Comparison used in the paper: bulk KTP under non-critical (walk-off-free) type-II phase matching, 20 mm, 25.1 kHz mW⁻¹ measured (Kim et al., Opt. Express 33, 22087 (2025)). BIBO at A in the walk-off-matched 10-mm geometry is 25–60 × lower; at 1 mm with a 30–60-µm waist, 5–8 × lower. That is the price of critical phase matching, stated as such in Sec. 4.5.

What went into the manuscript: one sentence in Sec. 4.5 (before the "Two qualifications" sentence) with the 10-mm and 1-mm numbers and the KTP comparison, plus the new reference `ling2008`.

## 3. Backmatter

Email confirmed; Funding and Acknowledgments sections removed (Optica treats both as optional). The only TODO left in the .tex is the Zenodo comment in Data availability.
