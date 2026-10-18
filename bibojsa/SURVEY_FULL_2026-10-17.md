# survey_full — the complete type-II landscape (17 Oct 2026)

Answers item C1 of `physics_check/report.md`: the survey before this date
traced each type-II interaction with the e-polarized daughter fixed as the
*shorter* one (`landscape.py` orders λ_i ≥ λ_s and the sequence string fixes
which daughter is e), so the half of every type-II phase-matching surface on
which the e daughter (x–z, y–z) or the z-polarized o daughter (x–y) is the
*longer* one was never scanned. Here the four partner sequences of LPR2010's
"≡" pairs — x–z `eoo`, x–z `oee`, y–z `eoo`, x–y `oee` — are added and the
whole survey re-run with the unchanged code path (2-nm pump steps, 320 signal
nodes, 1201 angular nodes). The original `survey/` is untouched.

## Code changes (all in `bibojsa/`)

| File | Change |
|---|---|
| `bibojsa/landscape.py` | `region_ok`: `eoo` counts as a "+" interaction (between the optic axes), `oee` as "−". Module docstring: paragraph on the two halves and the label flip. No numerical change for the eight original sequences (verified: 1160/1160 rows of `landscape_curves.csv` identical except the sign of θ_PM = ±90° at four exact-boundary points). |
| `bibojsa/deff.py` | `TYPE_LABELS` entries for the four partner sequences (same type labels). `normalize_oe_sequence` already mapped `eoo→oeo`, `oee→eoe`, so d_eff needs no change. |
| `examples/landscape_survey.py` | `INTERACTIONS` → 12 entries; optional output directory argument (`python examples/landscape_survey.py all survey_full`); degenerate members deduplicated across each "≡" pair (they are the same point). |
| `examples/landscape_headline_full.py` | New. Featured points of the new half (G, H, J, K, M below) through exactly the pipeline of `landscape_headline.py` / `window_tail.py`, plus `family_purity_new.csv` for the new families. |
| `examples/full_landscape_preview.py` | New. Merges the two halves by *physical* condition into `figdata/gvmcurve_{xz,xy}_typeII_{kp_eq_ke,epm,kp_eq_ko}.csv` (column `half`) and draws `fig_full_landscape_preview.png`. |
| `tests/test_landscape_both_halves.py` | New, 4 tests (108 pass): `region_ok` for the partners; the `eoo` 90° curve and the `oeo` 0° curve share point A; the `eoo` crossing at 640 nm (1.130 o / 1.476 e µm, θ = 50.25°) exists while `oeo` has none; the x–y `oee` point at 462 nm (0.618 e / 1.827 o µm, φ = 46.99°, |d_eff| = 1.66). |

A caveat found on the way: `window_tail.py` uses a fixed 1 rad/ps pixel. For
the new points the phase-matching width in the narrow arm is 2π/(L·300 fs/mm)
= 2.1 rad/ps at 10 mm and 0.7 rad/ps at 30 mm, so the 1-rad/ps grid is
under-resolved at L ≥ 20 mm (it gave a spurious 0.984 for H at 30 mm/2 nm).
`window_tail_new.csv` was therefore computed with a 0.25-rad/ps pixel and
`length_extension_new.csv` with 0.4 rad/ps (checked against 0.15 rad/ps:
identical to 4 decimals). The paper's own points are safe (A's narrow-arm
width is 1.2 rad/ps even at 30 mm), but the pixel should be made
scale-aware if `window_tail.py` is rerun for anything with a larger
group-delay scale.

## Coverage (landscape_coverage.csv)

2193 traced GVM points (was 1160): x–z 1024, x–y 63, y–z 1106 (all y–z
`oeo`/`eoo` with d_eff = 0). 1056 points have |d_eff| ≥ 0.5 pm/V, 358 of
them new. The eight original sequences reproduce their previous counts
exactly; the partner sequences add:

| plane / sequence | target 0° | 45° | 90° | note |
|---|---|---|---|---|
| x–y `oee` | 0 | 35 | 28 | **new families**, |d_eff| = d₁₄ sin 2φ up to 1.66 pm/V |
| y–z `eoo` | 98 | 360 | 208 | d_eff ≡ 0 (unchanged conclusion) |
| x–z `eoo` | 52 | 157 | 95 | continuations of the oeo families through A, D, C |
| x–z `oee` | 0 | 0 | 0 | nothing (eoe's partner) |

## The level curves by physical condition (figdata/, preview figure)

Because the (idler, signal) labels swap across each "≡" pair, the *same*
physical condition carries two target labels: in x–z, k_p′ = k_e′ is
`oeo`-0° ∪ `eoo`-90°, EPM is `oeo`-45° ∪ `eoo`-45°, k_p′ = k_o′ is
`oeo`-90° ∪ `eoo`-0°. Each physical level curve is one continuous curve
through its degenerate member; the degenerate members A, D, C are interior
points, not end points. Ranges (merged, from `landscape_curves.csv`):

| plane, condition | pump (nm) | e-wave (µm) | o-wave (µm) | |d_eff| (pm/V) | scale (fs/mm) | ρ_max (mrad) | best unfiltered P |
|---|---|---|---|---|---|---|---|
| x–z, k_p′ = k_e′ (was "0° family, 678–722") | **494–722** | 0.87–1.78 | 0.69–3.04 | 0.16–2.49 (≥ 1.5 for pump ≥ 558) | 89–287 | 9–85 | 0.935–0.954 (new half: 0.943–0.952) |
| x–z, EPM (was "774–882") | **556–882** | 1.03–2.49 | 0.72–3.05 | 0.28–2.49 | 24–175 | 15–83 | 0.78–0.86 (sinc value) |
| x–z, k_p′ = k_o′ (was "828–1117") | **828–1220** | 1.14–3.03 | 2.04–3.04 | 1.84–2.47 | 37–331 | 75–83 | 0.77–0.98 (above 1117 nm the e idler is beyond 2.6 µm for pump > 1170 nm) |
| x–y, k_p′ = k_o′ (new) | **438–492** | 0.57–0.68 | 1.78–1.86 | 0.30–1.66 (≥ 1.0 for 444–486) | 288–306 | 3–19 | **0.948–0.954** |
| x–y, EPM (new) | **484–552** | 0.59–0.71 | 2.50–2.60 | 0.44–1.66 (≥ 1.0 for 490–544) | 176–190 | 5–18 | 0.80–0.83 |

So: the x–y plane is *not* empty; the "three well-separated families" of
x–z oeo are three continuous curves whose pump ranges overlap (the EPM curve
runs under the k_p′ = k_e′ curve from 556 to 722 nm); and the 678-nm lower
end of the "0° family" was the 3.05-µm clip of one branch, not the end of
the curve.

## Candidate featured points (headline_points_new.csv; Table-2 convention)

All on the adaptive window, complex f, Gaussian filter on the narrow arm;
P_∞ and C from 6–24σ square windows (0.25-rad/ps pixel).

| | G | H | J | K | M |
|---|---|---|---|---|---|
| plane / sequence | x–z `eoo` | x–y `oee` | x–y `oee` | x–z `eoo` | x–z `eoo` |
| condition (label θ_PM) | k_p′ = k_e′ (90°) | k_p′ = k_o′ (90°) | EPM (45°) | EPM (45°) | k_p′ = k_o′ (0°) |
| pump (nm) | 532.0 | 464.0 | 532.0 | 800.0 | 1150.0 |
| λ_s / λ_i (µm) | 0.774 (o) / 1.703 (e) | 0.622 (e) / 1.824 (o) | 0.674 (e) / 2.532 (o) | 1.336 (o) / 1.994 (e) | 2.168 (o) / 2.449 (e) |
| cut | θ = 67.74° | φ = 44.92° | φ = 31.48° | θ = 48.23° | θ = 48.08° |
| |d_eff| (pm/V) | 1.21 | 1.66 | 1.48 | 2.13 | 2.14 |
| |Δ₃₂| / |Δ₃₁| (fs/mm) | 254.1 / 0 | 297.4 / 0 | 180.3 / 180.3 | 96.7 / 96.7 | 0 / 203.7 |
| L (mm), pump FWHM (nm) | 10, 4 | 10, 2 | 2, 2 | 2, 8 | 10, 24 |
| window half-width (σ) | 5.5 | 6.2 | 6.3 | 6.3 | 5.2 |
| P (unfiltered, on window) | **0.951** | **0.951** | 0.832 | 0.833 | 0.952 |
| P_∞ / C | 0.945 / 0.032 | 0.945 / 0.041 | 0.801 / 0.20 | 0.801 / 0.21 | — (clipped) |
| herald arm (narrow) | signal, 774 nm (o) | signal, 622 nm (e) | idler 2.53 µm | idler 1.99 µm | idler 2.45 µm |
| P₉₀ (filter FWHM) | 0.984 (2.3 nm) | 0.989 (1.3 nm) | 0.919 (94 nm) | 0.920 (106 nm) | 0.983 (29 nm) |
| P₇₀ (filter FWHM) | 0.993 (0.8 nm) | 0.995 (0.4 nm) | 0.983 (38 nm) | 0.982 (43 nm) | 0.992 (9.9 nm) |
| ρ_max (mrad) | 62.6 (e daughter) | 18.6 (e pump; e daughter 16.7) | 15.8 | 83.0 | 81.0 |
| d_eff/ρ (pm V⁻¹ mrad⁻¹) | 0.019 | **0.089** | 0.094 | 0.026 | 0.026 |
| λ_i below 2.6 µm | yes | yes | yes (2.53) | yes | yes |
| 30 mm | 0.969 (2 nm) | 0.968 (1 nm; the 2-nm grid value 0.947 is off-optimum) | | | |

Comparison with the paper's Table 2: A 0.949 / d_eff/ρ 0.026 / ρ 85; E 0.212 /
0.102 / 31. H matches A's purity with 3.4× A's walk-off figure of merit (87 %
of E's), heralds on a visible photon with a 1.3-nm filter, and its d_eff is
Kleinman-insensitive ((d₁₄+d₂₅)/2 = 1.65 with the measured set). G puts the
"0° family" physics on a 532-nm pump (frequency-doubled Nd/Yb). K is an EPM
member at 800 nm with both daughters below the absorption onset — a cleaner
Ti:sapphire EPM point than R (idler 2.81 µm). M continues C to longer pumps.

Unheralded g⁽²⁾ = 1 + P: H 1.951, G 1.951.

## What this changes in the manuscript (for decision)

1. Abstract / Sec. 1 / Sec. 3 / conclusions: "complete landscape",
   "only the x–z oeo interaction is useful", "nothing useful lies in the
   principal planes other than x–z", "x–y plane yields nothing", "three
   well-separated GVM families (678–1117 nm)", "1160 GVM points, 720 in
   x–z" — all need rewriting. Suggested reading: two useful type-II
   interactions, x–z oeo≡eoo (pump 494–1220 nm) and x–y eoe≡oee (438–552
   nm); the y–z plane and the x–z "−" type-II (eoe≡oee) stay negative.
2. Sec. 2.4 (level-curve method): state that each "≡" pair is two halves
   of one surface and that both are traced; the labels flip across the
   pair, so families are better named by the physical condition
   (pump matched to the e wave / to the o wave / EPM) than by θ_PM.
3. Figs. 1–3: redraw from `figdata/` with the full curves (one continuous
   curve per condition through its degenerate member; the x–y plane needs
   its own panel or the type-I panels make room). Fig. 3 (best purity vs
   pump) gains the 494–680-nm asymmetric branch and the two x–y families.
4. Table 2: add H (and G and/or K; R could be replaced by K). Sec. 4.1
   then has a second high-purity class: the x–y asymmetric family, 0.95
   at 10 mm with ρ < 20 mrad.
5. Secs. 4.3 and 4.5 (brightness): E is no longer "at least 3.4× any
   other" in d_eff/ρ once H (0.089) is featured; the walk-off framing
   (82–85 mrad "the price of BIBO's purity") is specific to the x–z cuts.
   The absolute-rate paragraph could add H (ρ 18.6 mrad → matched waist at
   10 mm is 105 µm instead of 482 µm; 21× higher proxy B than A).
6. Sec. 5.3 (walk-off compensation): the end-for-end flip (180° about the
   in-plane transverse axis) multiplies the lab-frame coupling by
   (−1)^(number of o waves) and a 180° rotation about the beam by (−1)³;
   these are the only two orientations that reverse the walk-off (a
   rotation about the o axis does not). For oeo/ooe the flip keeps the
   sign, as the paper says; for H (`oee`, one o wave) *both* reverse it,
   so a compensated H stack needs a π phase reset at each interface
   (Zondy-type design) — worth one sentence if H is featured.
7. Secs. 2.5 / 5.3 (off-plane scan): the loci of H should be checked off
   the x–y plane as A, D, C were (toward x–z the walk-off grows).
8. Data availability: `survey_full/` replaces `survey/` (or the full run
   is written over it) before the Zenodo release.

## Files

`landscape_curves.csv`, `landscape_degenerate.csv` (identical to `survey/`
for the old sequences), `landscape_coverage.csv`, `featured_new.csv`,
`headline_points_new.csv`, `headline_filter_new.csv`,
`headline_convergence_new.csv`, `window_tail_new.csv`,
`length_extension_new.csv` (G, H), `headline_jsa_{G,H,J,K,M}.csv`,
`family_purity_new.csv`, `figdata/gvmcurve_*.csv`,
`fig_full_landscape_preview.png`, run logs. `landscape_field.csv` (the
contour field) was not regenerated — it is no longer used by any figure.
