# Figure review — manuscript_oe.tex v0.5 (18 Sep 2026)

Scope: the six pgfplots figures in `manuscript/figures/`, their plot-ready
CSVs in `figures/data/`, the survey CSVs they derive from, and the numbers
quoted in the captions and the surrounding text. Everything below was
checked programmatically; the visual checks were made on the committed
`fig<N>.pdf` rasterised at 220–600 dpi and on the compiled `manuscript_oe.pdf`.

## Verdict

The data are right; one figure-drawing bug is not. Every one of the 45
regenerable CSVs in `figures/data/` is byte-identical to a fresh run of
`make_figure_data.py` on the survey, every hard-coded marker in the figure
sources matches the survey to the last digit, and every number the
captions and Secs. 3–4 quote from the figures reproduces from the CSVs.
But `make_figure_data.py` orders the level-curve points wrongly, and as a
result **Fig. 1 and Fig. 3 currently draw only 9 of the 43 points of the
oeo 0° family** — the family that contains the paper's best point A — and
10 of 62 points of the 45° family are missing too. Fix in Sec. 1; the
corrected script and the four rebuilt PDFs are next to this file.

## 1. Blocking: figs 1 and 3 drop and scramble the type-II families

`split_branches()` sorts each family by `(pump_nm, lambda_s_um)` and starts
a new branch wherever `lambda_s` jumps by more than 0.15 µm, then discards
every branch shorter than two points. Every x–z family folds back in pump
(the paper says so, Sec. 3.3), so `landscape_curves.csv` holds two rows per
pump wavelength on the folded part; sorting by pump interleaves the two
arms, `lambda_s` jumps at almost every step, nearly every "branch" is a
single point, and those points are thrown away. What survives is joined in
the wrong order.

Points lost per family (survey rows with |d_eff| ≥ 0.5 pm/V vs rows in
`gvmcurve_*.csv` / `deff_*.csv`):

| family | survey | plotted | lost | what is lost |
|---|---|---|---|---|
| oeo 0° | 43 | 9 | **34** | the whole interior, 686–716 nm on both arms; d_eff 2.23–2.49, Δ 94–168 fs/mm |
| oeo 45° | 62 | 52 | 10 | the upper arm past the fold, 874–882 nm |
| eeo 0° | 42 | 39 | 3 | 898–900 nm, the approach to the triple point F |
| all others | — | — | 0 | unaffected (no fold within the 2-nm step) |

Visible consequences in the committed PDFs (see `fig1_fig3_before_after.png`):

- Fig. 1, left panel: the θ_PM = 0° family is three disconnected stubs
  (one near 3.0 µm, one 1.16–2.2 µm drawn as a zig-zag between arms, one
  at 0.87 µm) instead of the closed hook through A and B that the caption
  describes ("each family carries two branches, which meet at its
  frequency-degenerate member").
- Fig. 3, both panels: the 0° family that the text describes as
  "89–169 fs/mm" and "2.22–2.49 pm/V" is two scribbles at 0.68 and 0.72 µm;
  the 45° family has a stray zig-zag near 0.88 µm at ~60 fs/mm.
- Fig. 4: the same interleaving makes the line through the oeo 0° and 45°
  points zig-zag between arms. Hidden by the markers at print size, but the
  purity CSVs have the same ordering problem.

Fix (in `make_figure_data_fixed.py`, diff in `make_figure_data.diff`):
order each family by `lambda_s_um` instead of pump. Along every traced x–z
family the signal wavelength increases monotonically — verified: ordered by
λ_s, consecutive points of every family are exactly one 2-nm pump step
apart, with no duplicates — so this recovers the curve with no
heuristics; break only on a real pump gap (> 10 nm), and never drop
points. The same ordering is applied to `purity_*.csv`. After the fix all
698 x–z points are plotted; the 15 CSVs that change are `gvmcurve_{eeo_0,
oeo_0, oeo_45}`, `deff_{eeo_0, oeo_0, oeo_45}`, `purity_{oeo_0, oeo_45}`
and the seven `filter_*` files (Sec. 2). `fig1`, `fig3`, `fig4`, `fig6`
rebuilt from the corrected data compile clean with pdflatex ×2; the
rebuilt PDFs are in this folder. Nothing in the figure `.tex` sources needs
to change.

Text check after the fix: the 0° family now visibly spans 678–722 nm with
|d_eff| 2.23–2.49 pm/V and Δ 89–169 fs/mm, as Secs. 3.3 and 4.1 state.

To apply: copy `make_figure_data_fixed.py` over `figures/make_figure_data.py`,
run `make data` then `make` (fig5 unchanged, no lualatex run needed), or
just copy the four PDFs from this folder into `figures/`. Update
`figures/README.md` accordingly and re-run `build.bat`.

## 2. Minor data/figure inconsistencies

**Fig. 6 does not actually reach transmission = 1.** `make_figure_data.py`
drops the rows with an empty `ratio` column, which are exactly the
unfiltered `transmission = 1.0` rows of `headline_filter.csv`. Each curve
therefore ends at T = 0.9998 (A), 0.9985 (D), 0.9946 (E) … while the
caption says "the right-hand end of each curve (transmission = 1) is the
unfiltered purity of Table 2". Invisible at print size (the last plotted
purity agrees with the unfiltered value to 3 decimals), but the statement
is not literally true of the data drawn. The fixed script keeps the row;
pgfplots accepts the empty cells in the unused columns, verified.

**Fig. 5 axis limits leave white frames around two panels.** Panel D is
drawn to ±120 rad/ps but the mesh only reaches ±114.9 (data ±113.9 plus
half a cell), so there is a ~5-rad/ps white border on all four sides.
Panel R is ±100 but the mesh spans x ∈ [−92.2, 93.9], y ∈ [−54.6, 92.2], so
besides the intended blank band at the bottom (the 3.05-µm Sellmeier limit
at δω_i = −53.9 rad/ps — confirmed, it is exactly where the data stop)
there are thin white margins on the top, left and right that the caption
does not account for. Suggested limits: D ±115, R ±92 on both axes (keeps
the equal-scale convention; the bottom band then runs from −92 to −54.6),
E ±509 instead of ±510 (a 0.8-rad/ps sliver, cosmetic). A, B, C are fine
(their windows crop the mesh, and less than 0.35 % of |f|² lies outside).

**Fig. 5 sidelobes in A, B, C are under-sampled.** The sinc side-lobe
period in the narrow direction is 2π/(Δ·L): 3.71 rad/ps at A, 4.58 at B,
3.78 at C. The 120×120 mesh that `JSA_COLS = 120` produces samples these at
2.1, 1.9 and 2.8 points per period — at or below Nyquist — which is why
the side-lobes in the committed PDF look like slightly tilted, beating
stripes rather than straight bands (visible at 600 dpi, and faintly at
print size). The 200×200 survey grid itself only manages 3.5, 3.2 and 4.6
points per period. Two cheap improvements: (i) crop each `headline_jsa`
grid to the panel's axis window *before* subsampling — for A, B, C that
leaves 132, 111 and 148 source points across the window, so no subsampling
is needed at all and the cell count stays under the current 14 400; (ii)
if the survey is ever re-run, compute `headline_jsa_{A,B,C}` on a finer
grid (400×400 over the same range). D, E, R are well resolved (side-lobe
periods 41–43 rad/ps).

**Dead data files.** Fig. 1 draws its boundary lines with `\addplot
{2.60}` / `{3.05}` and never reads `bound_degeneracy.csv`,
`bound_transparency.csv`, `bound_sellmeier.csv`; Fig. 2 hard-codes the
crossing coordinates and never reads `tpmdeg_cross.csv`; the `inwin` column
of `gvmcurve_*.csv` is unused. Harmless, but `figures/README.md` lists
`data/bound_*.csv` as Fig. 1 inputs. Either drop the four files from the
script or make the figures read them (Fig. 2 reading `tpmdeg_cross.csv`
would remove one place where numbers are duplicated by hand).

## 3. Presentation

- **Fig. 6 colour coding is inconsistent with Figs. 1, 3, 4.** There the
  GVM condition is blue/orange/green for 0°/45°/90°. In Fig. 6, A and B
  (0°) are blue and C (90°) green and D (45°) orange — consistent — but R
  (45°) is drawn in `cwarn` red, which Figs. 2 and 3 use for the *eeo*
  interaction, and E (*ooe*) is drawn in `cooe` green, nearly the same hue
  as C. Suggest R in `cg45` with a different dash and E in a hue not used
  for any GVM condition (`ceoe` purple, or black).
- **Fig. 6 vertical range.** The five type-II curves occupy P ∈ [0.81, 1]
  while E alone forces the axis down to 0.1; the region the text discusses
  (P₉₀, P₇₀ of A–D and R) is squeezed into the top fifth. An inset, or
  plotting E on a broken axis / in a small side panel, would let the main
  panel run 0.78–1.0.
- **Fig. 3, right panel:** the legend (south-east) sits on top of the
  eeo 0° curve's fall to the triple point at 0.90 µm and the tail of the
  ooe 90° curve. Free space is at 1.0–1.15 µm below 20 fs/mm; anchoring
  the legend at `(0.99,0.03)` with shorter legend samples, or moving it to
  the left panel (which has room at 0.95–1.15 µm, 0–1 pm/V), clears it.
- **Fig. 4:** the lines joining points are least useful in the two type-I
  panels, where the optimiser's scatter (ooe 0°: 0.03–0.27 over 0.74–0.82 µm)
  is real and the joins read as structure. `only marks` for those panels,
  or a thinner line, would be more honest; lines are fine in the oeo panel
  once the ordering fix is in.
- **Fig. 1, "A" label** now touches the 0° curve once the family is drawn in
  full; shift it left/down by ~0.01 µm, or anchor it south-east of the
  marker.
- **Fig. 5 caption** says F is omitted; the Fig. 6 caption and legend
  omit F too (the script does write `filter_F.csv`) without saying so.
  One clause in the Fig. 6 caption, or add F as a second thin green dashed
  line, whichever you prefer.
- **Effective font sizes after `\includegraphics` scaling** range from
  ~7.6 pt (Fig. 2, 0.76 linewidth) to ~10 pt (Figs. 3–4, footnotesize
  sources scaled ×1.24–1.26). Not wrong, but Figs. 3 and 4 look visibly
  larger-lettered than the rest on pages 7–8; setting their `bibopanel`
  fonts back to scriptsize, or including them at 0.8 linewidth, evens it out.

## 4. Numbers verified (all reproduce)

Hard-coded coordinates in the `.tex` sources vs `landscape_degenerate.csv`
/ `headline_points.csv`: A (682.241 nm, 1.364482 µm), B (722, 1.155982 /
1.923690), C (1116.841, 2.233682), D (869.583, 1.739166), E (818.696,
1.637392), F (896.760, 1.793520) — all exact. Fig. 2 crossings from
`tpmdeg_oeo.csv` interpolate to 682.2 / 869.6 / 1116.9 nm; the type-I
curves are −45.0000° at every point; the oeo curve is monotonic, −9.6° at
0.60 µm to 99.6° at 1.30 µm ("from −10° to beyond 90°" ✓). Fig. 5 panel
titles (pump, L, FWHM) match Table 2; E's window is 4.25× D's and 7.3× A's
("four times … seven times" ✓); every mesh is a full 120×120 with max
amplitude ≤ 1 and no NaNs. Fig. 6 curves interpolated at T = 0.9 / 0.7
reproduce every P₉₀, P₇₀ and filter FWHM in Table 2 to the quoted digits;
all seven curves are monotone in T. Fig. 4: the 0.8315 guide vs "P = 0.83"
✓; no plotted point has edge ≥ 0.15 ✓; 45° family median 0.8308 ✓;
ooe purities 0.028–0.339 ("0.03–0.34" ✓); eeo maximum 0.8897 at 836 nm
with a 2.919-µm idler ✓; 0° family 0.9347–0.9536 ("0.935–0.954" ✓).
Type-I bound: ooe 0.43–12.9 fs/mm ("0.4–13" ✓), eeo 0.49–53.6 ("0.5–54" ✓);
`gvm_scale_fs_mm` equals max(|d31|,|d32|) exactly in every row.

Three small text nits, none affecting a figure:

- Sec. 3.3: "38–166 fs/mm for the 90°" — the family minimum is 37.5 fs/mm
  (828 nm); "37–166".
- Sec. 4.1: the 90° family "passing 0.90 near 1.01 µm" — it is 0.8955 at
  1014 nm and 0.9057 at 1030 nm, so the crossing is nearer 1.02 µm.
- Sec. 4.2: "0.83 ± 0.02 at every pump wavelength along its 774–882 nm
  extent [Fig. 4]" — Fig. 4 shows the 45° family only from 794 nm (the
  774–792 members fail the edge < 0.15 cut), and the 802-nm member is
  0.791, outside ±0.02. "0.83 ± 0.02 over most of its extent" or "0.79–0.86,
  median 0.83" is what the figure shows.

## 5. Files in this folder

| file | what |
|---|---|
| `make_figure_data_fixed.py` | corrected script (ordering by λ_s; keeps the T = 1 filter row) |
| `make_figure_data.diff` | unified diff against the committed script |
| `fig1_gvm_landscape.pdf`, `fig3_deff.pdf`, `fig4_purity_families.pdf`, `fig6_filter.pdf` | rebuilt from the corrected data with the unchanged `.tex` sources |
| `fig1_fig3_before_after.png` | side-by-side of the committed and rebuilt Figs. 1 and 3 |

Not touched: `fig2`, `fig5` (their data are unaffected), the figure
sources, `figures/data/` (the corrected CSVs regenerate in seconds with
`make data` once the script is swapped).

---

## Status — applied 18 Sep 2026 (manuscript header now v0.5.1)

Everything in Secs. 1 and 2 and the three text nits of Sec. 4 have been
applied to `manuscript/` and the manuscript rebuilt (15 pages, no errors):

- `figures/make_figure_data.py` replaced (ordering by λ_s, T = 1 row kept,
  JSA grids cropped to the panel windows before subsampling, `JSA_COLS = 150`,
  per-panel mesh sizes written to `data/jsa_mesh.tex`).
- `figures/data/`: 15 CSVs regenerated + 7 `jsa_*.csv` + new `jsa_mesh.tex`.
- `fig1`, `fig3`, `fig5`, `fig6` sources edited (A label, legend, panel limits
  D ±115 / R ±92 / E ±509, mesh macros, colours); all five affected PDFs
  rebuilt; `fig2` untouched.
- `manuscript_oe.tex`: three text corrections, Fig. 6 caption clause on F,
  header comment; `manuscript_oe.pdf` rebuilt.
- `figures/README.md` updated.

Not applied (Sec. 3, judgement calls): Fig. 6 inset, Fig. 4 marks-only for
the type-I panels, evening out effective font sizes, and pruning the unused
`bound_*.csv` / `tpmdeg_cross.csv` outputs.
