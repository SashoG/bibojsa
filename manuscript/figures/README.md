# Manuscript figures — pgfplots, drawn from the survey data

Every figure is LaTeX/pgfplots source, not a picture. The chain is:

```
bibojsa/survey/*.csv  --make_figure_data.py-->  data/*.csv  --fig<N>.tex-->  fig<N>.pdf
   (the survey)            (plot-ready slices)     (pgfplots)     (vector PDF)
```

| figure | source file | data files | manuscript label |
|---|---|---|---|
| 1 | `fig1_gvm_landscape.tex` | `data/gvmcurve_xz_{ke,epm,ko}.csv`, `data/gvmcurve_xy_{ko,epm}.csv`, `data/gvmcurve_{ooe,eeo}_{0,45,90}.csv` (the 2.60/3.05 µm idler limits are drawn as constants) | `fig:landscape` |
| 2 | `fig3_deff.tex` | `data/deff_*.csv` (same families) | `fig:deff` |
| 3 | `fig4_purity_families.tex` | `data/purity_*.csv` (same families) | `fig:familypurity` |
| 4 | `fig5_jsa_maps.tex` | `data/jsa_{A,H,G,D,K,E}.csv`, `data/jsa_mesh.tex` | `fig:jsa` |
| 5 | `fig6_filter.tex` | `data/filter_{A,B,C,D,G,H,K,E}.csv` | `fig:filter` |

File names keep their v0.9 numbers (`fig3_` ... `fig6_`); the manuscript
numbers are one lower since Fig. 2 (`fig2_tpm_degenerate.tex`, the
Sellmeier comparison of the degenerate points) was retired in v0.10; its
source and data files were removed on 18 Oct 2026 (the built PDF is kept
for reference).

Since v0.11 (17 Oct 2026) the type-II data files are named by the plane
and the physical GVM condition (`xz_ke` = pump matched to the e wave,
`xz_epm` = extended phase matching, `xz_ko` = pump matched to the o wave;
`xy_ko`, `xy_epm` for the x–y plane) and each carries both halves of the
surface (`seq` column: `oeo`/`eoo`, `oee` in x–y). The type-I files keep
the θ_PM naming (`ooe_0`, …).

## Building

```
make                         # builds the five PDFs from data/
python make_figure_data.py   # only if the survey changed: re-derives data/
```

`data/` (47 CSVs plus `jsa_mesh.tex`; only files a figure reads) is kept in the folder. It is derived output —
`make_figure_data.py` regenerates it from `../../bibojsa/survey` in a few
seconds — but it is the exact input the five figures were drawn from, so
it is version-controlled with them rather than reconstructed on demand.
Re-derive it only when the survey data change.

Where the numbers come from, in order:

| layer | location | what it is |
|---|---|---|
| survey output | `../../bibojsa/survey/` | the computation itself, complete: the 2193 traced GVM points of all twelve principal-plane interactions (both halves of every type-II surface), the degenerate members, the headline points A–M with their JSA grids, filter scans and window-tail fits; `survey/README.txt` documents every file |
| plot-ready slices | `data/*.csv` | this folder: the columns each `\addplot` actually reads |

fig1, fig3, fig4 and fig6 compile with `pdflatex` (run twice — the shared
legends travel through the `.aux` file). **fig5 needs `lualatex`**: its
six |JSA| meshes (up to 150×150 cells) exceed pdflatex's static memory
and take a few minutes. It also needs the `luatex85` compatibility
package (TeX Live `texlive-luatex`); if it is missing, `standalone`
cannot crop and the output is a 2-page letter-size PDF with the word
"output" on page 1. On Overleaf, set the compiler to LuaLaTeX if you
compile the figures there.

## Two ways to use the figures in the manuscript

`manuscript_oe.tex` has a toggle near the top:

- `\figsfromdatafalse` (default) — the manuscript does
  `\includegraphics` of the pre-built `fig<N>.pdf`. Fast, engine-
  independent; this is the mode for submission.
- `\figsfromdatatrue` — the manuscript `\input`s the figure sources and
  redraws everything from `data/*.csv` at compile time (via the
  `standalone` package). Compile the manuscript twice with **lualatex**
  in this mode. Use it while iterating on figure style so the
  manuscript always reflects the current data.

## Editing

- Global style (colors, line styles, panel sizes, markers):
  `bibofig-preamble.tex`. The GVM conditions are colour-coded there:
  `cg0` / `cg45` / `cg90` for the three GVM conditions (since v0.11 named
  physically: pump matched to the e wave, EPM, pump matched to the o
  wave; on the type-I families they keep their θ_PM = 0°, 45°, 90°
  meaning); the interactions use a separate palette (`coeo` / `cooe` /
  `ceeo`, purple / teal / brown, and `cxy` for the x–y plane) so that one
  colour never means both an interaction and a GVM condition. fig4 draws members with `inwin = 0` (idler beyond the
  2.6-µm absorption onset) as open markers via `scatter/classes`; fig6
  plots the impurity 1 − P on a log axis.
- Per-figure layout (panels, ranges, legends, annotations): the
  `fig<N>.tex` files — plain pgfplots, every `\addplot` names the CSV
  it draws.
- JSA panels (fig5) are plotted in angular-frequency **detuning** with
  equal scales on both axes, so the ridge orientation on the page is
  comparable to θ_PM. Do not switch them back to wavelength axes: a 45°
  frequency-space ridge appears at ~9° in wavelength space.
- The per-panel axis limits in `fig5_jsa_maps.tex` are hand-set to keep
  the two axes equal, and `make_figure_data.py` crops each survey grid to
  the same window (`JSA_WINDOW`) before subsampling, so the side-lobes
  keep the survey's resolution. Change the two together. The mesh
  dimensions per panel are written to `data/jsa_mesh.tex`, which
  `fig5_jsa_maps.tex` inputs — never hard-code `mesh/cols`.
- The point order in `data/gvmcurve_*.csv`, `deff_*.csv`, `purity_*.csv`
  is along the level curve (sorted by signal wavelength), not by pump:
  every x–z family folds back in pump, so sorting by pump interleaves
  its two arms. Fixed 18 Sep 2026; the old ordering had silently dropped
  most of the oeo 0° family from figs 1 and 3.
- `data/` is derived output; never edit it by hand.

`make_figure_data.py` reads `../../bibojsa/survey` (the layout of the
project folder) and writes `data/`.
