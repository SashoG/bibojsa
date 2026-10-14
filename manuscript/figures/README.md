# Manuscript figures — pgfplots, drawn from the survey data

Every figure is LaTeX/pgfplots source, not a picture. The chain is:

```
bibojsa/survey/*.csv  --make_figure_data.py-->  data/*.csv  --fig<N>.tex-->  fig<N>.pdf
   (the survey)            (plot-ready slices)     (pgfplots)     (vector PDF)
```

| figure | source file | data files | manuscript label |
|---|---|---|---|
| 1 | `fig1_gvm_landscape.tex` | `data/gvmcurve_*.csv` (the `bound_*.csv` files are written but the boundaries are drawn as constants) | `fig:landscape` |
| 2 | `fig2_tpm_degenerate.tex` | `data/tpmdeg_*.csv` | `fig:tpmdeg` |
| 3 | `fig3_deff.tex` | `data/deff_*.csv` | `fig:deff` |
| 4 | `fig4_purity_families.tex` | `data/purity_*.csv` | `fig:familypurity` |
| 5 | `fig5_jsa_maps.tex` | `data/jsa_{A,B,C,D,E,R}.csv`, `data/jsa_mesh.tex` | `fig:jsa` |
| 6 | `fig6_filter.tex` | `data/filter_{A,B,C,D,E,R}.csv` | `fig:filter` |

## Building

```
make                         # builds all six PDFs from data/
python make_figure_data.py   # only if the survey changed: re-derives data/
```

`data/` (48 CSVs plus `jsa_mesh.tex`, ~3 MB) is kept in the folder. It is derived output —
`make_figure_data.py` regenerates it from `../../bibojsa/survey` in a few
seconds — but it is the exact input the six figures were drawn from, so
it is version-controlled with them rather than reconstructed on demand.
Re-derive it only when the survey data change.

Where the numbers come from, in order:

| layer | location | what it is |
|---|---|---|
| survey output | `../../bibojsa/survey/*.csv` | the computation itself: 1160 traced GVM points, the degenerate members, per-point purity and convergence grids |
| complete survey | `../data_file_1.zip` | the whole survey tree, including `landscape_field.csv` and the 162 per-point files under `points2/` that are too bulky to keep loose |
| plot-ready slices | `data/*.csv` | this folder: the columns each `\addplot` actually reads |

fig1–fig4 and fig6 compile with `pdflatex` (run twice — the shared
legends travel through the `.aux` file). **fig5 needs `lualatex`**: its
six |JSA| meshes (up to 150×150 cells) exceed pdflatex's static memory
and take a few minutes. On Overleaf, set the compiler to LuaLaTeX if you
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
  `cg0` / `cg45` / `cg90` for θ_PM = 0°, 45°, 90°; the interactions use a
  separate palette (`coeo` / `cooe` / `ceeo`, since 14 Oct 2026 purple /
  teal / brown) so that one colour never means both an interaction and a
  GVM condition (fig3 right panel: colour = interaction, line style = GVM
  condition). fig4 draws members with `inwin = 0` (idler beyond the
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
