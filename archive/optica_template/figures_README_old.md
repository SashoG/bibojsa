# Manuscript figures — pgfplots, drawn from the survey data

Every figure is LaTeX/pgfplots source, not a picture. The chain is:

```
bibojsa/survey/*.csv  --make_figure_data.py-->  data/*.csv  --fig<N>.tex-->  fig<N>.pdf
   (archived data)        (plot-ready slices)      (pgfplots)      (vector PDF)
```

| figure | source file | data files | manuscript label |
|---|---|---|---|
| 1 | `fig1_pm_loci.tex` | `data/loci_*.csv`, `data/crossings_*.csv` | `fig:loci` |
| 2 | `fig2_theta_pm.tex` | `data/loci_*.csv`, `data/crossings_*.csv` | `fig:thetapm` |
| 3 | `fig3_deff.tex` | `data/loci_*.csv` | `fig:deff` |
| 4 | `fig4_summary_ranking.tex` | `data/ranking.csv` | `fig:ranking` |
| 5 | `fig5_jsa_maps.tex` | `data/jsa_point{45,31,42,2,28,14}.csv` | `fig:jsa` |
| 6 | `fig6_purity_scans.tex` | `data/purity_point{45,31,42}_vs{L,bw}.csv` | `fig:purityscans` |

## Building

```
make data   # only if the survey CSVs changed (re-derives data/*.csv)
make        # builds all six PDFs
```

fig1–fig4 and fig6 compile with `pdflatex` (run twice — the shared
legends travel through the `.aux` file). **fig5 needs `lualatex`** and
takes a few minutes: its six 140×140 |JSA| meshes exceed pdflatex's
static memory. On Overleaf, set the compiler to LuaLaTeX if you compile
the figures there.

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
  `bibofig-preamble.tex`.
- Per-figure layout (panels, ranges, legends, annotations): the
  `fig<N>.tex` files — plain pgfplots, every `\addplot` names the CSV
  it draws.
- If the survey data change, rerun `make data` (it prints the fig4
  ranking labels; if the ranking order changed, paste the new list into
  `fig4_summary_ranking.tex`'s `yticklabels`).
- `data/` is derived output; never edit it by hand.

`make_figure_data.py` defaults to reading `../../bibojsa/survey`
(i.e. the layout of the project folder); override with
`--survey PATH --out data`.
