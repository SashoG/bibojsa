# bibojsa — group-velocity matching and photon-pair purity in BiB₃O₆

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Code and data underlying

> A. Gaydardzhiev, *Group-velocity-matched photon pairs in BiB₃O₆: a map of
> nondegenerate purity loci across all twelve interactions*, submitted to
> Optics Express (2026).

`bibojsa` computes collinear spontaneous parametric down-conversion in
biaxial crystals: Sellmeier dispersion and group indices, biaxial
phase matching from the Fresnel equation of wave normals, the effective
nonlinearity of point group 2 (the full monoclinic tensor), classification
of the group-velocity-matching (GVM) conditions as level curves in pump
wavelength, and the joint spectral amplitude with its Schmidt decomposition
and heralded-photon purity. The survey data, the plot-ready figure slices
and the pgfplots figure sources of the manuscript are included, so every
number and every figure in the paper can be regenerated from this
repository.

## Layout

```
bibojsa/                 the Python project
  pyproject.toml         pip install -e bibojsa   (numpy, scipy, matplotlib)
  bibojsa/               package: sellmeier, biaxial, phasematch, deff, gvm,
                         jsa, landscape, materials/{bibo,kdp,ktp}
  tests/                 108 pytest tests (KDP and KTP literature checks,
                         BIBO Sellmeier, d_eff, GVM, EPM/JSA milestones,
                         both halves of the type-II surfaces)
  examples/              the survey scripts that wrote survey/
  notebooks/             validation and survey notebooks (.ipynb + jupytext .py)
  survey/                THE DATA — every CSV read by the manuscript;
                         survey/README.txt documents each file and column
  SURVEY_FULL_2026-10-17.md  what the 17 Oct 2026 full survey added and why
manuscript/
  figures/               pgfplots sources (fig1, fig3–fig6 = manuscript
                         Figs. 1–5), Makefile, make_figure_data.py, and
                         data/ (the plot-ready slices)
  exploration_2026-09-21/  independent-code checks reported in Sec. 4.5 and
                         Sec. 5.3: off-principal-plane GVM loci and the
                         absolute pair rate at point A
CITATION.cff, .zenodo.json   citation metadata
LICENSE                  MIT (code)
LICENSE-DATA             CC BY 4.0 (data: survey/, figures/data/, exploration CSVs)
```

The manuscript text itself is not in this repository.

## Reproducing the paper

```bash
git clone https://github.com/SashoG/bibojsa
cd bibojsa/bibojsa
pip install -e .[dev]
python -m pytest tests                     # 108 tests, ~1 min

python examples/landscape_survey.py all    # survey/landscape_*.csv (~hours)
python examples/full_landscape_preview.py  # merged curves, survey/figdata/
python examples/landscape_points.py        # survey/points2/, points2_summary.csv
python examples/landscape_headline.py all  # Table 1, headline_*.csv, family_purity.csv
python examples/window_tail.py             # window_tail.csv, length_extension.csv
python examples/add_1nm_rows.py            # 1-nm pump rows of length_extension.csv (G, H)

cd ../manuscript/figures
python make_figure_data.py                 # survey → data/*.csv (seconds)
make                                       # figure PDFs (fig5 needs lualatex + luatex85)
```

The committed `survey/` tree is the exact output the manuscript was drawn
from (regenerated 17 October 2026, manuscript v0.11, 1-nm length rows added
20 October 2026: all twelve
principal-plane interactions, i.e. both halves of every type-II
phase-matching surface, x–y plane included), so the long survey steps can
be skipped; `make_figure_data.py` and `make` reproduce the five figures
from it in minutes. The eight-interaction run of 18 September 2026
(manuscript v0.6–v0.10) is not part of this repository; its rows are a
strict subset of the present tables.

## Conventions (short form; full text in `bibojsa/survey/README.txt`)

Interactions are labelled by the LPR2010 o/e sequence of
(λ₁, λ₂, λ₃) = (idler, signal, pump), λ₁ ≥ λ₂ > λ₃. θ is the crystal-cut
angle from *z* in the *y–z* and *x–z* planes and from *x* in the *x–y*
plane. θ_PM is the JSA ridge orientation measured from the signal axis:
0° means k′_p = k′_s, 90° means k′_p = k′_i, +45° is symmetric GVM
(extended phase matching). Purities are Schmidt purities of the complex
JSA, and every filtered purity is reported together with the fraction of
pairs the filter transmits. The idler-transparency flag uses BIBO's
1 cm⁻¹ absorption onset at 2.60 µm, not the 3.09 µm Sellmeier validity
limit.

Since 17 October 2026 a type-II sequence also fixes which half of the
phase-matching surface a row lies on: `oeo` (x–z) is the half on which the
e daughter is the shorter one, `eoo` the half on which it is the longer
one, and likewise `eoe`/`oee` for the z-polarized wave in the x–y plane;
the GVM level curves are continuous across the two halves
(`bibojsa/survey/README.txt`).

## Citing

Please cite the paper and, for the code and data, the Zenodo record
(`CITATION.cff` has both in machine-readable form):

```
Gaydardzhiev, A. (2026). bibojsa: group-velocity matching and photon-pair
purity in BiB3O6 (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX
```

## Licence

Code (everything under `bibojsa/bibojsa/`, `tests/`, `examples/`,
`notebooks/`, `manuscript/figures/*.tex`, `*.py`, `Makefile`) is released
under the MIT licence (`LICENSE`). Data (`bibojsa/survey/`,
`manuscript/figures/data/`, and the CSV/PNG/TXT outputs in
`manuscript/exploration_2026-09-21/`) are released under
Creative Commons Attribution 4.0 International (`LICENSE-DATA`).

Sellmeier coefficients and nonlinear tensor values are taken from the
published sources cited in the module docstrings (Miyata–Umemura–Kato,
Hellwig et al., Umemura et al.; Petrov et al., Laser Photonics Rev. 4, 53
(2010)).
