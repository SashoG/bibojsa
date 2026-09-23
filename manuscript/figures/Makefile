# Build the manuscript figures (vector PDFs) from the survey data.
#
#   make data     regenerate figures/data/*.csv from ../../bibojsa/survey
#   make          build all six figure PDFs
#   make clean    remove aux files
#
# fig1..fig4, fig6 compile with pdflatex (twice: legends go via .aux).
# fig5 needs lualatex (6 meshes x 40000 cells exceed pdflatex's memory);
# it takes a few minutes -- that is normal.

PDFFIGS = fig1_gvm_landscape fig2_tpm_degenerate fig3_deff \
          fig4_purity_families fig6_filter
LUAFIGS = fig5_jsa_maps

all: $(PDFFIGS:%=%.pdf) $(LUAFIGS:%=%.pdf)

$(PDFFIGS:%=%.pdf): %.pdf: %.tex bibofig-preamble.tex
	pdflatex -interaction=nonstopmode $*.tex
	pdflatex -interaction=nonstopmode $*.tex

$(LUAFIGS:%=%.pdf): %.pdf: %.tex bibofig-preamble.tex
	lualatex -interaction=nonstopmode $*.tex

data:
	python3 make_figure_data.py

clean:
	rm -f *.aux *.log *.out

.PHONY: all data clean
