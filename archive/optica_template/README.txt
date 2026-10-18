ARCHIVE — files no longer used by the manuscript
================================================

Nothing in here is referenced by manuscript_oe.tex, build.bat or the
figure sources. It is kept only so the provenance of the folder stays
clear. Deleting the whole archive folder would break nothing.

Reviewed 2026-09-14.

From the Optica universal template package
------------------------------------------
These shipped with the template download and belong to Optica's own
demonstration document, not to this paper.

  Optica-template.tex   Optica's sample manuscript, showing how to use
                        the class. Useful as a reference for markup we
                        have not used (multi-author blocks, tables,
                        supplementary-material syntax); not compiled by
                        anything here.
  opticafig1.pdf        the two sample figures that Optica-template.tex
  opticafig2.pdf        includes. Ours live in figures/.
  sample.bib            the sample bibliography for that template. Our
                        bibliography is references.bib.
  opticajnlnt.bst       the "no title" variant of the BibTeX style, which
                        omits article titles from the reference list.
                        optica-article.cls hard-codes
                        \bibliographystyle{opticajnl} (line 186), so this
                        variant is never loaded. Move it back and set the
                        style by hand only if a journal asks for
                        title-less references.

Superseded project files
------------------------
  figures_README_old.md the figures/README.md from 24 July 2026. It
                        documents the pre-v0.3 figure set — fig1_pm_loci,
                        fig2_theta_pm, fig4_summary_ranking,
                        fig6_purity_scans — which the September rewrite
                        replaced with fig1_gvm_landscape,
                        fig2_tpm_degenerate, fig4_purity_families and
                        fig6_filter. figures/README.md now describes the
                        current set.

What was deliberately NOT archived
----------------------------------
  optica-article.cls, opticajnl.bst, styles/opticajournal.sty
        required to compile; the class errors without them.
  jabbrv.sty, jabbrv-ltwa-en.ldf, jabbrv-ltwa-all.ldf
        required: optica-article.cls does \RequirePackage{jabbrv}
        (line 252) to abbreviate journal names in the bibliography.
  manuscript_oe.aux/.bbl/.blg/.log
        transient build output. build.bat deletes and regenerates them
        on every run, so they need no housekeeping.
  data_file_1.zip
        the Data File 1 upload for Prism (the whole bibojsa/survey tree).
  submission_build.zip
        the flat submission package. Regenerate it whenever
        manuscript_oe.tex changes, or it will go stale.
