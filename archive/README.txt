ARCHIVE - superseded files of the BIBO_JSA project
==================================================

Reorganised 14 October 2026 (manuscript at v0.9). Nothing in here is
read by manuscript_oe.tex, build.bat, the figure Makefile, the tests or
the example scripts, and nothing in here is part of the public git
repository (the root .gitignore whitelists only the release files, so
this folder is ignored automatically). Deleting it would break nothing;
it is kept for provenance.

notes_2026-07/      early working notes, superseded by the revision
                    notes and the README
  BIBO_JSA_state_of_the_art.md   literature survey, 7 Jul 2026
  claude_code_prompt_bibojsa.md  original brief for the bibojsa code, 8 Jul
  PROGRESS.md                    July progress log (milestones 1-5, the
                                 fixed-pump survey); last edited 24 Jul

reviews/            all six review rounds, answered by manuscript
                    v0.3-v0.9. No review is open as of v0.9.
  REVIEW_2026-09-04.md     first review (v0.2)
  REVIEW_2026-09-14.md     second review (v0.3)
  REVIEW_2026-09-14b.md    third review (v0.4, novelty/physics)
  REVIEW_2026-09-18.md     fourth review (v0.5.1, OE novelty + numerics)
  REVIEW_2026-09-25.md     fifth review (answered by v0.8)
  REVIEW_2026-10-14.md     sixth review (answered by v0.9)
  review_2026-09-18/       figure review working files: before/after
                           PDFs, independent_check.py + output,
                           make_figure_data_fixed.py and its diff (fix
                           already merged into figures/make_figure_data.py).
                           SUBMISSION_CHECKLIST_OE.md still cites
                           independent_check.py at its old path.
  review_2026-10-14/       numerics audit of the sixth review
                           (NUMERICS_AUDIT.md, audit_part1/2 scripts +
                           output) and the fig4/fig6 proposals, which v0.9
                           adopted into manuscript/figures/ (the adopted
                           versions differ only in comments and label
                           placement). Moved in a second pass, 14 Oct.

superseded_data/
  data_file_1.zip          Prism "Data File 1" of 18 Sep. STALE since
                           v0.9 (headline CSVs regenerated 14 Oct) and
                           superseded by the Zenodo release of
                           bibojsa/survey/. Kept in case a referee asks
                           for a supplementary copy - regenerate from
                           bibojsa/survey/ if so.
                           (manuscript/figures/README.md still points to
                           ../data_file_1.zip; that row changes to the
                           Zenodo DOI in v0.10.)
  bibojsa_github_v1.0.0.zip  22 Sep repository snapshot used once for the
                           first push to github.com/SashoG/bibojsa (step 1
                           of github_release/ZENODO_RELEASE_STEPS.md, done).
                           Predates the v0.9 code/data changes; the
                           release is now cut from the git repository.
  figures_png/             July matplotlib PNGs of the fixed-pump survey
                           (fig1_pm_loci ... fig6_purity_scans_top3),
                           numbered differently from the paper's figures.
                           Was bibojsa/survey/figures_png/ (gitignored).

optica_template/    was manuscript/archive/ (archived 14 Sep): Optica's
                    demo files (Optica-template.tex, opticafig1/2.pdf,
                    sample.bib, opticajnlnt.bst) and the July
                    figures/README. Its own README.txt explains each.

papers/
  13. Independent high-purity photons ... .pdf
                    arXiv:1712.07140v2 preprint of Graffitti et al.
                    (8 pp, includes the supplementary material). The
                    published Optica 5, 514 (2018) version stays in
                    papers/ as "13.Independent ...".

Deliberately NOT archived
-------------------------
- bibojsa/survey/survey_*.csv, convergence_all.csv, points/ and the
  examples that wrote them: superseded fixed-pump survey, but point R of
  Table 2 is drawn from them and they are in the public repository.
- Everything tracked by git (bibojsa/, manuscript/figures/,
  manuscript/exploration_2026-09-21/, README, LICENSEs, CITATION.cff,
  .zenodo.json).
- optica-article.cls, opticajnl.bst, styles/, jabbrv* : required to build.
- manuscript_oe.aux/.bbl/.blg/.log : transient, build.bat rewrites them.
- submission_build.zip : current (regenerated from v0.9).
- github_release/ZENODO_RELEASE_STEPS.md : steps 2-6 still open.
- BIBO_JSA.code-workspace : VS Code workspace file.
