# AI interaction log

Record of the use of a generative AI assistant in this project, kept for
transparency and for the AI-use statement of the manuscript. The author
(D. B. Provete) designed the study, took every decision recorded here and is
responsible for the content. Entries are appended at the end of each working
session.

**Tool.** Claude (Anthropic), model Claude Opus 5.5 (`claude-opus-5-5`), used
through Claude Code (command-line agent) with access to this repository, a
local R installation and, when authorized, an AWS account and web search tools.

**How the assistant was used**
- Writing and refactoring R code (simulation, tests, pipeline), always with
  unit tests and two independent automated code reviews per change (coding
  standards; agreement with the plan), whose findings were fixed before use.
- Running simulations (locally and on AWS, after the author's approval of each
  run and its estimated cost) and summarizing their outputs.
- Drafting internal reports (in Portuguese) and manuscript text (in English).
- Searching and summarizing the literature. Rule: no reference, DOI or number is
  cited unless verified in a real search, a PDF supplied by the author or the
  project's notes; unverified items are marked `[CITAR: ...]`.

**What the assistant did not do**
- Choose the research question, the hypotheses or the final interpretation.
- Add references that were not verified.

**Verification.** Numbers reported in documents were checked against the
result files by an independent review step; corrections are recorded in the
entries below and in `docs/MEMORIA_PROJETO.md` (which holds the detailed log of
decisions, commands and results). The full history of changes is in git.

## Entries

### 2026-10-16 — Reassessment and plan
- Assessed whether the 2016 chapter overlapped with Duarte et al. 2016
  (literature search, reading of the published article and appendices);
  conclusion: no overlap in the core question (tree shape).
- With the author, reframed the study (validity of inference on P under
  controlled tree shape; method comparison as second aim) and drafted
  `docs/PLANO_SIMULACAO.md`. Project moved to this git repository.

### 2026-10-17 — Implementation
- Implemented the tree-shape grid (γ × Colless, calibrated to empirical anuran
  clades from Portik et al. 2023), the static layer (Peres-Neto et al. 2012
  protocol), the D_P tests (ADONIS, Mantel, GLM on PCPS axes), the PCC
  (re-implemented from the paper) and the process-based layer B (own
  implementation of the Thompson et al. 2020 equations). Validation against
  published results; pilots.

### 2026-10-18 — Full runs (AWS)
- Ran layer A (`full`, `axes`) and layer B on AWS spot instances after the
  author's approval; fixed a failure in layer B (pool-tree rejection sampling)
  and re-ran the failed chunks. Figure 1 (schematic) drafted and approved by
  the author.

### 2026-10-19 — Analysis, objective 2 and formal models
- Status report (`docs/resumo_projeto.qmd`); GAM of layer B on the realized
  shape of the analysed tree (revised three times after reviews: per-stratum
  fits, no interpolation across gaps, pruning effect separated from the pool).
- Author's decisions: evoPCA and phylofactorization included; planted-clade
  recovery metric; no second simulator; regional pools from the author's
  empirical metacommunities; acknowledgements (the simulation idea is the
  author's).
- Implemented the planted-clade objective (scores per node for PCPS, PCC,
  evoPCA and binomial/compositional phylofactorization, validated against the
  original packages), ran it locally (864,000 rows).
- Formal beta-binomial models (`brms`) of type I error, power and clade
  recovery, with tree shape as 3-level factors (accepted by the author);
  convergence and posterior predictive checks reported. Hypothesis verdicts
  were revised twice after reviews flagged over-claiming.
- Started this log at the author's request.
- Author's decisions for the manuscript: one article with layers A and B, layer
  B as a second main objective (methods under metacommunity assembly
  processes) with its own formal model; regional pools from the author's
  empirical metacommunities, names harmonized (AmphiNom); Methods and Results
  as .docx following DiRenzo et al. 2023 and Williams et al. 2024.
- Fitted the layer B formal model; reviews found missing N interactions (fit
  misses concentrated at N = 50 and 200) and a stale cached fit reused after a
  sampler change; both fixed and the model refitted.
- Built the regional pools: names from the author's 14 empirical tadpole
  metacommunities resolved (codes, spelling) and harmonized with the ASW
  (AmphiNom), matched to the Portik tree; three nested pools (60, 81, 123
  species). The author set the rules (exclude undetermined taxa, no manual
  review) and resolved the doubtful names; reviews caught bugs in code
  resolution and an over-claimed bootstrap interval, fixed before reporting.

### 2026-10-20 — Case study and manuscript draft
- Applied the tests and clade-scoring methods to three published tadpole
  metacommunities chosen by the author's criteria (≥ 10 ponds, depth
  measured); the author identified the sources (Provete et al. 2014; Both et
  al. 2011; Prado & Rossa-Feres 2014) and pointed to his 2018 preprint, whose
  gradient (canopy cover) was added for Bocaina.
- Drafted Methods, Results, Box 1 and Supporting Information in the author's
  voice, following DiRenzo et al. 2023 and Williams et al. 2024; every number
  was checked by independent review passes, which found wrong values and
  over-claims that were corrected. References were verified on Crossref;
  unverified ones are marked [CITAR: ...].

### 2026-10-20 — Manuscript restructure after the author's read of the docx
- At the author's request: formulas and subscripts converted to native Word
  equations (LaTeX math rendered by pandoc as OMML); Methods and Results
  shortened (about 2,800 and 1,500 words) by moving generator, Layer B,
  node-score, model-formula, prior, computation, validation, Layer B null-rate,
  abundance-recovery and convergence details to Appendices S1–S2 of the
  Supporting Information; Results subheadings removed.
- Figures: readable palettes (YlGnBu with adaptive label colour; Okabe–Ito for
  Fig. 4 and Box 1), several short biplot arrows in the RDA icon of Fig. 1, and
  Box 1 points drawn at the bootstrap median with the point estimate as an
  open diamond.
- Two review passes (Standards, Spec) checked that no number was lost in the
  move; they caught a wrong pointer (N = 200 vs N = 100), six dropped facts and
  missing SI references, all fixed. Appendix S1 now states three model details
  taken from the code (σ = √2 σ_T, multinomial immigration, final state only).

### 2026-10-20 — Submission material, anonymization and Discussion
- Figures exported to the Wiley specifications (vector PDF and 600 dpi PNG,
  80–180 mm, Figure_N) in manuscript/figures/submission/.
- GitHub release v0.1.0 archived on Zenodo (concept DOI 10.5281/zenodo.23500104),
  after the author enabled the integration.
- Numbered abstract, keywords, title page and data statement drafted; the
  author's name and identifying links removed from the manuscript for
  double-anonymous review; tools/make_review_bundle.sh builds an anonymized
  zip of the repository for an anonymous Zenodo deposit to be made by the author.
- Discussion drafted from three literature reviews run by subagents with
  Scholar Gateway, Undermind, Consensus and OpenAlex (docs/lit_review/), plus
  two Scopus/LeapSpace reports supplied by the author; every reference was
  verified on Crossref. An accuracy review corrected 18 over-claims or
  attributions beyond what the sources say.

### 2026-10-20 — Mock review and responses
- A subagent with no access to the conversation wrote a mock MEE review
  (docs/review/MEE_mock_review.md; major revision). The author chose which
  items to address; text-only items were answered at once (title, γ × N,
  analysed-tree shape in Layer B, clade-recovery ranges, Procrustes variant).
- New analyses requested by the author: prior sensitivity (priorsense),
  prior predictive checks, H1 equivalence, paired contrasts, PCPS2 arch
  effect, RDA PCPS in Layer A (which changed one recommendation), clade age,
  δ = 0 baseline (27-min local run), and the lack of fit of the mean structure
  (random effect per condition). All numbers were taken from the scripts'
  CSV outputs and checked before entering the text.
- Two new simulations designed with the author's approval and run on AWS spot
  instances: optima with a fixed marginal distribution (separates tree shape
  effects on the optima and on the analysis) and a spatial Layer B (torus,
  autocorrelated gradient, kernel dispersal, torus-translation null validated
  in a local pilot, incidence and abundance). Spot interruptions were handled
  by aws/supervise.sh, which resumes from the S3 checkpoints.

### 2026-10-21 — AWS results, terminology, accuracy review and submission files
- Both AWS runs finished after several spot interruptions handled by
  aws/supervise.sh. With optima of fixed marginal distribution, the loss of
  power with low γ persisted (γ acts mainly through the analysis); in the
  spatial process-based simulations, mass effects reduced power much less, and
  the PCC kept its power with abundance but not incidence data. These results
  were summarized by analysis/36–37 and written into the text and Appendix S2.
- At the author's request, "Layer A/B" became "static" and "process-based
  simulations" throughout the reader-facing text and Fig. 1.
- A third accuracy review (subagent) found over-claims of scope (uniform
  dispersal, the PCC with abundance data, γ acting through the analysis) and
  literature attributions beyond what was read; all were corrected, and the
  text was trimmed to 7,996 words.
- Built a preprint (tools/build_preprint.py; single docx with the title page,
  not anonymized, citing the public repository) and a submission copy of the
  cover letter with reviewers' e-mails found on institutional pages
  (tools/render_cover_letter.py); e-mails are kept out of the public repository.
  Two reviewers' e-mails were not found and were left blank rather than guessed.

### 2026-10-21 — AI disclosure statement (AIdIT)
- At the author's request, the free-form AI paragraph of the title page was
  replaced by an "AI disclosure statement" in the AIdIT v0.1 format (Drobniak et
  al. 2026), built from the template of the official app
  (github.com/szymekdr/AI_reporting, code/app.R) and placed after the
  references of the main text. The author confirmed the engines, areas of use,
  dates and add-on answers, asked for the singular ("The author declares"), and
  asked to state compliance with Portaria CNPq nº 2.664/2026 (details taken from
  the author's page provetelab.org/ia-escrita-cientifica) and to point to this
  log on GitHub; the link appears only in the preprint, since the submission
  text is anonymized.

### 2026-10-22 — Codex methodological review and correction plan
- Codex (OpenAI, GPT-6-based assistant) reviewed the simulation code, manuscript and existing outputs at the author’s request, adopting the perspective of a simulated Methods in Ecology and Evolution reviewer. All existing testthat tests passed with the local project dependencies; the full simulation and Bayesian models were not rerun.
- Identified eight methodological concerns, distinguishing confirmed implementation choices and output checks from questions about null hypotheses and causal interpretation. Proposed targeted controls, sensitivity analyses and revisions, with a preliminary estimate of 50–85 technical work hours plus benchmark-dependent computation.
- At the author’s explicit request, saved the critique and plan in docs/review/CODEX_critica_e_plano_2026-10-22.md, identified Codex as its author and committed the documentation. No proposed simulation corrections were implemented, and no scientific decisions accepting the proposed corrections were inferred.

### 2026-10-22 — Point 1 of the Codex review: estimands of S2σ/S2h (Claude)
- Claude (Opus 5.5, Claude Code) checked the eight points of the Codex review against the code and judged point 1 the most important. Under S2σ and S2h, the taxa-shuffle null hypothesis (species exchangeable on the tree) is false, so the null hypothesis of the max test is false too, and a rejection there is not a type I error. It is an error only relative to "optima without phylogenetic signal", which no tip permutation can isolate.
- Re-read the existing diagnostic output (data/derived/diag_s2sigma/diag_s2sigma.csv; no new simulation). The correlation of the CWM of σ with E² was the same in S2 and S2σ (0.87–0.93). Signal in σ instead raised the linear correlation of PCPS2 with E (low γ: 0.29 → 0.44 at N = 100 and 0.34 → 0.59 at N = 200). The quadratic GLM rejected less often than the linear one, which corrected an earlier statement by Claude in the same conversation and a sentence in Appendix S2.
- Proposal in docs/review/ponto1_proposta_estimandos.md. At the author's request it was then applied to the main text (Scenarios and estimands with a new column, abstract items 3–4, results, Fig. 2 caption and colour-scale label, discussion), Appendix S2 and the cover letter, with H2 kept as originally stated (author's choice). No number changed. The author approved the plan and asked for the commit; the scientific decision on the reframing is the author's.
- Word limits (same day): after the reframing, the main text was about 8,450 words and the abstract 367. Claude proposed 19 cuts (docs/review/cortes_palavras_2026-10-22.md): the null hypotheses moved to Methods compared, a sentence repeated in the Results, a shorter Discussion paragraph on S2σ, and supporting sentences whose citations appeared only there (Diniz-Filho et al. 2012, Molina-Venegas & Roquet 2014, Ai & Ellwood 2022, whose references were removed). The author approved all of them. The result was about 7,965 words with references (same rule as the earlier 7,996 count, AI disclosure statement excluded) and an abstract of 346 words. No result number was changed. The preprint and the submission copy of the cover letter were rebuilt.
- Points 4 and 5 of the Codex review, text only (same day, at the author's request). Point 4: the conclusion that γ acts "mainly through the analysis" was limited to "not explained by the marginal distribution of the optima alone", because in the fixed-marginal run related species still had similar optima (main text and Appendix S1/S2). Point 5: the main text now states that "PCPS" results refer to the GLM on single axes without Procrustes, and the RDA PCPS power contrast moved from Appendix S2 into the Results. A leftover "null scenario" (H1 paragraph) became "scenario from S0 to S2h". Small compensating cuts kept the text at about 7,966 words.
- Point 7 of the Codex review (same day, at the author's request): Claude wrote analysis/38_raw_contrasts.R. It recomputes every model contrast from the raw rates with Monte Carlo intervals that respect the design (unpaired across grid cells and archetypes, paired by tree for S2σ/S2h against S2) and lists the conditions where the beta-binomial and random-effect models differ. Two review sub-agents (standards and spec) checked it. The fixes were Jeffreys-based standard errors for rates of 0 or 1, a sign-aware agreement criterion, per-condition Monte Carlo intervals, and clearer names. Findings: the model contrasts behind the main claims agree with the raw ones, except that the beta-binomial shrinks the S2σ increases (e.g., GLM PCPS2 at N = 200 with low γ: raw +0.30 against model +0.21) and understates the largest imbalance loss of power (raw 0.13 against 0.09). Its per-condition medians fall outside the raw 95% interval in about a third of the power and clade conditions (by up to 0.13), which matters for the figures, which show those medians. No manuscript text was changed at this step.
- Point 7 applied (same day; the author chose option A). Figs 2–5 now show raw rates instead of beta-binomial medians: Fig. 4 shows raw rates with 95% Jeffreys intervals, and analysis/22_manuscript_figures.R reads aggregated.csv. In the text, Methods treats raw rates and Monte Carlo contrasts as the primary evidence. Results report the S2σ increases as raw paired contrasts, the raw H1 range, a PERMANOVA power that did not decrease with low γ, and a maximum imbalance loss of 0.13. The sentence that presented the quasibinomial dispersion as a justification was removed, because that dispersion is lack of fit of the means. Appendix S2 gained a section on raw contrasts and model fit. To keep the text consistent with the new figures, the other model-based numbers in the abstract and results (γ effect on power, mass effects, competition, clade recovery) were also replaced by raw values. The text is at about 7,982 words and the abstract at 346.
- Appendix S2, "Details of the main results" (same day, at the author's request): the remaining model-based numbers were replaced by raw values (raw contrasts from raw_contrasts.csv, raw rates from the aggregated files) so that they match the figures. Static power by N, Mantel and PERMANOVA ranges, the baselines and losses of the process-based simulations, the range of the competition effect, and the clade-recovery contrasts changed by at most a few hundredths. The one relevant change was the S2σ increase of GLM PCPS2 with abundance data: 0.04 (paired raw) instead of 0.08 (model).
- Point 6 of the Codex review (same day, plan docs/review/ponto6_plano.md approved by the author: N = 100, incidence, h × 0.2, up to 4 workers). Claude added `objective = "interaction"` to phylofactor_node_scores() and `h_factor` to sim_metacommunity_clade(), both with tests, and with defaults that keep the main run unchanged. It wrote analysis/39_clade_sensitivity.R: 24 conditions × 1,000 replicates on the trees and seeds of the main run, 35 min on 4 workers. The run reproduced the main run exactly (hit, rank and Jaccard). Two review sub-agents checked the code. Their spec findings were fixed in the summary step (counts of NA, intervals at 0/1, a ranking table). Results: evoPCA and PCPS do not depend on the number of axes (evoPCA always chose axis 1), and the PCC changes by at most 0.06. The full binomial objective of phylofactorization recovered a clade that differed only in prevalence in 0.07–0.12 of replicates. The interaction-only objective did so in 0.001–0.004, with equal or slightly higher recovery of planted clades. No manuscript text was changed at this step.
- Point 6 applied (same day, approved by the author). Main text: one sentence in Methods (the sensitivity analysis), one in Results (axes and the phylofactorization objective). The Discussion and recommendation (v) now require the phylofactorization objective to be restricted to the interaction with the gradient. Four compensating cuts (two supporting sentences and two repetitions) keep the text at about 7,982 words. Appendix S1 gained the method paragraph (rank and Jaccard moved there) and Appendix S2 the section "Sensitivity of clade recovery". The docx, preprint and submission copy of the cover letter were rebuilt.

### 2026-10-23 — Contributions to treestats and PCPS (Claude)
- At the author's request, Claude prepared two contributions with roxygen help and patches that apply cleanly to the current GitHub clones. The first is `sim_tree_shape()` for treestats, from R/trees.R; 15 tests pass; building treestats locally failed at the linker for lack of the CRAN gfortran toolchain, so the function was tested against the installed treestats. The second is `permut.torus()` and a `permutation.site` argument for PCPS's pcps.sig(), matrix.p.sig() and matrix.p.null(). It was checked on a patched local build: the results are unchanged by default, and with the argument there is one null statistic per permutation. The workflow for the pull requests is in contrib/FLUXO_PRS.md. No reference for the torus-translation test was verified, so the help leaves a TODO placeholder instead of a citation.
- 2026-10-26, correction prompted by a literature search for the treeshapesim package (Consensus, Undermind, Scholar Gateway; the laser source code was checked). The branching-time process of the tree generator, a pure birth in which each lineage splits at rate k^c, is the pure-birth case of the density-dependent exponential (DDX) model of the R package laser (per-lineage rate r·k^(−x), so c = −x). Appendix S1 called the generator "not a published model"; at the author's request, the Methods and Appendix S1 now cite Rabosky (2006, laser; DOI checked on Crossref) and say that only the combination with β-splitting topologies is our design. Compensating cuts in the Methods overview keep the text at about 7,983 words.
- 2026-11-06, attribution correction prompted by a question from a Claude session on the comdynamite package. Claude read Peres-Neto et al. (2012; PDF supplied by the author), Peres-Neto & Kembel (2015; Appendix A) and Duarte et al. (2016; Appendix S2). The Poisson draw and h ~ U(0, 30) come from the two later studies, not from the 2012 protocol, so the Methods now cite them. Code comments, the plan (dated note) and the reporting guide were corrected; the simulations and results did not change. The author approved the text before it was applied.
