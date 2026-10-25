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

### 2026-10-25 — treestats issue #35
- Read the maintainer's reply to the author's proposal (treestats#35) to add the
  tree-shape simulator to treestats. The maintainer declined (out of scope) and
  suggested a separate package. The author decided to release it as a separate
  package; the assistant recorded the decision and listed the open choices.
- The author chose the name treeshapesim, GitHub plus r-universe only, and to
  keep the local copy in pcps-sim. The assistant built the package (code and
  tests copied from R/trees.R, roxygen documentation, README) and ran R CMD
  check, and drafted a reply to the issue for the author to post.
- The author created github.com/diogoprov/treeshapesim and pushed the package (the GitHub integration could not create repositories).
- Opened treeshapesim#1 (pkgdown site via GitHub Actions/Pages, keywords) and diogoprov.r-universe.dev#2 (registry entry); drafted four hex-sticker concepts.
