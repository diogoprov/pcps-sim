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
