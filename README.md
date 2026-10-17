# pcps-sim

Research compendium for a simulation study on the validity of inference based on
phylogeny-weighted species composition (phylogenetic fuzzy weighting, matrix **P**,
and Principal Coordinates of Phylogenetic Structure, PCPS), and on how tree shape
affects it. A secondary aim compares PCPS with phylogenetic community composition
(Peres-Neto & Kembel 2015), evoPCA (Pavoine 2016) and phylofactorization.

Author: Diogo B. Provete.

Status: simulations in progress. See `docs/PLANO_SIMULACAO.md` (simulation plan)
and `docs/MEMORIA_PROJETO.md` (project log), both in Portuguese.

## Repository structure

```
pcps-sim/
├── R/                             functions (sourced by the pipeline and tests)
│   ├── pfw.R                      matrix P, D_P, ADONIS and Mantel with site/taxa shuffles
│   ├── pcps_glm.R                 PCPS axes, GLM on PCPS axes, joint D_P tests
│   ├── pcc.R                      phylogenetic community composition (Peres-Neto & Kembel 2015)
│   ├── trees.R                    tree generator (beta-splitting topology, lineage-rate times)
│   ├── tree_grid.R                calibration of the gamma x Colless grid
│   ├── simulate_static.R          static Gaussian-niche metacommunities (validation)
│   ├── simulate_layerA.R          layer A scenarios (Pagel's lambda; signal in niche breadth and height)
│   ├── validation.R               replication of published simulations
│   └── pipeline.R                 design profiles and replicate chunks for targets
├── analysis/                      numbered scripts, run in order
│   ├── 01_validate_layerA.R       replication of Duarte et al. (2016) and Peres-Neto & Kembel (2015)
│   ├── 02_compare_published.R     comparison with published rejection rates
│   ├── 03_pnk_tree_spread.R       spread of PCC rejection rates across Yule trees
│   ├── 04_pnk_tree_shape.R        PCC rejection rates vs tree shape
│   ├── 05_tree_shape_empirical.R  gamma and Colless of anuran clades
│   ├── 06_tree_grid.R             simulated tree grid
│   ├── 07_pilot_tree_shape.R      tree-shape pilot
│   ├── 08_tree_stats.R            treestats statistics for grid and empirical trees
│   └── run_pipeline.R             runs the targets pipeline for one design profile
├── _targets.R                     targets pipeline (profiles: smoke, bench, pilot_glm, full)
├── tests/testthat/                unit tests
├── data/
│   ├── raw/                       immutable inputs (Portik et al. 2023 amphibian timetree)
│   ├── published/                 rejection rates transcribed from published tables
│   └── derived/                   simulation outputs (not tracked)
├── aws/                           scripts to run the pipeline on an EC2 spot instance
├── Dockerfile                     container used on AWS (rocker/r-ver + renv)
├── renv.lock                      package versions
├── docs/                          simulation plan, project log, method specifications
├── manuscript/                    Quarto manuscript (in preparation)
├── legacy/                        material from the 2012–2016 version of the project (read-only)
└── references/                    literature (not tracked)
```

## Reproducing the analyses

```r
renv::restore()
# quick end-to-end check
system("Rscript analysis/run_pipeline.R smoke 4")
# full study (about 1,100 CPU-hours; aws/README.md describes running it on AWS)
system("Rscript analysis/run_pipeline.R full 16")
```

Outputs are written to `data/derived/<profile>/results.csv`, and the targets store
to `_targets_<profile>/`.
