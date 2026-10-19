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
│   ├── simulate_layerB.R          layer B process-based dynamics (Thompson et al. 2020), pool calibration
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
│   ├── 09_diag_s2sigma.R          diagnosis of the GLM inflation in scenario S2sigma
│   ├── 10_calibrate_layerB.R      layer B pool sizes giving realized N = 50/100/200
│   ├── 11_pilot_layerB.R          layer B pilot
│   ├── 12_bench_blas.R            BLAS timing (Docker)
│   ├── 13_fig1_schematic.R        Figure 1 (simulation schematic)
│   └── run_pipeline.R             runs the targets pipeline for one design profile
├── _targets.R                     targets pipeline, layers A and B (profiles in R/pipeline.R)
├── tests/testthat/                unit tests
├── data/
│   ├── raw/                       immutable inputs (Portik et al. 2023 timetree; not redistributed, see data/raw/README.md)
│   ├── published/                 rejection rates transcribed from published tables
│   ├── calibrated/                calibrated design inputs (layer B pool sizes)
│   └── derived/                   simulation outputs (not tracked)
├── aws/                           scripts to run the pipeline on an EC2 spot instance
├── Dockerfile                     container used on AWS (rocker/r-ver + renv)
├── renv.lock                      package versions
├── docs/                          simulation plan, project log, method specifications,
│                                  status report (resumo_projeto.qmd; quarto render)
├── manuscript/                    Quarto manuscript (in preparation) and figures
├── legacy/                        material from the 2012–2016 version of the project (not tracked)
└── references/                    literature (not tracked)
```

## Reproducing the analyses

Place the Portik et al. (2023) timetree in `data/raw/` first (see
`data/raw/README.md`).

```r
renv::restore()
# quick end-to-end check
system("Rscript analysis/run_pipeline.R smoke 4")
# full study (about 1,100 CPU-hours; aws/README.md describes running it on AWS)
system("Rscript analysis/run_pipeline.R full 16")
```

Outputs are written to `data/derived/<profile>/results.csv`, and the targets store
to `_targets_<profile>/`.

## License

Code is released under the MIT License (`LICENSE`). Text, figures and simulation
outputs are released under CC BY 4.0 (`LICENSE.md`).

## Citation

See `CITATION.cff`.
