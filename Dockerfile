# Container for the simulation pipeline (targets + crew), used on AWS EC2.
# Build from the repository root:  docker build -t pcps-sim .
FROM rocker/r-ver:4.6.1

# System requirements reported by renv::restore() (packages may be built from
# source when no Linux binary matches the lockfile version)
RUN apt-get update && apt-get install -y --no-install-recommends \
      cmake libcurl4-openssl-dev libglpk-dev libgmp3-dev libicu-dev libssl-dev \
      libuv1-dev libx11-dev libxml2-dev pandoc zlib1g-dev libpng-dev \
      libfontconfig1-dev libfreetype6-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /project
# rocker/r-ver links R to multithreaded OpenBLAS; with one crew worker per core,
# its threads only compete with the other workers (single-thread OpenBLAS is
# 2-3x faster on the tests, analysis/12_bench_blas.R)
ENV OPENBLAS_NUM_THREADS=1 \
    OMP_NUM_THREADS=1
ENV RENV_PATHS_LIBRARY=/opt/renv/library \
    RENV_CONFIG_SANDBOX_ENABLED=FALSE \
    MAKEFLAGS="-j8"
COPY renv.lock .Rprofile ./
COPY renv/activate.R renv/settings.json renv/
# rocker points `repos` to Linux binaries on Posit Package Manager. Only the
# packages the simulation pipeline uses (and their dependencies, from the
# lockfile) are restored: the lockfile also holds packages for the local
# analyses (clade recovery, formal models, figures), some from Bioconductor
# or needing graphics libraries, which the pipeline never loads.
RUN R -q -e 'renv::restore(repos = getOption("repos"), prompt = FALSE, packages = c("targets", "tarchetypes", "crew", "qs2", "ape", "phytools", "geiger", "treestats", "vegan", "RSpectra", "dplyr", "tidyr", "tibble", "purrr", "readr", "data.table"))'

COPY R/ R/
COPY analysis/ analysis/
COPY data/raw/ data/raw/
COPY data/calibrated/ data/calibrated/
COPY _targets.R dependencies.R ./

CMD ["Rscript", "analysis/run_pipeline.R", "full"]
