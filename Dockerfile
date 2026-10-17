# Container for the simulation pipeline (targets + crew), used on AWS EC2.
# Build from the repository root:  docker build -t pcps-sim .
FROM rocker/r-ver:4.6.1

# Runtime libraries of the binary R packages (Posit Package Manager, Ubuntu)
RUN apt-get update && apt-get install -y --no-install-recommends \
      libglpk40 libgmp10 libxml2 libcurl4 libssl3 zlib1g \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /project
ENV RENV_PATHS_LIBRARY=/opt/renv/library \
    RENV_CONFIG_SANDBOX_ENABLED=FALSE
COPY renv.lock .Rprofile ./
COPY renv/activate.R renv/settings.json renv/
# rocker points `repos` to Linux binaries on Posit Package Manager
RUN R -q -e 'renv::restore(repos = getOption("repos"), prompt = FALSE)'

COPY R/ R/
COPY analysis/ analysis/
COPY data/raw/ data/raw/
COPY _targets.R dependencies.R ./

CMD ["Rscript", "analysis/run_pipeline.R", "full"]
