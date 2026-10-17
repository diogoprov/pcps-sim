# Run the targets pipeline for one profile, with its own store.
# Usage: Rscript analysis/run_pipeline.R <profile> <workers>
#   e.g. Rscript analysis/run_pipeline.R pilot_glm 16

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[[1]] else "pilot_glm"
workers <- if (length(args) >= 2) args[[2]] else as.character(parallel::detectCores() - 1)
Sys.setenv(PCPS_PROFILE = profile, PCPS_WORKERS = workers)
store <- paste0("_targets_", profile)

t0 <- Sys.time()
targets::tar_make(store = store, reporter = "balanced")
message(profile, " finished in ", format(Sys.time() - t0))
print(targets::tar_read(rejection, store = store), n = Inf)
