# Pool sizes of layer B: tree size whose richness after the dynamics is close
# to N = 50/100/200, for each archetype x scenario (decision of 2026-10-17).
# Scenarios that differ only in env_original have the same dynamics, so they
# share one calibration.
# Usage: Rscript analysis/10_calibrate_layerB.R <workers> <reps>

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[[1]]) else 6L
reps <- if (length(args) >= 2) as.integer(args[[2]]) else 10L

for (f in list.files("R", full.names = TRUE)) source(f)
future::plan(future::multisession, workers = workers)
# one job at a time: the cost varies a lot among jobs (target 200, mass
# effects), and static chunks leave most workers idle at the end
opts <- furrr::furrr_options(seed = 20261017L, globals = TRUE, scheduling = Inf,
                             packages = c("dplyr", "tibble"))
# tracked: an input of the layer-B design
dir.create(dirname(layerB_pool_file()), recursive = TRUE, showWarnings = FALSE)

traits <- scenarios_layerA() |>
  dplyr::distinct(lambda, sigma_signal, h_signal)
jobs <- tidyr::crossing(archetypes_layerB(), traits, target = c(50, 100, 200))

t0 <- Sys.time()
cal <- furrr::future_pmap(jobs, \(...) {
  job <- list(...)  # one row: archetype (dispersal, competition) and traits
  calibrate_pool_size(job$target, arch = job, sc = job, reps = reps)
}, .options = opts) |>
  dplyr::bind_rows()
message("Calibration: ", format(Sys.time() - t0))

pool_sizes <- dplyr::bind_cols(jobs, dplyr::select(cal, -target)) |>
  dplyr::inner_join(scenarios_layerA(), by = c("lambda", "sigma_signal", "h_signal"),
                    relationship = "many-to-many") |>
  dplyr::select(archetype, scenario, target, n_pool, realized, iterations, converged,
                dplyr::everything())
readr::write_csv(pool_sizes, layerB_pool_file())
print(as.data.frame(dplyr::distinct(pool_sizes, archetype, lambda, sigma_signal,
                                    h_signal, target, n_pool, realized, iterations,
                                    converged)))
