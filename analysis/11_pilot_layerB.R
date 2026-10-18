# Pilot of layer B: type I error and power of the D_P tests (ADONIS, Mantel,
# GLM on PCPS axes, with Procrustes and RDA) under the process-based dynamics,
# for each archetype x scenario, with the pool sizes calibrated so that the
# realized richness is close to N (analysis/10_calibrate_layerB.R). Yule trees,
# a new one for each replicate; the realized N and the shape of the pruned
# tree are recorded.
# Usage: Rscript analysis/11_pilot_layerB.R <n_reps> <workers> <nperm> [pool_sizes.csv]

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[[1]]) else 50L
workers <- if (length(args) >= 2) as.integer(args[[2]]) else 6L
nperm <- if (length(args) >= 3) as.integer(args[[3]]) else 499L
pool_file <- if (length(args) >= 4) args[[4]] else
  file.path("data", "derived", "layerB", "pool_sizes.csv")

for (f in list.files("R", full.names = TRUE)) source(f)
future::plan(future::multisession, workers = workers)
opts <- furrr::furrr_options(seed = 20261018L, globals = TRUE,
                             packages = c("dplyr", "tibble"))
out_dir <- file.path("data", "derived", "pilot_layerB")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pools <- readr::read_csv(pool_file, show_col_types = FALSE) |>
  dplyr::select(archetype, scenario, target, n_pool)
jobs <- pools |>
  dplyr::inner_join(archetypes_layerB(), by = "archetype") |>
  dplyr::inner_join(scenarios_layerA(), by = "scenario") |>
  tidyr::crossing(rep = seq_len(n_reps))

t0 <- Sys.time()
res <- furrr::future_pmap(jobs, \(...) {
  job <- list(...)
  mc <- sim_metacommunity_B(sim_tree_shape(job$n_pool), job$lambda, job$sigma_signal,
                            job$h_signal, job$env_original, dispersal = job$dispersal,
                            competition = job$competition)
  t1 <- Sys.time()
  analyse_metacommunity(mc, c("pfw", "glm"), nperm, extra_axes = TRUE) |>
    dplyr::mutate(archetype = job$archetype, scenario = job$scenario,
                  target = job$target, n_pool = job$n_pool, rep = job$rep,
                  n_realized = ncol(mc$comm),
                  gamma = treestats::gamma_statistic(mc$tree),
                  colless = treestats::colless(mc$tree, normalization = "yule"),
                  test_seconds = as.numeric(Sys.time() - t1, units = "secs"),
                  .before = 1)
}, .options = opts) |>
  dplyr::bind_rows()
message("Pilot layer B: ", format(Sys.time() - t0))
readr::write_csv(res, file.path(out_dir, "pilot_layerB.csv"))

res |>
  dplyr::group_by(archetype, scenario, target, test) |>
  dplyr::summarise(rej_max = mean(p_max <= 0.05), n_realized = mean(n_realized),
                   .groups = "drop") |>
  tidyr::pivot_wider(names_from = test, values_from = rej_max) |>
  as.data.frame() |>
  print(digits = 2)
