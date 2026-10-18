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

pools <- readr::read_csv(pool_file, show_col_types = FALSE)
# cells whose target richness is out of reach are dropped (decision of
# 2026-10-18: mass effects with phylogenetic competition levels off at ~150
# species whatever the pool size, so N = 200 is infeasible there); a cell
# within 15% of its target is kept, with its realized N recorded per replicate
infeasible <- abs(pools$realized - pools$target) / pools$target > 0.15
if (any(infeasible)) {
  message("Dropped as infeasible: ",
          paste(unique(paste(pools$archetype[infeasible], pools$target[infeasible])),
                collapse = ", "))
  pools <- pools[!infeasible, ]
}
if ("converged" %in% names(pools) && !all(pools$converged))
  warning("pool sizes not converged for ", sum(!pools$converged), " rows")
pools <- dplyr::select(pools, archetype, scenario, target, n_pool)
jobs <- pools |>
  dplyr::inner_join(archetypes_layerB(), by = "archetype") |>
  dplyr::inner_join(scenarios_layerA(), by = "scenario") |>
  tidyr::crossing(rep = seq_len(n_reps))

t0 <- Sys.time()
# a failed replicate (e.g. too few species left) is dropped, not fatal; the
# summary reports the number of replicates kept
one_rep <- purrr::possibly(\(...) {
  job <- list(...)
  t0 <- Sys.time()
  mc <- sim_metacommunity_B(sim_tree_shape(job$n_pool), job$lambda, job$sigma_signal,
                            job$h_signal, job$env_original, dispersal = job$dispersal,
                            competition = job$competition)
  t1 <- Sys.time()
  analyse_metacommunity(mc, c("pfw", "glm"), nperm, extra_axes = TRUE) |>
    dplyr::mutate(archetype = job$archetype, scenario = job$scenario,
                  target = job$target, n_pool = job$n_pool, rep = job$rep,
                  n_realized = ncol(mc$comm), !!!as.list(tree_shape(mc$tree)),
                  sim_seconds = as.numeric(t1 - t0, units = "secs"),
                  test_seconds = as.numeric(Sys.time() - t1, units = "secs"),
                  .before = 1)
})
res <- furrr::future_pmap(jobs, one_rep, .options = opts) |>
  dplyr::bind_rows()
message("Pilot layer B: ", format(Sys.time() - t0))
readr::write_csv(res, file.path(out_dir, sprintf("pilot_layerB_reps%d_perm%d.csv",
                                                 n_reps, nperm)))

rejection_rates(res, archetype, scenario, target, test) |>
  dplyr::select(archetype, scenario, target, test, n, p_max) |>
  tidyr::pivot_wider(names_from = test, values_from = p_max) |>
  as.data.frame() |>
  print(digits = 2)
message("Seconds per replicate (mean): simulation ",
        round(mean(res$sim_seconds[res$test == "adonis"]), 2), ", tests ",
        round(mean(res$test_seconds[res$test == "adonis"]), 2))
