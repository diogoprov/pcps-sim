# Pilot of the tree-shape experiment: type I error and power of PFW (Mantel,
# ADONIS on sqrt Bray-Curtis) and PCC in three grid cells, N = 50/100/200.
# A new tree from the cell for each replicate.
# Usage: Rscript analysis/07_pilot_tree_shape.R <n_reps> <workers> <nperm>

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[[1]]) else 200L
workers <- if (length(args) >= 2) as.integer(args[[2]]) else 16L
nperm <- if (length(args) >= 3) as.integer(args[[3]]) else 999L

for (f in list.files("R", full.names = TRUE)) source(f)
future::plan(future::multisession, workers = workers)
opts <- furrr::furrr_options(seed = 20261018L, globals = TRUE,
                             packages = c("dplyr", "tibble"))
out_dir <- file.path("data", "derived", "pilot")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

grid <- readRDS(file.path("data", "derived", "trees", "tree_grid_k200.rds"))
cells <- c("g_yule__c_yule", "g_low__c_yule", "g_yule__c_high")
stopifnot(n_reps <= 200)
trees <- dplyr::filter(grid, cell %in% cells)

jobs <- tidyr::crossing(scenarios_4(), cell = cells, n_spp = c(50, 100, 200),
                        rep = seq_len(n_reps))
tree_of <- \(cell, n_spp, rep)
  trees$trees[[which(trees$cell == cell & trees$n == n_spp)]][[rep]]

t0 <- Sys.time()
res <- furrr::future_pmap(
  jobs,
  \(scenario, phylo, env, cell, n_spp, rep) {
    tr <- tree_of(cell, n_spp, rep)
    rep_pilot(tr, phylo, env, nperm = nperm) |>
      dplyr::mutate(scenario = scenario, cell = cell, n_spp = n_spp, rep = rep,
                    gamma = treestats::gamma_statistic(tr),
                    colless = treestats::colless(tr, normalization = "yule"))
  },
  .options = opts
) |> dplyr::bind_rows()
message("Pilot: ", format(Sys.time() - t0))
readr::write_csv(res, file.path(out_dir, sprintf("pilot_reps%d_perm%d.csv", n_reps, nperm)))

rejection_rates(res, test, scenario, cell, n_spp) |>
  dplyr::arrange(test, scenario, n_spp, cell) |>
  print(n = Inf)
