# Pilot of objective 2: recovery of a planted clade by PCPS, PCC, evoPCA and
# binomial phylofactorization (R/clade_recovery.R, R/clade_methods.R), in
# three grid cells, N = 50/100/200, basal or recent clades and two shifts,
# incidence data. A new tree from the cell for each replicate. Measures the
# time per replicate to size the full run.
# Usage: Rscript analysis/16_pilot_clade.R <n_reps> <workers>

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[[1]]) else 50L
workers <- if (length(args) >= 2) as.integer(args[[2]]) else 12L

for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
future::plan(future::multisession, workers = workers)
opts <- furrr::furrr_options(seed = 20261019L, globals = TRUE, scheduling = Inf,
                             packages = c("dplyr", "tibble"))
out_dir <- file.path("data", "derived", "pilot_clade")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

grid <- readRDS(file.path("data", "derived", "trees", "tree_grid_k200.rds"))
cells <- c("g_yule__c_yule", "g_low__c_yule", "g_yule__c_high")
trees <- dplyr::filter(grid, cell %in% cells)
stopifnot(n_reps <= min(lengths(trees$trees)))

jobs <- tidyr::crossing(cell = cells, n_spp = c(50, 100, 200), depth = c("old", "young"),
                        shift = c(1, 2), rep = seq_len(n_reps))
tree_of <- \(cell, n_spp, rep)
  trees$trees[[which(trees$cell == cell & trees$n == n_spp)]][[rep]]

t0 <- Sys.time()
res <- furrr::future_pmap(
  jobs,
  \(cell, n_spp, depth, shift, rep) {
    run_clade_rep(tree_of(cell, n_spp, rep), shift, depth) |>
      dplyr::mutate(cell = cell, n_spp = n_spp, depth = depth, shift = shift, rep = rep)
  },
  .options = opts
) |> dplyr::bind_rows()
message("Pilot: ", format(Sys.time() - t0))
readr::write_csv(res, file.path(out_dir, sprintf("pilot_clade_reps%d.csv", n_reps)))

res |>
  dplyr::group_by(method, n_spp, depth, shift) |>
  dplyr::summarise(hit = mean(hit, na.rm = TRUE), jaccard = mean(jaccard, na.rm = TRUE),
                   rank = stats::median(rank, na.rm = TRUE), .groups = "drop") |>
  print(n = Inf)
res |>
  dplyr::distinct(cell, n_spp, depth, shift, rep, seconds) |>
  dplyr::group_by(n_spp) |>
  dplyr::summarise(mean_s = mean(seconds), max_s = max(seconds)) |>
  print()
