# Objective 2, full run: recovery of a planted clade by PCPS, PCC, evoPCA and
# phylofactorization (R/clade_recovery.R, R/clade_methods.R) in the 9 grid
# cells x N = 50/100/200, basal or recent clades and shifts 0.5/1/2, with
# incidence data; abundance only for N = 100. Run locally (Diogo, 2026-10-19:
# ~2 h on 12 cores, cheaper than rebuilding the AWS image with Bioconductor).
#
# Phase 1 draws `n_reps` new trees per cell and N with the calibrated
# (beta, c) of data/derived/trees/tree_grid_k200.rds (the cells of layer A).
# Phase 2 runs one chunk per cell x N x data type x depth x shift. Paired
# design: replicate r uses tree r and a seed set by (tree, depth, r), so for a
# given replicate the planted clade, optima and gradient are the same across
# shifts and data types, and only the factor of interest changes. Files carry
# `n_reps` in their names and are written atomically, so a run that stops can
# be resumed without changing the results.
# Usage: Rscript analysis/17_run_clade.R <workers> [n_reps] [out_dir]

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[[1]]) else 12L
n_reps <- if (length(args) >= 2) as.integer(args[[2]]) else 1000L
out_dir <- if (length(args) >= 3) args[[3]] else file.path("data", "derived", "clade")
stopifnot(n_reps <= 1000)  # the replicate seeds below assume it

for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
future::plan(future::multisession, workers = workers)
opts <- furrr::furrr_options(seed = NULL, globals = TRUE, scheduling = Inf,
                             packages = c("dplyr", "tibble"))
tree_dir <- file.path(out_dir, "trees")
chunk_dir <- file.path(out_dir, "chunks")
dir.create(tree_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
base_seed <- 20261019L
depths <- c("old", "young")

write_atomic <- function(write, object, path) {
  tmp <- paste0(path, ".tmp")
  write(object, tmp)
  file.rename(tmp, path)
}

grid <- readRDS(file.path("data", "derived", "trees", "tree_grid_k200.rds")) |>
  dplyr::select(n, cell, beta, c, gamma_target, colless_target) |>
  dplyr::mutate(tree_id = dplyr::row_number(),
                tree_file = file.path(tree_dir, sprintf("%s_n%d_k%d.rds", cell, n, n_reps)))

# Phase 1: trees -------------------------------------------------------------------
t0 <- Sys.time()
todo <- dplyr::filter(grid, !file.exists(tree_file))
invisible(furrr::future_pmap(todo, \(n, cell, beta, c, gamma_target, colless_target,
                                     tree_id, tree_file) {
  set.seed(base_seed + tree_id)
  tr <- sim_tree_cell(n_reps, n, beta, c, gamma_target, colless_target)
  stopifnot(length(tr) == n_reps)
  write_atomic(saveRDS, tr, tree_file)
  NULL
}, .options = opts))
message("Trees: ", format(Sys.time() - t0))

# Phase 2: chunks ------------------------------------------------------------------
jobs <- tidyr::crossing(grid, data_type = c("incidence", "abundance"), depth = depths,
                        shift = c(0.5, 1, 2)) |>
  dplyr::filter(data_type == "incidence" | n == 100) |>
  dplyr::mutate(chunk_file = file.path(chunk_dir, sprintf("%s_n%d_%s_%s_s%g_k%d.csv", cell,
                                                          n, data_type, depth, shift,
                                                          n_reps)))
todo <- dplyr::filter(jobs, !file.exists(chunk_file))
message(nrow(todo), " of ", nrow(jobs), " chunks to run")
t0 <- Sys.time()
invisible(furrr::future_pmap(
  dplyr::select(todo, n, cell, tree_id, tree_file, data_type, depth, shift, chunk_file),
  \(n, cell, tree_id, tree_file, data_type, depth, shift, chunk_file) {
    trees <- readRDS(tree_file)
    incidence <- data_type == "incidence"
    res <- purrr::map(seq_len(n_reps), \(r) {
      set.seed(base_seed + 1000L * (2L * tree_id + match(depth, depths)) + r)
      run_clade_rep(trees[[r]], shift, depth, incidence) |>
        dplyr::mutate(rep = r)
    }) |>
      dplyr::bind_rows() |>
      dplyr::mutate(cell = cell, n_spp = n, data_type = data_type, depth = depth,
                    shift = shift, .before = 1)
    write_atomic(readr::write_csv, res, chunk_file)
    NULL
  },
  .options = opts
))
message("Chunks: ", format(Sys.time() - t0))

res <- purrr::map(jobs$chunk_file, readr::read_csv, show_col_types = FALSE) |>
  dplyr::bind_rows()
stopifnot(nrow(res) == 4 * n_reps * nrow(jobs))
readr::write_csv(res, file.path(out_dir, "results.csv"))
message("Written: ", file.path(out_dir, "results.csv"), " (", nrow(res), " rows)")
