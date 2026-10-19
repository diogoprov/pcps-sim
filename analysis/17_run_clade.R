# Objective 2, full run: recovery of a planted clade by PCPS, PCC, evoPCA and
# phylofactorization (R/clade_recovery.R, R/clade_methods.R) in the 9 grid
# cells x N = 50/100/200, basal or recent clades and shifts 0.5/1/2, with
# incidence data; abundance only for N = 100. Run locally.
#
# Phase 1 draws 1000 trees per cell and N (calibrated (beta, c) of
# data/derived/trees/tree_grid_k200.rds). Phase 2 runs one chunk per cell x N
# x data type x depth x shift; replicate r uses tree r of its cell, so the
# trees are shared across conditions (paired design). Each tree file and chunk
# is written to disk with its own seed and skipped if present, so a run that
# stops can be resumed without changing the results.
# Usage: Rscript analysis/17_run_clade.R <workers> [n_reps] [out_dir]

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[[1]]) else 12L
n_reps <- if (length(args) >= 2) as.integer(args[[2]]) else 1000L
out_dir <- if (length(args) >= 3) args[[3]] else file.path("data", "derived", "clade")

for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
future::plan(future::multisession, workers = workers)
opts <- furrr::furrr_options(seed = NULL, globals = TRUE, scheduling = Inf,
                             packages = c("dplyr", "tibble"))
tree_dir <- file.path(out_dir, "trees")
chunk_dir <- file.path(out_dir, "chunks")
dir.create(tree_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
base_seed <- 20261019L

grid <- readRDS(file.path("data", "derived", "trees", "tree_grid_k200.rds")) |>
  dplyr::select(n, cell, beta, c, gamma_target, colless_target) |>
  dplyr::mutate(tree_id = dplyr::row_number(),
                tree_file = file.path(tree_dir, sprintf("%s_n%d.rds", cell, n)))

# Phase 1: trees -------------------------------------------------------------------
t0 <- Sys.time()
todo <- dplyr::filter(grid, !file.exists(tree_file))
invisible(furrr::future_pmap(todo, \(n, cell, beta, c, gamma_target, colless_target,
                                     tree_id, tree_file) {
  set.seed(base_seed + tree_id)
  tr <- sim_tree_cell(n_reps, n, beta, c, gamma_target, colless_target)
  stopifnot(length(tr) == n_reps)
  saveRDS(tr, tree_file)
}, .options = opts))
message("Trees: ", format(Sys.time() - t0))

# Phase 2: chunks ------------------------------------------------------------------
jobs <- tidyr::crossing(grid, data_type = c("incidence", "abundance"),
                        depth = c("old", "young"), shift = c(0.5, 1, 2)) |>
  dplyr::filter(data_type == "incidence" | n == 100) |>
  dplyr::mutate(chunk_id = dplyr::row_number(),
                chunk_file = file.path(chunk_dir, sprintf("%s_n%d_%s_%s_s%g.csv", cell, n,
                                                          data_type, depth, shift)))
todo <- dplyr::filter(jobs, !file.exists(chunk_file))
message(nrow(todo), " of ", nrow(jobs), " chunks to run")
t0 <- Sys.time()
invisible(furrr::future_pmap(
  dplyr::select(todo, n, cell, tree_file, data_type, depth, shift, chunk_id, chunk_file),
  \(n, cell, tree_file, data_type, depth, shift, chunk_id, chunk_file) {
    set.seed(base_seed + 1000L * chunk_id)
    trees <- readRDS(tree_file)
    incidence <- data_type == "incidence"
    res <- purrr::map(seq_len(n_reps), \(r) {
      t1 <- Sys.time()
      mc <- sim_metacommunity_clade(trees[[r]], shift = shift, depth = depth,
                                    incidence = incidence)
      recover_clade_all(mc, incidence = incidence) |>
        dplyr::mutate(rep = r, rel_age = mc$planted$rel_age,
                      n_planted = sum(mc$tree$tip.label %in% mc$planted$tips),
                      n_kept = ape::Ntip(mc$tree), !!!as.list(tree_shape(mc$tree)),
                      seconds = as.numeric(Sys.time() - t1, units = "secs"))
    }) |>
      dplyr::bind_rows() |>
      dplyr::mutate(cell = cell, n_spp = n, data_type = data_type, depth = depth,
                    shift = shift, .before = 1)
    readr::write_csv(res, chunk_file)
  },
  .options = opts
))
message("Chunks: ", format(Sys.time() - t0))

res <- purrr::map(jobs$chunk_file, readr::read_csv, show_col_types = FALSE) |>
  dplyr::bind_rows()
stopifnot(nrow(res) == 4 * n_reps * nrow(jobs))
readr::write_csv(res, file.path(out_dir, "results.csv"))
message("Written: ", file.path(out_dir, "results.csv"), " (", nrow(res), " rows)")
