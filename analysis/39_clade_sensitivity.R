# Sensitivity of clade recovery (point 6 of the Codex review, 2026-10-22;
# plan in docs/review/ponto6_plano.md). On the metacommunities of the main run
# (analysis/17_run_clade.R: same trees and seeds), N = 100 and incidence data,
# in three cells (Yule/Yule, low gamma, high Colless), basal and recent clades:
#   - PCPS, PCC and evoPCA with the axis chosen among the first 2, 4 or 8;
#   - phylofactorization with the full objective (phylo + phylo:env) and with
#     the interaction only (phylo:env);
#   - shifts 0, 1 and 2, and a prevalence control (shift 0, maximum abundance of
#     the planted clade x 0.2), where a hit is a false recovery.
# All variants are computed on the same metacommunity, so their contrasts are
# paired by replicate. With 4 axes and the full objective, the shift scenarios
# reproduce the main run exactly, which the script checks.
# Writes data/derived/clade_sensitivity/results.csv,
# data/derived/models/clade_sensitivity.csv (rates with 95% Jeffreys intervals,
# replicates without a score, and paired contrasts against the main-run variant)
# and data/derived/models/clade_sensitivity_ranking.csv (order of the four
# methods for every number of axes and phylofactor objective). In the
# prevalence control, planted species are often lost (h x 0.2), so `n` there is
# smaller than in the shift scenarios.
# Usage: Rscript analysis/39_clade_sensitivity.R [workers] [n_reps]

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[[1]]) else 4L
n_reps <- if (length(args) >= 2) as.integer(args[[2]]) else 1000L
stopifnot(n_reps <= 1000)  # the trees and seeds of the main run

for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
future::plan(future::multisession, workers = workers)
opts <- furrr::furrr_options(seed = NULL, globals = TRUE, scheduling = Inf,
                             packages = c("dplyr", "tibble"))
out_dir <- file.path("data", "derived", "clade_sensitivity")
chunk_dir <- file.path(out_dir, "chunks")
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
z <- stats::qnorm(0.975)

# as in analysis/17_run_clade.R
base_seed <- 20261019L
depths <- c("old", "young")
write_atomic <- function(write, object, path) {
  tmp <- paste0(path, ".tmp")
  write(object, tmp)
  file.rename(tmp, path)
}
grid <- readRDS(file.path("data", "derived", "trees", "tree_grid_k200.rds")) |>
  dplyr::select(n, cell) |>
  dplyr::mutate(tree_id = dplyr::row_number(),
                tree_file = file.path("data", "derived", "clade", "trees",
                                      sprintf("%s_n%d_k1000.rds", cell, n))) |>
  dplyr::filter(n == 100, cell %in% c("g_yule__c_yule", "g_low__c_yule", "g_yule__c_high"))

scenarios <- tibble::tibble(scenario = c("shift0", "shift1", "shift2", "prevalence"),
                            shift = c(0, 1, 2, 0), h_factor = c(1, 1, 1, 0.2))
axes <- c(2, 4, 8)

# node scores of every variant on one metacommunity
variant_scores <- function(mc) {
  comm <- mc$comm
  tree <- mc$tree
  env <- mc$env
  ord <- purrr::map(stats::setNames(axes, axes), \(a) list(
    pcps = node_scores(tree, pcps_species_scores(comm, tree, env, a)),
    pcc = node_scores(tree, pcc_species_scores(comm, tree, env, a)),
    evopca = node_scores(tree, evopca_species_scores(comm, tree, env, a))
  ))
  ord <- unlist(unname(purrr::imap(ord, \(s, a) stats::setNames(s, paste0(names(s), "_", a)))),
                recursive = FALSE)
  c(ord, list(
    phylofactor_full = phylofactor_node_scores(comm, tree, env, "full"),
    phylofactor_interaction = phylofactor_node_scores(comm, tree, env, "interaction")
  ))
}

jobs <- tidyr::crossing(grid, depth = depths, scenarios) |>
  dplyr::mutate(chunk_file = file.path(chunk_dir, sprintf("%s_%s_%s_k%d.csv", cell, depth,
                                                          scenario, n_reps)))
todo <- dplyr::filter(jobs, !file.exists(chunk_file))
message(nrow(todo), " of ", nrow(jobs), " chunks to run")
t0 <- Sys.time()
invisible(furrr::future_pmap(
  dplyr::select(todo, cell, tree_id, tree_file, depth, scenario, shift, h_factor, chunk_file),
  \(cell, tree_id, tree_file, depth, scenario, shift, h_factor, chunk_file) {
    trees <- readRDS(tree_file)
    res <- purrr::map(seq_len(n_reps), \(r) {
      set.seed(base_seed + 1000L * (2L * tree_id + match(depth, depths)) + r)
      mc <- sim_metacommunity_clade(trees[[r]], shift = shift, depth = depth, h_factor = h_factor)
      in_clade <- colnames(mc$comm) %in% mc$planted$tips
      prev <- colMeans(mc$comm > 0)
      dplyr::bind_rows(lapply(variant_scores(mc), recovery_metrics, tree = mc$tree,
                              planted_tips = mc$planted$tips), .id = "variant") |>
        dplyr::mutate(rep = r, prev_in = mean(prev[in_clade]), prev_out = mean(prev[!in_clade]))
    }) |>
      dplyr::bind_rows() |>
      dplyr::mutate(cell = cell, depth = depth, scenario = scenario, .before = 1)
    write_atomic(readr::write_csv, res, chunk_file)
    NULL
  },
  .options = opts
))
message("Chunks: ", format(Sys.time() - t0))

res <- purrr::map(jobs$chunk_file, readr::read_csv, show_col_types = FALSE) |>
  dplyr::bind_rows()
n_variants <- 3 * length(axes) + 2
stopifnot(nrow(res) == n_variants * n_reps * nrow(jobs))
readr::write_csv(res, file.path(out_dir, "results.csv"))

# check: 4 axes and the full objective reproduce the main run
main <- readr::read_csv(file.path("data", "derived", "clade", "results.csv"),
                        show_col_types = FALSE) |>
  dplyr::filter(n_spp == 100, data_type == "incidence", shift %in% c(0, 1, 2)) |>
  dplyr::transmute(cell, depth, scenario = paste0("shift", shift), rep,
                   variant = ifelse(method == "phylofactor", "phylofactor_full",
                                    paste0(method, "_4")),
                   hit_main = hit, rank_main = rank, jaccard_main = jaccard)
check <- dplyr::inner_join(res, main, by = c("cell", "depth", "scenario", "rep", "variant"))
stopifnot(nrow(check) == 4 * n_reps * nrow(dplyr::filter(jobs, scenario != "prevalence")),
          identical(is.na(check$hit), is.na(check$hit_main)),
          all(check$hit == check$hit_main, na.rm = TRUE),
          identical(is.na(check$rank), is.na(check$rank_main)),
          all(check$rank == check$rank_main, na.rm = TRUE),
          isTRUE(all.equal(check$jaccard, check$jaccard_main)))
message("Main run reproduced: ", nrow(check), " rows")

# rates (95% Jeffreys interval) and paired contrasts against the main-run variant
reference <- \(variant) ifelse(grepl("^phylofactor", variant), "phylofactor_full",
                               sub("_[0-9]+$", "_4", variant))
rates <- res |>
  dplyr::group_by(cell, depth, scenario, variant) |>
  dplyr::summarise(n_na = sum(is.na(hit)), k = sum(hit, na.rm = TRUE), n = sum(!is.na(hit)),
                   prev_in = mean(prev_in, na.rm = TRUE), prev_out = mean(prev_out),
                   .groups = "drop") |>
  dplyr::mutate(rate = k / n,
                lower = ifelse(k == 0, 0, stats::qbeta(0.025, k + 0.5, n - k + 0.5)),
                upper = ifelse(k == n, 1, stats::qbeta(0.975, k + 0.5, n - k + 0.5)))
# paired difference of two hit rates: discordant pairs b (only the variant) and c
# (only the reference), with 0.5 added to each so that no discordance does not
# give an interval of zero width
paired <- res |>
  dplyr::mutate(ref = reference(variant)) |>
  dplyr::inner_join(dplyr::select(res, cell, depth, scenario, rep, ref = variant, hit_ref = hit),
                    by = c("cell", "depth", "scenario", "rep", "ref")) |>
  dplyr::filter(!is.na(hit), !is.na(hit_ref)) |>
  dplyr::group_by(cell, depth, scenario, variant) |>
  dplyr::summarise(b = sum(hit & !hit_ref), c = sum(!hit & hit_ref), m = dplyr::n(),
                   .groups = "drop") |>
  dplyr::mutate(diff = (b - c) / m,
                se = sqrt((b + c + 1) - ((b - c)^2) / m) / m,
                diff_lower = diff - z * se, diff_upper = diff + z * se) |>
  dplyr::select(-b, -c, -m, -se)
sensitivity <- dplyr::left_join(rates, paired, by = c("cell", "depth", "scenario", "variant"))
readr::write_csv(sensitivity, file.path("data", "derived", "models", "clade_sensitivity.csv"))

# order of the four methods for every number of axes and phylofactor objective
ranking <- tidyr::crossing(axes = axes, objective = c("full", "interaction")) |>
  purrr::pmap(\(axes, objective) {
    keep <- c(paste0(c("pcps", "pcc", "evopca"), "_", axes), paste0("phylofactor_", objective))
    dplyr::filter(rates, variant %in% keep) |>
      dplyr::mutate(method = sub("_.*$", "", variant)) |>
      dplyr::group_by(cell, depth, scenario) |>
      dplyr::summarise(order = paste(method[order(-rate)], collapse = " > "), .groups = "drop") |>
      dplyr::mutate(axes = axes, objective = objective)
  }) |>
  dplyr::bind_rows()
readr::write_csv(ranking, file.path("data", "derived", "models", "clade_sensitivity_ranking.csv"))
print(dplyr::count(ranking, scenario, order), n = Inf)
print(sensitivity |>
        dplyr::group_by(scenario, variant) |>
        dplyr::summarise(rate = mean(rate), diff = mean(diff), .groups = "drop") |>
        tidyr::pivot_wider(names_from = scenario, values_from = c(rate, diff)),
      n = Inf, width = Inf)
