# Peres-Neto & Kembel 2015 used a single Yule tree, so their power and their
# phylo_only column-test rates depend on it. Spread of those rates across trees.
# Usage: Rscript analysis/03_pnk_tree_spread.R <n_trees> <n_reps> <workers>

args <- commandArgs(trailingOnly = TRUE)
n_trees <- if (length(args) >= 1) as.integer(args[[1]]) else 20L
n_reps <- if (length(args) >= 2) as.integer(args[[2]]) else 200L
workers <- if (length(args) >= 3) as.integer(args[[3]]) else 16L

for (f in list.files("R", full.names = TRUE)) source(f)
future::plan(future::multisession, workers = workers)
opts <- furrr::furrr_options(seed = 20261017L, globals = TRUE,
                             packages = c("dplyr", "tibble"))
out_dir <- file.path("data", "derived", "validation")

set.seed(2015)
trees <- replicate(n_trees, ape::rphylo(50, birth = 1, death = 0), simplify = FALSE)
grid <- tidyr::crossing(
  dplyr::filter(scenarios_4(), scenario %in% c("both", "phylo_only")),
  tree_id = seq_len(n_trees), mu_tol = c(5, 10), squared = c(FALSE, TRUE),
  rep = seq_len(n_reps)
)
t0 <- Sys.time()
res <- furrr::future_pmap(
  grid,
  \(scenario, phylo, env, tree_id, mu_tol, squared, rep)
    rep_pnk2015(trees[[tree_id]], phylo, env, mu_tol = mu_tol, squared = squared) |>
    dplyr::mutate(scenario = scenario, tree_id = tree_id, mu_tol = mu_tol,
                  squared = squared, rep = rep),
  .options = opts
) |> dplyr::bind_rows()
message("PNK tree spread: ", format(Sys.time() - t0))
readr::write_csv(res, file.path(out_dir, sprintf("pnk2015_trees%d_reps%d.csv", n_trees, n_reps)))

per_tree <- rejection_rates(res, squared, scenario, mu_tol, tree_id)
per_tree |>
  dplyr::summarise(
    dplyr::across(c(p_taxa, p_max), list(
      mean = mean, min = min, q10 = \(x) quantile(x, 0.1), q90 = \(x) quantile(x, 0.9),
      max = max)),
    .by = c(squared, scenario, mu_tol)
  ) |>
  print(width = Inf)
