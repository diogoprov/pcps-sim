# Replicates of the published simulation designs, used to validate the pipeline.

#' Scenario definitions shared by both papers: does the phylogeny structure the
#' optima, and is the original gradient (TRUE) or a random one (FALSE) analysed?
scenarios_4 <- function() {
  tibble::tribble(
    ~scenario,    ~phylo, ~env,
    "neither",    FALSE,  FALSE,
    "phylo_only", TRUE,   FALSE,
    "env_only",   FALSE,  TRUE,
    "both",       TRUE,   TRUE
  )
}

#' One replicate of Duarte et al. 2016 (Table 1): PFW, Mantel and ADONIS
#'
#' New Yule tree per replicate; Grafen's rho = 1 (signal) or 1e-4 (none);
#' incidence data; 100 sites; mu_tol = 10.
rep_duarte2016 <- function(n_spp, phylo, env_original, mu_tol = 10,
                           n_sites = 100, nperm = 999, sqrt_adonis = FALSE) {
  tree <- sim_tree_duarte(n_spp)
  u <- sim_optima_grafen(tree, rho = if (phylo) 1 else 1e-4)
  env <- stats::runif(n_sites, 0, 100)
  comm <- sim_comm_static(u, env, mu_tol = mu_tol, incidence = TRUE)
  env_used <- if (env_original) env else stats::runif(n_sites, 0, 100)
  ok <- drop_empty_sites(comm, env_used)
  test_pfw(ok$comm, tree, ok$env, nperm = nperm, sqrt_adonis = sqrt_adonis) |>
    dplyr::mutate(n_dropped = ok$n_dropped)
}

#' One replicate of Peres-Neto & Kembel 2015 (Table 1): PCC, one gradient
#'
#' Same tree in all replicates (`tree` argument); 100 sites x 50 species;
#' incidence data; 99 permutations, as in the paper.
rep_pnk2015 <- function(tree, phylo, env_original, mu_tol,
                        n_sites = 100, nperm = 99, squared = FALSE) {
  u <- sim_optima_pnk(tree, signal = phylo)
  env <- stats::runif(n_sites, 0, 100)
  comm <- sim_comm_static(u, env, mu_tol = mu_tol, incidence = TRUE)
  env_used <- if (env_original) env else stats::runif(n_sites, 0, 100)
  ok <- drop_empty_sites(comm, env_used)
  test_pcc(ok$comm, tree, ok$env, nperm = nperm, squared = squared) |>
    dplyr::mutate(n_dropped = ok$n_dropped)
}

#' Rejection rates at alpha with binomial 95% CI
rejection_rates <- function(res, ..., alpha = 0.05) {
  res |>
    dplyr::group_by(...) |>
    dplyr::summarise(
      n = dplyr::n(),
      dplyr::across(c(p_site, p_taxa, p_max), \(p) mean(p <= alpha)),
      .groups = "drop"
    )
}

#' One replicate of the tree-shape pilot: PFW (Mantel, ADONIS) and PCC on the
#' same metacommunity, for a tree drawn from a grid cell
#'
#' Optima evolve by Brownian motion on the actual branch lengths
#' (`sim_optima_pnk`), not on Grafen's lengths, which would erase gamma.
rep_pilot <- function(tree, phylo, env_original, mu_tol = 10, n_sites = 100,
                      nperm = 999) {
  u <- sim_optima_pnk(tree, signal = phylo)
  env <- stats::runif(n_sites, 0, 100)
  comm <- sim_comm_static(u, env, mu_tol = mu_tol, incidence = TRUE)
  env_used <- if (env_original) env else stats::runif(n_sites, 0, 100)
  ok <- drop_empty_sites(comm, env_used)
  pfw <- test_pfw(ok$comm, tree, ok$env, nperm = nperm, sqrt_adonis = TRUE)
  pcc <- test_pcc(ok$comm, tree, ok$env, nperm = nperm)
  dplyr::bind_rows(pfw, pcc) |>
    dplyr::mutate(n_dropped = ok$n_dropped)
}
