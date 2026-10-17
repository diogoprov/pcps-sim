# Layer A of the new study: static metacommunities (Peres-Neto, Leibold & Dray
# 2012 protocol) on trees from the shape grid, with the eight scenarios of
# docs/PLANO_SIMULACAO.md, section 3.

#' Scenario table (docs/PLANO_SIMULACAO.md, section 3)
#'
#' lambda: Pagel's lambda of the optima; sigma_signal / h_signal: niche breadth
#' or maximum abundance evolve by BM (lambda = 1); env_original: analyse the
#' gradient that generated the data (TRUE) or an independent one (FALSE).
scenarios_layerA <- function() {
  tibble::tribble(
    ~scenario, ~lambda, ~sigma_signal, ~h_signal, ~env_original, ~estimates,
    "S0",      0,       FALSE,         FALSE,     FALSE,         "type_I",
    "S1",      1,       FALSE,         FALSE,     FALSE,         "type_I",
    "S2",      0,       FALSE,         FALSE,     TRUE,          "type_I_phylo",
    "S2sigma", 0,       TRUE,          FALSE,     TRUE,          "type_I_phylo",
    "S2h",     0,       FALSE,         TRUE,      TRUE,          "type_I_phylo",
    "S3a",     0.5,     FALSE,         FALSE,     TRUE,          "power",
    "S3b",     1,       FALSE,         FALSE,     TRUE,          "power",
    "S3sigma", 1,       TRUE,          FALSE,     TRUE,          "power"
  )
}

#' Trait on the tree under Brownian motion with Pagel's lambda
#'
#' Multivariate normal with covariance lambda * C + (1 - lambda) * diag(C),
#' C = vcv(tree); lambda = 0 gives independent species.
sim_trait_lambda <- function(tree, lambda = 1) {
  cv <- ape::vcv.phylo(tree)
  cv <- lambda * cv + (1 - lambda) * diag(diag(cv))
  x <- drop(crossprod(chol(cv), stats::rnorm(ape::Ntip(tree))))
  stats::setNames(x, tree$tip.label)
}

#' Impose a target marginal distribution on a trait by ranks
#'
#' Keeps the phylogenetic ordering of `x` and replaces its values by a sorted
#' sample from the target, so that scenarios with and without signal share the
#' same marginal distribution of sigma or h.
rank_to_marginal <- function(x, target) {
  stats::setNames(sort(target)[rank(x, ties.method = "random")], names(x))
}

#' One layer-A metacommunity
#'
#' L_ij = h_j exp(-(E_i - u_j)^2 / (2 sigma_j^2)), Poisson counts.
#' u_j: BM with Pagel's lambda, rescaled to [-1, 101]. sigma_j ~ N(mu_tol, 10)
#' truncated at > 0; h_j ~ U(0, 30); either drawn independently or ranked by a BM
#' trait. Empty sites are dropped, and species absent from all sites are dropped
#' and pruned from the tree (as `PCPS` does and as with real data).
#' @return list with comm (sites x species, tip order), tree, env (analysed),
#'   n_dropped_sites, n_dropped_spp.
sim_metacommunity_A <- function(tree, lambda, sigma_signal, h_signal, env_original,
                                incidence = TRUE, n_sites = 100, mu_tol = 10) {
  k <- ape::Ntip(tree)
  u <- rescale_to(sim_trait_lambda(tree, lambda))
  sigma <- numeric(0)
  while (length(sigma) < k) {
    s <- stats::rnorm(2 * k, mu_tol, 10)
    sigma <- c(sigma, s[s > 0])
  }
  sigma <- sigma[seq_len(k)]
  h <- stats::runif(k, 0, 30)
  if (sigma_signal) sigma <- rank_to_marginal(sim_trait_lambda(tree, 1), sigma)
  if (h_signal) h <- rank_to_marginal(sim_trait_lambda(tree, 1), h)

  env <- stats::runif(n_sites, 0, 100)
  lambda_ij <- sweep(outer(env, u, "-")^2, 2, 2 * sigma^2, "/")
  lambda_ij <- sweep(exp(-lambda_ij), 2, h, "*")
  comm <- matrix(stats::rpois(length(lambda_ij), lambda_ij), nrow = n_sites,
                 dimnames = list(NULL, tree$tip.label))
  if (incidence) comm <- (comm > 0) * 1

  env_used <- if (env_original) env else stats::runif(n_sites, 0, 100)
  ok <- drop_empty_sites(comm, env_used)
  present <- colSums(ok$comm) > 0
  tree_obs <- if (all(present)) tree else ape::keep.tip(tree, tree$tip.label[present])
  list(comm = ok$comm[, tree_obs$tip.label, drop = FALSE], tree = tree_obs,
       env = ok$env, n_dropped_sites = ok$n_dropped, n_dropped_spp = sum(!present))
}

#' Run all methods on one metacommunity
#' @param methods subset of c("pfw", "glm", "pcc").
analyse_metacommunity <- function(mc, methods = c("pfw", "glm", "pcc"),
                                  nperm = 999, nperm_pcc = nperm) {
  out <- list()
  if (all(c("pfw", "glm") %in% methods)) {
    out$dp <- test_dp_all(mc$comm, mc$tree, mc$env, axes = 1:2, nperm = nperm)
  } else if ("pfw" %in% methods) {
    out$pfw <- test_pfw(mc$comm, mc$tree, mc$env, nperm = nperm, sqrt_adonis = TRUE)
  } else if ("glm" %in% methods) {
    out$glm <- test_pcps_glm(mc$comm, mc$tree, mc$env, axes = 1:2, nperm = nperm)
  }
  if ("pcc" %in% methods)
    out$pcc <- test_pcc(mc$comm, mc$tree, mc$env, nperm = nperm_pcc)
  dplyr::bind_rows(out) |>
    dplyr::mutate(n_dropped_sites = mc$n_dropped_sites, n_dropped_spp = mc$n_dropped_spp)
}
