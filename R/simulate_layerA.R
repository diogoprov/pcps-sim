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

#' Species traits of layers A and B
#'
#' u_j: BM with Pagel's lambda, rescaled to [-1, 101]. sigma_j ~ N(mu_tol, 10)
#' truncated at > 0; h_j ~ U(h_range); sigma and h either drawn independently or
#' ranked by a BM trait.
#' @param optima "bm" (the BM values rescaled) or "uniform": the BM values are
#'   replaced by ranks by evenly spaced values on [-1, 101], so the marginal
#'   distribution of the optima is the same whatever the shape of the tree and
#'   only their arrangement on the tree differs (reviewer check, 2026-10-20).
#' @return data.frame (u, sigma, h), one row per tip in tip order.
sim_species_traits <- function(tree, lambda, sigma_signal, h_signal, mu_tol = 10,
                               h_range = c(0, 30), optima = c("bm", "uniform")) {
  optima <- match.arg(optima)
  k <- ape::Ntip(tree)
  u <- sim_trait_lambda(tree, lambda)
  u <- if (optima == "uniform") rank_to_marginal(u, seq(-1, 101, length.out = k)) else
    rescale_to(u)
  sigma <- numeric(0)
  while (length(sigma) < k) {
    s <- stats::rnorm(2 * k, mu_tol, 10)
    sigma <- c(sigma, s[s > 0])
  }
  sigma <- sigma[seq_len(k)]
  h <- stats::runif(k, h_range[1], h_range[2])
  if (sigma_signal) sigma <- rank_to_marginal(sim_trait_lambda(tree, 1), sigma)
  if (h_signal) h <- rank_to_marginal(sim_trait_lambda(tree, 1), h)
  data.frame(u = u, sigma = sigma, h = h)
}

#' Metacommunity as analysed, for layers A and B
#'
#' The environment analysed is the one that generated `comm` or, if not
#' `env_original`, a new U(0, 100) draw. Empty sites are dropped, and species
#' absent from all sites are dropped and pruned from the tree (as `PCPS` does
#' and as with real data).
#' @param new_env generator of the independent gradient (spatially structured
#'   in the spatial version of layer B).
#' @return list with comm (sites x species, tip order), tree, env (analysed),
#'   n_dropped_sites, n_dropped_spp, traits of the species kept, and the
#'   analysed gradient at all sites (`env_all`) with the sites kept
#'   (`kept_sites`), used by restricted site permutations.
finish_metacommunity <- function(comm, tree, env, env_original, traits,
                                 new_env = function(n) stats::runif(n, 0, 100)) {
  env_used <- if (env_original) env else new_env(length(env))
  ok <- drop_empty_sites(comm, env_used)
  present <- colSums(ok$comm) > 0
  tree_obs <- if (all(present)) tree else ape::keep.tip(tree, tree$tip.label[present])
  list(comm = ok$comm[, tree_obs$tip.label, drop = FALSE], tree = tree_obs,
       env = ok$env, n_dropped_sites = ok$n_dropped, n_dropped_spp = sum(!present),
       traits = traits[present, ], env_all = env_used, kept_sites = rowSums(comm) > 0)
}

#' One layer-A metacommunity
#'
#' L_ij = h_j exp(-(E_i - u_j)^2 / (2 sigma_j^2)), Poisson counts, with traits
#' from [sim_species_traits()] and E ~ U(0, 100) on `n_sites` sites.
#' @param optima "bm" or "uniform", passed to [sim_species_traits()].
#' @return see [finish_metacommunity()].
sim_metacommunity_A <- function(tree, lambda, sigma_signal, h_signal, env_original,
                                incidence = TRUE, n_sites = 100, mu_tol = 10,
                                optima = "bm") {
  tr <- sim_species_traits(tree, lambda, sigma_signal, h_signal, mu_tol, optima = optima)
  env <- stats::runif(n_sites, 0, 100)
  finish_metacommunity(sim_comm_traits(tr, env, tree$tip.label, incidence), tree, env,
                       env_original, tr)
}

#' Static community from species traits (Gaussian response, Poisson counts)
#'
#' @param traits data.frame with `u`, `sigma` and `h`, one row per species.
#' @param env environmental value of each site.
#' @param species species labels, in the row order of `traits`.
#' @return sites x species matrix of counts, or of 0/1 if `incidence`.
sim_comm_traits <- function(traits, env, species, incidence = TRUE) {
  lambda_ij <- sweep(outer(env, traits$u, "-")^2, 2, 2 * traits$sigma^2, "/")
  lambda_ij <- sweep(exp(-lambda_ij), 2, traits$h, "*")
  comm <- matrix(stats::rpois(length(lambda_ij), lambda_ij), nrow = length(env),
                 dimnames = list(NULL, species))
  if (incidence) comm <- (comm > 0) * 1
  comm
}

#' Run all methods on one metacommunity
#' @param methods subset of c("pfw", "glm", "pcc").
#' @param extra_axes passed to [test_dp_all()] (Procrustes GLM and RDA on the
#'   PCPS axes; needs both "pfw" and "glm").
#' @param site_env optional matrix of permuted gradients (one row per
#'   permutation) that replaces the site shuffle, e.g. torus translations in
#'   the spatial version of layer B ([test_dp_all()], [test_pcc()]).
analyse_metacommunity <- function(mc, methods = c("pfw", "glm", "pcc"),
                                  nperm = 999, nperm_pcc = nperm, extra_axes = FALSE,
                                  site_env = NULL) {
  stopifnot(!extra_axes || all(c("pfw", "glm") %in% methods))
  stopifnot(is.null(site_env) || all(c("pfw", "glm") %in% methods))
  out <- list()
  if (all(c("pfw", "glm") %in% methods)) {
    out$dp <- test_dp_all(mc$comm, mc$tree, mc$env, axes = 1:2, nperm = nperm,
                          extra_axes = extra_axes, site_env = site_env)
  } else if ("pfw" %in% methods) {
    out$pfw <- test_pfw(mc$comm, mc$tree, mc$env, nperm = nperm, sqrt_adonis = TRUE)
  } else if ("glm" %in% methods) {
    out$glm <- test_pcps_glm(mc$comm, mc$tree, mc$env, axes = 1:2, nperm = nperm)
  }
  if ("pcc" %in% methods)
    out$pcc <- test_pcc(mc$comm, mc$tree, mc$env, nperm = nperm_pcc, site_env = site_env)
  dplyr::bind_rows(out) |>
    dplyr::mutate(n_dropped_sites = mc$n_dropped_sites, n_dropped_spp = mc$n_dropped_spp)
}
