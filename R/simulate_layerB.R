# Layer B: process-based metacommunity dynamics (docs/PLANO_SIMULACAO.md, 3b),
# after the model of Thompson et al. 2020 (Ecol. Lett. 23: 1314-1329,
# doi:10.1111/ele.13568): density-independent Gaussian growth along the
# environment, Beverton-Holt competition, Poisson demographic stochasticity and
# binomial emigration. Written from the published equations, not from the
# authors' code (neither the Julia code nor mcomsimr has a licence).
#
# Differences from the original model, all deliberate:
# - sigma is the SD of the Gaussian response, exp(-(E - u)^2 / (2 sigma^2)), as
#   in layer A. Thompson et al. write exp(-((u - E) / (2 sigma_T))^2), so
#   sigma = sqrt(2) * sigma_T.
# - alpha[j, k] is the effect of species k on species j (the usual
#   convention); the original code uses the transpose.
# - the environment is static and has no spatial structure, and by default
#   dispersal is uniform among sites (sites remain exchangeable, so the site
#   shuffle stays valid; plan, care 1 of layer B).
# - only the final state is kept (a single survey), not a temporal mean.

#' Competition matrix
#'
#' Intraspecific coefficients 1, interspecific in [0, inter_max), times `scale`
#' (0.05 in Thompson et al. 2020, so that a species alone at its optimum with
#' r = 5 has N* = (5 - 1) / 0.05 = 80). Symmetric.
#' - `none`: no interspecific competition.
#' - `random`: interspecific coefficients ~ U(0, inter_max).
#' - `phylo`: the same U(0, inter_max) values, assigned by rank of phylogenetic
#'   similarity (1 - d / max d), so that close relatives compete more
#'   (competition-relatedness) with the same distribution of coefficients as
#'   `random`, whatever the tree.
competition_matrix <- function(tree, type = c("none", "random", "phylo"),
                               inter_max = 0.5, scale = 0.05) {
  type <- match.arg(type)
  k <- ape::Ntip(tree)
  a <- matrix(0, k, k, dimnames = list(tree$tip.label, tree$tip.label))
  up <- upper.tri(a)
  vals <- if (type == "none") numeric(sum(up)) else stats::runif(sum(up), 0, inter_max)
  if (type == "phylo") {
    d <- ape::cophenetic.phylo(tree)[tree$tip.label, tree$tip.label]
    s <- (1 - d / max(d))[up]
    vals <- sort(vals)[rank(s, ties.method = "random")]
  }
  a[up] <- vals
  a <- a + t(a)
  diag(a) <- 1
  scale * a
}

#' Metacommunity dynamics
#'
#' Each step: N' ~ Poisson(r N / (1 + N alpha')), then emigrants ~
#' Binomial(N', a_j); the emigrants of each species are redistributed
#' multinomially among sites with probabilities proportional to
#' disp_mat %*% emigrants (Thompson et al. 2020). With the default uniform
#' kernel an emigrant goes to any other site with equal probability.
#' @param r sites x species matrix of density-independent growth rates.
#' @param alpha species x species competition matrix ([competition_matrix()]).
#' @param dispersal emigration probability, one value or one per species.
#' @param disp_mat sites x sites dispersal kernel (columns = origin); NULL for
#'   uniform dispersal among the other sites.
#' @param seed_steps steps at which Poisson(seed_mean) individuals of every
#'   species are added to every site (initial colonisation).
#' @return sites x species abundances at the last step.
mc_dynamics <- function(r, alpha, dispersal, n_steps = 1000, disp_mat = NULL,
                        seed_steps = seq(10, 100, 10), seed_mean = 0.5) {
  m <- nrow(r)
  k <- ncol(r)
  if (is.null(disp_mat)) disp_mat <- (1 - diag(m)) / (m - 1)
  p_em <- rep(rep_len(dispersal, k), each = m)
  at <- t(alpha)
  n <- matrix(stats::rpois(m * k, seed_mean), m, k)
  for (step in seq_len(n_steps)) {
    if (step %in% seed_steps) n <- n + stats::rpois(m * k, seed_mean)
    lam <- r * n / (1 + n %*% at)
    n <- matrix(stats::rpois(m * k, lam), m, k)
    em <- matrix(stats::rbinom(m * k, n, p_em), m, k)
    tot <- colSums(em)
    if (any(tot > 0)) {
      w <- disp_mat %*% em
      imm <- matrix(0L, m, k)
      for (j in which(tot > 0)) imm[, j] <- stats::rmultinom(1, tot[j], w[, j])
      n <- n - em + imm
    }
  }
  dimnames(n) <- dimnames(r)
  n
}

#' One layer-B metacommunity
#'
#' Traits as in layer A ([sim_metacommunity_A()]): optima u by BM with Pagel's
#' lambda rescaled to [-1, 101]; sigma ~ N(mu_tol, 10) truncated at > 0; the
#' maximum growth rate r_max ~ U(rmax_range) takes the role of h. Sigma and
#' r_max are drawn independently or ranked by a BM trait (signal). E ~ U(0, 100)
#' on `n_sites` sites. Same treatment of empty sites and absent species as in
#' layer A, and the same output.
#' @param dispersal emigration probability (archetypes: species sorting at
#'   intermediate, mass effects at high dispersal).
#' @param competition type of competition matrix ([competition_matrix()]).
sim_metacommunity_B <- function(tree, lambda, sigma_signal, h_signal, env_original,
                                dispersal = 0.01, competition = "none",
                                incidence = TRUE, n_sites = 100, mu_tol = 10,
                                rmax_range = c(2, 8), n_steps = 1000, ...) {
  k <- ape::Ntip(tree)
  u <- rescale_to(sim_trait_lambda(tree, lambda))
  sigma <- numeric(0)
  while (length(sigma) < k) {
    s <- stats::rnorm(2 * k, mu_tol, 10)
    sigma <- c(sigma, s[s > 0])
  }
  sigma <- sigma[seq_len(k)]
  rmax <- stats::runif(k, rmax_range[1], rmax_range[2])
  if (sigma_signal) sigma <- rank_to_marginal(sim_trait_lambda(tree, 1), sigma)
  if (h_signal) rmax <- rank_to_marginal(sim_trait_lambda(tree, 1), rmax)
  alpha <- competition_matrix(tree, competition)

  env <- stats::runif(n_sites, 0, 100)
  r <- exp(-sweep(outer(env, u, "-")^2, 2, 2 * sigma^2, "/"))
  r <- sweep(r, 2, rmax, "*")
  colnames(r) <- tree$tip.label
  comm <- mc_dynamics(r, alpha, dispersal, n_steps = n_steps, ...)
  if (incidence) comm <- (comm > 0) * 1

  env_used <- if (env_original) env else stats::runif(n_sites, 0, 100)
  ok <- drop_empty_sites(comm, env_used)
  present <- colSums(ok$comm) > 0
  tree_obs <- if (all(present)) tree else ape::keep.tip(tree, tree$tip.label[present])
  list(comm = ok$comm[, tree_obs$tip.label, drop = FALSE], tree = tree_obs,
       env = ok$env, n_dropped_sites = ok$n_dropped, n_dropped_spp = sum(!present),
       traits = data.frame(u = u, sigma = sigma, h = rmax)[present, ])
}
