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
# - the environment is static and has no spatial structure, and dispersal is
#   uniform among sites (sites remain exchangeable, so the site shuffle stays
#   valid; plan, care 1 of layer B).
# - only the final state is kept (a single survey), not a temporal mean.

#' Competition matrix
#'
#' Intraspecific coefficients 1, interspecific in [0, inter_max), times `scale`
#' (0.05 in Thompson et al. 2020, so that a species alone at its optimum with
#' r = 5 has N* = (5 - 1) / 0.05 = 80). `inter_max` is therefore relative to
#' intraspecific competition. Symmetric. The default 0.1 is weaker than the
#' "stabilising" scenario of Thompson et al. (0.5), which excludes most
#' species of a 100-tip pool (docs/MEMORIA_PROJETO.md).
#' - `none`: no interspecific competition.
#' - `random`: interspecific coefficients ~ U(0, inter_max).
#' - `phylo`: the same U(0, inter_max) values, assigned by rank of phylogenetic
#'   similarity (1 - d / max d), so that close relatives compete more
#'   (competition-relatedness) with the same distribution of coefficients as
#'   `random`, whatever the tree.
competition_matrix <- function(tree, type = c("none", "random", "phylo"),
                               inter_max = 0.1, scale = 0.05) {
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
#' Binomial(N', a_j). The emigrants of each species form a pool that is
#' redistributed multinomially among sites with probabilities proportional to
#' the number of emigrants from the other sites (Thompson et al. 2020 with a
#' uniform kernel), so an individual may return to its own site.
#' @param r sites x species matrix of density-independent growth rates.
#' @param alpha species x species competition matrix ([competition_matrix()]).
#' @param dispersal emigration probability, one value or one per species.
#' @param seed_steps steps at which Poisson(seed_mean) individuals of every
#'   species are added to every site (initial colonisation).
#' @return sites x species abundances at the last step.
mc_dynamics <- function(r, alpha, dispersal, n_steps = 1000,
                        seed_steps = seq(10, 100, 10), seed_mean = 0.5) {
  m <- nrow(r)
  k <- ncol(r)
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
      w <- sweep(-em, 2, tot, "+")
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
#' Traits from [sim_species_traits()], with the maximum growth rate r_max ~
#' U(rmax_range) in the role of h (so `h_signal` gives signal in r_max). E ~
#' U(0, 100) on `n_sites` sites. Output as in layer A ([finish_metacommunity()]).
#' @param dispersal emigration probability (archetypes: species sorting at
#'   intermediate, mass effects at high dispersal).
#' @param competition type of competition matrix ([competition_matrix()]).
sim_metacommunity_B <- function(tree, lambda, sigma_signal, h_signal, env_original,
                                dispersal = 0.01, competition = "none",
                                incidence = TRUE, n_sites = 100, mu_tol = 10,
                                rmax_range = c(2, 8), n_steps = 1000) {
  tr <- sim_species_traits(tree, lambda, sigma_signal, h_signal, mu_tol, rmax_range)
  alpha <- competition_matrix(tree, competition)
  env <- stats::runif(n_sites, 0, 100)
  r <- exp(-sweep(outer(env, tr$u, "-")^2, 2, 2 * tr$sigma^2, "/"))
  r <- sweep(r, 2, tr$h, "*")
  colnames(r) <- tree$tip.label
  comm <- mc_dynamics(r, alpha, dispersal, n_steps = n_steps)
  if (incidence) comm <- (comm > 0) * 1
  finish_metacommunity(comm, tree, env, env_original, tr)
}

#' Layer-B archetypes: dispersal regime x competition
archetypes_layerB <- function() {
  tidyr::crossing(
    tibble::tibble(regime = c("sorting", "mass_effects"), dispersal = c(0.01, 0.5)),
    competition = c("none", "phylo")
  ) |>
    dplyr::mutate(archetype = paste(regime, competition, sep = "__"), .before = 1)
}

#' Mean number of species left after the dynamics, for a pool of `n_pool` tips
#' (Yule trees, one per replicate)
#' @param arch one row of [archetypes_layerB()].
#' @param sc one row of [scenarios_layerA()].
realized_richness <- function(n_pool, arch, sc, reps = 10, n_steps = 1000) {
  mean(vapply(seq_len(reps), \(r) {
    mc <- sim_metacommunity_B(sim_tree_shape(n_pool), sc$lambda, sc$sigma_signal,
                              sc$h_signal, TRUE, dispersal = arch$dispersal,
                              competition = arch$competition, n_steps = n_steps)
    ncol(mc$comm)
  }, 0))
}

#' Pool size whose realized richness after the dynamics is close to `target`
#'
#' Starts at n = target and rescales n by target / realized until the mean
#' realized richness is within `tol` (relative) of the target.
#' @return tibble with the pool size, its mean realized richness and the
#'   number of iterations.
calibrate_pool_size <- function(target, arch, sc, reps = 10, tol = 0.05,
                                max_iter = 5) {
  n <- target
  for (it in seq_len(max_iter)) {
    m <- realized_richness(n, arch, sc, reps)
    if (abs(m - target) / target <= tol) break
    n <- round(n * target / m)
  }
  tibble::tibble(target = target, n_pool = n, realized = m, iterations = it)
}

