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
#'
#' Species absent from every site after the last seeding step can never come
#' back, so from then on they are dropped from the computations (the product
#' with alpha dominates the cost and is quadratic in the number of species).
#' @param r sites x species matrix of density-independent growth rates.
#' @param alpha species x species competition matrix ([competition_matrix()]).
#' @param dispersal emigration probability, one value or one per species.
#' @param seed_steps steps at which Poisson(seed_mean) individuals of every
#'   species are added to every site (initial colonisation).
#' @return sites x species abundances at the last step.
mc_dynamics <- function(r, alpha, dispersal, n_steps = 1000,
                        seed_steps = seq(10, 100, 10), seed_mean = 0.5,
                        disp_matrix = NULL) {
  m <- nrow(r)
  k <- ncol(r)
  dispersal <- rep_len(dispersal, k)
  alive <- seq_len(k)
  r_alive <- r
  p_em <- rep(dispersal, each = m)
  at <- t(alpha)
  n <- matrix(stats::rpois(m * k, seed_mean), m, k)
  last_seed <- max(c(0, seed_steps))
  for (step in seq_len(n_steps)) {
    if (step %in% seed_steps) n <- n + stats::rpois(m * k, seed_mean)
    keep <- colSums(n) > 0
    if (step > last_seed && !all(keep)) {
      alive <- alive[keep]
      n <- n[, keep, drop = FALSE]
      r_alive <- r_alive[, keep, drop = FALSE]
      at <- at[keep, keep, drop = FALSE]
      k <- length(alive)
      p_em <- rep(dispersal[alive], each = m)
    }
    lam <- r_alive * n / (1 + n %*% at)
    n <- matrix(stats::rpois(m * k, lam), m, k)
    em <- matrix(stats::rbinom(m * k, n, p_em), m, k)
    tot <- colSums(em)
    if (!is.null(disp_matrix)) {
      # spatial version: emigrants of site s reach site t with probability
      # disp_matrix[s, t]; immigrants drawn as Poisson around the expected flow
      if (any(tot > 0))
        n <- n - em + matrix(stats::rpois(m * k, crossprod(disp_matrix, em)), m, k)
    } else if (any(tot > 0)) {
      w <- sweep(-em, 2, tot, "+")
      imm <- matrix(0L, m, k)
      for (j in which(tot > 0)) imm[, j] <- stats::rmultinom(1, tot[j], w[, j])
      n <- n - em + imm
    }
  }
  out <- matrix(0L, m, ncol(r), dimnames = dimnames(r))
  out[, alive] <- n
  out
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
                                rmax_range = c(2, 8), n_steps = 1000, space = NULL) {
  tr <- sim_species_traits(tree, lambda, sigma_signal, h_signal, mu_tol, rmax_range)
  alpha <- competition_matrix(tree, competition)
  if (!is.null(space)) n_sites <- nrow(space$dist)
  new_env <- if (is.null(space)) function(n) stats::runif(n, 0, 100) else
    function(n) spatial_env(space)
  env <- new_env(n_sites)
  r <- exp(-sweep(outer(env, tr$u, "-")^2, 2, 2 * tr$sigma^2, "/"))
  r <- sweep(r, 2, tr$h, "*")
  colnames(r) <- tree$tip.label
  comm <- mc_dynamics(r, alpha, dispersal, n_steps = n_steps,
                      disp_matrix = if (is.null(space)) NULL else space$disp)
  if (incidence) comm <- (comm > 0) * 1
  finish_metacommunity(comm, tree, env, env_original, tr, new_env = new_env)
}

#' Spatial landscape for layer B (reviewer check, 2026-10-20)
#'
#' Sites on a `side` x `side` grid wrapped as a torus (no edges). The gradient
#' is spatially autocorrelated: independent normal deviates smoothed by a
#' Gaussian kernel of `env_bandwidth` cells and given the U(0, 100) marginal
#' of the non-spatial version by ranks. Emigrants reach the other sites with
#' probability proportional to exp(-d / `disp_scale`) (d in cells, toroidal).
#' Site permutations that keep the spatial autocorrelation of the gradient are
#' the torus translations and the eight symmetries of the square
#' (`site_perms`, side^2 x 8 - 1 non-identity rearrangements).
#' @return list with dist (toroidal distances), kernel (smoothing weights),
#'   disp (row-stochastic dispersal matrix, zero diagonal) and site_perms
#'   (matrix, one permutation of the site indices per row).
spatial_landscape <- function(side = 10, env_bandwidth = 2, disp_scale = 2) {
  xy <- expand.grid(x = seq_len(side) - 1, y = seq_len(side) - 1)
  dx <- abs(outer(xy$x, xy$x, "-"))
  dy <- abs(outer(xy$y, xy$y, "-"))
  d <- sqrt(pmin(dx, side - dx)^2 + pmin(dy, side - dy)^2)
  disp <- exp(-d / disp_scale)
  diag(disp) <- 0
  disp <- disp / rowSums(disp)
  maps <- list(\(x, y) cbind(x, y), \(x, y) cbind(y, x), \(x, y) cbind(-x, y),
               \(x, y) cbind(x, -y), \(x, y) cbind(-x, -y), \(x, y) cbind(-y, x),
               \(x, y) cbind(y, -x), \(x, y) cbind(-y, -x))
  perms <- list()
  for (f in maps) for (sx in seq_len(side) - 1) for (sy in seq_len(side) - 1) {
    p <- f(xy$x, xy$y)
    p <- cbind((p[, 1] + sx) %% side, (p[, 2] + sy) %% side)
    perms[[length(perms) + 1]] <- p[, 1] + side * p[, 2] + 1
  }
  perms <- unique(do.call(rbind, perms))
  perms <- perms[!apply(perms, 1, \(p) all(p == seq_len(side^2))), , drop = FALSE]
  list(dist = d, kernel = exp(-d^2 / (2 * env_bandwidth^2)), disp = disp,
       site_perms = perms)
}

#' Spatially autocorrelated gradient on a [spatial_landscape()]
spatial_env <- function(space) {
  z <- drop(space$kernel %*% stats::rnorm(nrow(space$kernel)))
  rank_to_marginal(z, stats::runif(length(z), 0, 100))
}

#' Gradients permuted by the torus rearrangements, at the sites kept
#'
#' Row b holds the gradient at the kept sites after rearrangement b of all
#' sites, for the `site_env` argument of [analyse_metacommunity()].
torus_site_env <- function(mc, space) {
  t(apply(space$site_perms, 1, \(p) mc$env_all[p][mc$kept_sites]))
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
#'
#' Yule trees, one per replicate: the calibration assumes that realized
#' richness depends little on tree shape (to be checked in the grid cells,
#' where the realized N of each replicate is recorded).
#' @param arch list or row with `dispersal` and `competition` ([archetypes_layerB()]).
#' @param sc list or row with `lambda`, `sigma_signal` and `h_signal` ([scenarios_layerA()]).
#' @param reps number of replicate metacommunities.
realized_richness <- function(n_pool, arch, sc, reps = 10) {
  mean(vapply(seq_len(reps), \(i) {
    mc <- sim_metacommunity_B(sim_tree_shape(n_pool), sc$lambda, sc$sigma_signal,
                              sc$h_signal, env_original = TRUE,
                              dispersal = arch$dispersal, competition = arch$competition)
    ncol(mc$comm)
  }, 0))
}

#' Pool size whose realized richness after the dynamics is close to `target`
#'
#' Starts at n = target and rescales n by target / realized until the mean
#' realized richness is within `tol` (relative) of the target, for at most
#' `max_iter` simulated pool sizes.
#' @return tibble with the last simulated pool size, its mean realized
#'   richness, the number of iterations and whether it converged.
calibrate_pool_size <- function(target, arch, sc, reps = 10, tol = 0.05,
                                max_iter = 5) {
  n <- target
  for (it in seq_len(max_iter)) {
    m <- realized_richness(n, arch, sc, reps)
    converged <- abs(m - target) / target <= tol
    if (converged || it == max_iter) break
    n <- round(n * target / max(m, 1))
  }
  tibble::tibble(target = target, n_pool = n, realized = m, iterations = it,
                 converged = converged)
}

#' Path of the calibrated pool sizes (written by analysis/10_calibrate_layerB.R)
layerB_pool_file <- function() file.path("data", "calibrated", "layerB_pool_sizes.csv")

#' Calibrated pool sizes of the feasible cells
#'
#' Reads the table written by analysis/10_calibrate_layerB.R and drops the cells
#' whose mean realized richness is more than `tol` (relative) off the target
#' (decision of 2026-10-18: mass effects with phylogenetic competition levels
#' off at ~150 species for pools of 2300-2900 tips, so N = 200 is infeasible
#' there). Cells within `tol` are kept although not converged, with their
#' realized N recorded per replicate.
#' @return tibble archetype, scenario, target, n_pool.
feasible_pools <- function(path, tol = 0.15) {
  pools <- readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::mutate(feasible = abs(realized - target) / target <= tol)
  dropped <- dplyr::filter(pools, !feasible)
  if (nrow(dropped) > 0)
    message("Dropped as infeasible (realized N more than ", 100 * tol,
            "% off target): ", nrow(dropped), " rows, ",
            paste(unique(paste(dropped$archetype, dropped$target)), collapse = ", "),
            "; realized ", paste(round(range(dropped$realized)), collapse = "-"))
  pools <- dplyr::filter(pools, feasible)
  if (!all(pools$converged))
    message("Kept although not converged (within ", 100 * tol, "% of target): ",
            sum(!pools$converged), " rows")
  dplyr::select(pools, archetype, scenario, target, n_pool)
}
