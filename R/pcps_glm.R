# Linear models on PCPS axes (Duarte et al. 2016; `PCPS::pcps.sig` with
# `PCPS::FUN.GLM`), one axis at a time, with site and taxa shuffles.

#' PCPS axes: PCoA of D_P (square-rooted Bray-Curtis of matrix P)
#'
#' Same space as `PCPS::pcps()` (which drops negative eigenvalues); the scaling
#' of the axes is irrelevant for the F statistic.
#' @param axes which axes to return.
pcps_axes <- function(comm, q, axes = 1:2) {
  g <- gower_centre(d_p(matrix_p(comm, q), sqrt = TRUE))
  # leading eigenpairs only (Lanczos); ~10x faster than eigen() for n = 100
  eig <- RSpectra::eigs_sym(g, max(axes))
  sweep(eig$vectors[, axes, drop = FALSE], 2, sqrt(pmax(eig$values[axes], 0)), "*")
}

#' F statistic of the regression of each column of y on x (one predictor)
f_axes <- function(y, x) {
  n <- length(x)
  e <- x - mean(x)
  yc <- sweep(y, 2, colMeans(y))
  ss_reg <- drop(crossprod(e, yc))^2 / sum(e^2)
  ss_tot <- colSums(yc^2)
  ss_reg / ((ss_tot - ss_reg) / (n - 2))
}

#' GLM (Gaussian, identity) of PCPS axes on an environmental variable
#'
#' Site shuffle permutes `env`. Taxa shuffle permutes species across the tips
#' (rows and columns of Q), recomputes P and the PCPS, and refits. `pcps.sig`
#' first rotates the null axes onto the observed ones by symmetric Procrustes;
#' for a single axis that is a sign and scale change, which leaves F unchanged,
#' so it is skipped here (tests/testthat/test-pcps_glm.R).
#' @return tibble with one row per axis.
test_pcps_glm <- function(comm, tree, env, axes = 1:2, nperm = 999) {
  q <- phylo_q(tree)
  comm <- comm[, tree$tip.label, drop = FALSE]
  y <- pcps_axes(comm, q, axes)
  obs <- f_axes(y, env)

  site <- replicate(nperm, f_axes(y, env[sample.int(length(env))]))
  k <- ncol(comm)
  taxa <- replicate(nperm, {
    j <- sample.int(k)
    f_axes(pcps_axes(comm, q[j, j], axes), env)
  })
  site <- matrix(site, nrow = length(axes))
  taxa <- matrix(taxa, nrow = length(axes))

  p_site <- vapply(seq_along(axes), \(a) perm_p(obs[a], site[a, ]), 0)
  p_taxa <- vapply(seq_along(axes), \(a) perm_p(obs[a], taxa[a, ]), 0)
  tibble::tibble(
    test = paste0("glm_pcps", axes),
    stat = obs,
    p_site = p_site,
    p_taxa = p_taxa,
    p_max = pmax(p_site, p_taxa)
  )
}

#' F of the multivariate regression (RDA) of the columns of y on x
f_rda <- function(y, x) {
  n <- length(x)
  e <- x - mean(x)
  yc <- sweep(y, 2, colMeans(y))
  ss_reg <- sum(drop(crossprod(e, yc))^2) / sum(e^2)
  ss_tot <- sum(yc^2)
  ss_reg / ((ss_tot - ss_reg) / (n - 2))
}

#' ADONIS, Mantel and GLM on PCPS axes with shared permutations
#'
#' Same tests as [test_pfw()] (sqrt Bray-Curtis for both) and
#' [test_pcps_glm()], but each site or taxa permutation is used by all of them,
#' so P and D_P are computed once per taxa permutation.
#'
#' With `extra_axes = TRUE`, two more tests on the PCPS axes:
#' - `glm_pcps<a>_proc`: GLM on each axis after symmetric Procrustes of the
#'   null axes onto the observed ones, all `axes` at once (`PCPS::pcps.sig`
#'   with `choices = axes`). Site shuffle leaves the axes unchanged, so its p
#'   equals that of the plain GLM.
#' - `rda_pcps`: RDA of all `axes` on env (invariant to rotations of the axes,
#'   so it needs no Procrustes).
#' It also adds the share of each axis in the total variation (`eig_share<a>`).
test_dp_all <- function(comm, tree, env, axes = 1:2, nperm = 999,
                        extra_axes = FALSE) {
  q <- phylo_q(tree)
  comm <- comm[, tree$tip.label, drop = FALSE]
  de <- as.matrix(stats::dist(env))
  pcoa_axes <- function(g) {
    eig <- RSpectra::eigs_sym(g, max(axes))
    sweep(eig$vectors[, axes, drop = FALSE], 2, sqrt(pmax(eig$values[axes], 0)), "*")
  }
  dp <- d_p(matrix_p(comm, q), sqrt = TRUE)
  g <- gower_centre(dp)
  y <- pcoa_axes(g)
  stats_dp <- function(dp, x, de_x, g = gower_centre(dp), y_null = pcoa_axes(g)) {
    out <- c(db_stat(g, x)[["F"]], mantel_stat(dp, de_x), f_axes(y_null, x))
    if (extra_axes) {
      y_proc <- stats::fitted(vegan::procrustes(y, y_null, symmetric = TRUE))
      out <- c(out, f_axes(y_proc, x), f_rda(y_null, x))
    }
    out
  }
  obs <- stats_dp(dp, env, de, g, y)

  site <- replicate(nperm, {
    i <- sample.int(length(env))
    out <- c(db_stat(g, env[i])[["F"]], mantel_stat(dp, de[i, i]), f_axes(y, env[i]))
    if (extra_axes) out <- c(out, f_axes(y, env[i]), f_rda(y, env[i]))
    out
  })
  k <- ncol(comm)
  taxa <- replicate(nperm, {
    j <- sample.int(k)
    stats_dp(d_p(matrix_p(comm, q[j, j]), sqrt = TRUE), env, de)
  })
  p_site <- vapply(seq_along(obs), \(a) perm_p(obs[a], site[a, ]), 0)
  p_taxa <- vapply(seq_along(obs), \(a) perm_p(obs[a], taxa[a, ]), 0)
  tests <- c("adonis", "mantel", paste0("glm_pcps", axes))
  if (extra_axes) tests <- c(tests, paste0("glm_pcps", axes, "_proc"), "rda_pcps")
  out <- tibble::tibble(
    test = tests,
    stat = c(db_stat(g, env)[["R2"]], obs[-1]),
    p_site = p_site,
    p_taxa = p_taxa,
    p_max = pmax(p_site, p_taxa)
  )
  if (extra_axes) {
    # share of each PCPS in the total variation of D_P (the eigenvalue share of
    # the old tree-shape design, legacy notes 2012); trace(G) is the sum of
    # the eigenvalues, all non-negative for sqrt Bray-Curtis
    share <- colSums(y^2) / sum(diag(g))
    for (a in seq_along(axes)) out[[paste0("eig_share", axes[a])]] <- share[a]
  }
  out
}
