# Node scores of each method for the recovery of the planted clade
# (R/clade_recovery.R). PCPS, PCC and evoPCA give species scores, turned into
# node scores by [node_scores()]; phylofactorization gives an objective per
# split. All on the metacommunity as analysed (`mc` from
# [sim_metacommunity_clade()]). The ordinations are unsupervised, so each takes
# the axis most correlated with E among its first four as the axis "associated
# with E"; phylofactorization uses E directly in its objective.

#' Axis most correlated with the gradient among the first `n_axes`
best_axis <- function(axes, env) which.max(abs(stats::cor(axes, env)))

#' PCPS: correlation of each species (column of P) with the PCPS axis most
#' correlated with E (species scores as in `PCPS::pcps`)
pcps_species_scores <- function(comm, tree, env, n_axes = 4) {
  q <- phylo_q(tree)
  axes <- pcps_axes(comm, q, seq_len(n_axes))
  a <- axes[, best_axis(axes, env)]
  p <- matrix_p(comm, q)
  stats::setNames(drop(stats::cor(p, a)), tree$tip.label)
}

#' PCC: values of the phylogenetic eigenvector whose PCC column is most
#' correlated with E, among the first `n_axes` eigenvectors
#'
#' The analogue of the PCPS and evoPCA axes: one phylogenetic axis chosen by
#' its association with E. Combining all eigenvectors by their regression
#' coefficients would nearly span the whole species space, so the tree would
#' cancel out and the score would rank species by E alone.
pcc_species_scores <- function(comm, tree, env, n_axes = 4) {
  comm <- comm[, tree$tip.label, drop = FALSE]
  pev <- phylo_eigenvectors(tree)[, seq_len(n_axes), drop = FALSE]
  pev_std <- weighted_std(pev, colSums(comm))
  pcc <- (comm %*% pev_std) / rowSums(comm)
  stats::setNames(pev_std[, best_axis(pcc, env)], tree$tip.label)
}

#' evoPCA (Hellinger): correlation of each species with the evoPCA axis most
#' correlated with E (Pavoine 2016, `adiv::evopcahellinger`)
#'
#' Species scores, as for PCPS and PCC, from the species columns of the
#' transformed table. The branch columns are not used: with incidence data
#' (`abundance = FALSE`) a branch is present wherever any descendant is, so the
#' branch above a planted clade saturates and small subclades always correlate
#' better; and the loadings (`co`) scale with branch length, so long basal
#' branches would dominate.
evopca_species_scores <- function(comm, tree, env, n_axes = 4, incidence = TRUE) {
  e <- adiv::evopcahellinger(tree, as.data.frame(comm), scannf = FALSE, nf = n_axes,
                             abundance = !incidence)
  a <- best_axis(as.matrix(e$li), env)
  stats::setNames(drop(stats::cor(as.matrix(e$tab)[, tree$tip.label], e$li[, a])),
                  tree$tip.label)
}

#' Binomial phylofactorization: objective of every split
#'
#' The objective of `phylofactor::gpf(algorithm = "phylo")` with
#' cbind(Successes, Failures) ~ phylo * env and the default `pvDeviance`: the
#' sequential deviance of the terms with `phylo` (phylo, then phylo:env after
#' env). Computed on presences aggregated by site and side of the split, which
#' leaves the deviance differences unchanged (the model has no species effect),
#' so every node gets a score, not only the best one. Only internal nodes are
#' scored, as for the other methods. Note that the phylo term also rewards a
#' difference in prevalence unrelated to E, and that the logit is linear in E
#' although responses are unimodal; with the planted shift the response of the
#' clade is monotone along E.
phylofactor_node_scores <- function(comm, tree, env) {
  tips <- clade_tips(tree)
  n_site <- nrow(comm)
  k <- ncol(comm)
  # a split that separates presences perfectly fits probabilities of 0 or 1;
  # the deviance is still the right objective, so only that warning is muted;
  # a fit that does not converge gives no score
  dev <- function(y, n, x) {
    fit <- withCallingHandlers(
      stats::glm.fit(x, cbind(y, n - y), family = stats::binomial()),
      warning = \(w) if (grepl("fitted probabilities", conditionMessage(w)) ||
                         grepl("probabilidades ajustadas", conditionMessage(w)))
        invokeRestart("muffleWarning")
    )
    if (fit$converged) fit$deviance else NA_real_
  }
  score <- vapply(tips, \(t) {
    in_g <- colnames(comm) %in% t
    y <- c(rowSums(comm[, in_g, drop = FALSE]), rowSums(comm[, !in_g, drop = FALSE]))
    n <- rep(c(sum(in_g), k - sum(in_g)), each = n_site)
    g <- rep(c(1, 0), each = n_site)
    e <- rep(env, 2)
    one <- rep(1, 2 * n_site)
    d0 <- dev(y, n, cbind(one))
    d1 <- dev(y, n, cbind(one, g))
    d2 <- dev(y, n, cbind(one, g, e))
    d3 <- dev(y, n, cbind(one, g, e, g * e))
    (d0 - d1) + (d2 - d3)
  }, numeric(1))
  tibble::tibble(node = as.integer(names(tips)), score = unname(score))
}

#' Compositional phylofactorization (abundance): objective of every split
#'
#' The objective of `phylofactor::PhyloFactor(choice = "var")` with
#' Data ~ X: the explained variance (regression sum of squares over the number
#' of sites) of the isometric log-ratio balance of the split regressed on E.
#' Zeros are replaced by `pseudocount` before the log, as `PhyloFactor` needs
#' positive data.
phylofactor_ilr_node_scores <- function(comm, tree, env, pseudocount = 0.5) {
  comm <- comm[, tree$tip.label, drop = FALSE]
  lg <- log(comm + pseudocount)
  tips <- clade_tips(tree)
  score <- vapply(tips, \(t) {
    in_g <- colnames(comm) %in% t
    r <- sum(in_g)
    s <- sum(!in_g)
    y <- sqrt(r * s / (r + s)) *
      (rowMeans(lg[, in_g, drop = FALSE]) - rowMeans(lg[, !in_g, drop = FALSE]))
    sum((stats::lm.fit(cbind(1, env), y)$fitted.values - mean(y))^2) / length(y)
  }, numeric(1))
  tibble::tibble(node = as.integer(names(tips)), score = unname(score))
}

#' Recovery of the planted clade by every method
#'
#' @param mc output of [sim_metacommunity_clade()].
#' @param incidence whether `mc$comm` is 0/1; phylofactorization is binomial
#'   for incidence and compositional (log-ratio) for abundance.
#' @return tibble with one row per method: `method`, `hit`, `rank`, `jaccard`.
recover_clade_all <- function(mc, incidence = TRUE) {
  comm <- mc$comm
  tree <- mc$tree
  env <- mc$env
  scores <- list(
    pcps = node_scores(tree, pcps_species_scores(comm, tree, env)),
    pcc = node_scores(tree, pcc_species_scores(comm, tree, env)),
    evopca = node_scores(tree, evopca_species_scores(comm, tree, env, incidence = incidence)),
    phylofactor = if (incidence) phylofactor_node_scores(comm, tree, env)
                  else phylofactor_ilr_node_scores(comm, tree, env)
  )
  dplyr::bind_rows(lapply(scores, recovery_metrics, tree = tree,
                          planted_tips = mc$planted$tips), .id = "method")
}
