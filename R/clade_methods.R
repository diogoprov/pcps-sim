# Node scores of each method for the recovery of the planted clade
# (R/clade_recovery.R). PCPS, PCC and evoPCA give species scores, turned into
# node scores by [node_scores()]; phylofactorization gives an objective per
# split. All on the metacommunity as
# analysed (`mc` from [sim_metacommunity_clade()]).

#' Axis most correlated with the gradient among the first `n_axes`
best_axis <- function(axes, env) which.max(abs(stats::cor(axes, env)))

#' PCPS: correlation of each species (column of P) with the PCPS axis most
#' correlated with E (species scores as in `PCPS::pcps`)
pcps_species_scores <- function(comm, tree, env, n_axes = 4) {
  p <- matrix_p(comm, phylo_q(tree))
  axes <- pcoa_axes(gower_centre(d_p(p, sqrt = TRUE)), seq_len(n_axes))
  a <- axes[, best_axis(axes, env)]
  stats::setNames(drop(stats::cor(p, a)), tree$tip.label)
}

#' PCC: species values of the phylogenetic direction fitted to E
#'
#' The PCC columns (sites x eigenvectors) are regressed on E as in
#' [pcc_stat()]; the species score is the combination of the standardized
#' eigenvectors weighted by the regression coefficients.
pcc_species_scores <- function(comm, tree, env) {
  pev <- phylo_eigenvectors(tree)
  w_spp <- colSums(comm)
  w_site <- rowSums(comm)
  pev_std <- weighted_std(pev, w_spp)
  pcc <- (comm %*% pev_std) / w_site
  e_std <- weighted_std(env, w_site)
  beta <- crossprod(e_std, w_site * pcc) / sum(w_site * e_std^2)
  stats::setNames(drop(pev_std %*% t(beta)), tree$tip.label)
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
evopca_species_scores <- function(comm, tree, env, n_axes = 4) {
  e <- adiv::evopcahellinger(tree, as.data.frame(comm), scannf = FALSE, nf = n_axes,
                             abundance = FALSE)
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
#' leaves the deviance differences unchanged, so every node gets a score.
phylofactor_node_scores <- function(comm, tree, env) {
  tips <- clade_tips(tree)
  n_site <- nrow(comm)
  k <- ncol(comm)
  # a split that separates presences perfectly fits probabilities of 0 or 1;
  # the deviance is still the right objective, so the warning is muted
  dev <- function(y, n, x) {
    suppressWarnings(stats::glm.fit(x, cbind(y, n - y), family = stats::binomial()))$deviance
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

#' Recovery of the planted clade by every method
#'
#' @param mc output of [sim_metacommunity_clade()] (incidence data).
#' @return tibble with one row per method: `method`, `hit`, `rank`, `jaccard`.
recover_clade_all <- function(mc) {
  comm <- mc$comm
  tree <- mc$tree
  env <- mc$env
  scores <- list(
    pcps = node_scores(tree, pcps_species_scores(comm, tree, env)),
    pcc = node_scores(tree, pcc_species_scores(comm, tree, env)),
    evopca = node_scores(tree, evopca_species_scores(comm, tree, env)),
    phylofactor = phylofactor_node_scores(comm, tree, env)
  )
  dplyr::bind_rows(lapply(scores, recovery_metrics, tree = tree,
                          planted_tips = mc$planted$tips), .id = "method")
}
