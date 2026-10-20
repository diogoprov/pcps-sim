# Phylogenetic community composition (PCC; Peres-Neto & Kembel 2015).
# Specification and reconstructed equations: docs/ESPEC_PCC.md.

#' Phylogenetic eigenvectors (PCoA of the patristic distance matrix)
#'
#' @param squared if FALSE (as printed in Peres-Neto & Kembel 2015, p. 712),
#'   B = -PD / 2; if TRUE, B = -PD^2 / 2 (classic PCoA).
#' @return species x (k - 1) matrix of eigenvectors with positive eigenvalues.
phylo_eigenvectors <- function(tree, squared = FALSE) {
  pd <- ape::cophenetic.phylo(tree)[tree$tip.label, tree$tip.label]
  b <- -0.5 * if (squared) pd^2 else pd
  b <- sweep(b, 1, rowMeans(b))
  b <- sweep(b, 2, colMeans(b))
  eig <- eigen(b, symmetric = TRUE)
  keep <- eig$values > sqrt(.Machine$double.eps) * max(eig$values)
  v <- eig$vectors[, keep, drop = FALSE]
  rownames(v) <- tree$tip.label
  v
}

#' Weighted standardization of the columns of x (weighted mean 0, variance 1)
#'
#' The variance denominator is the sum of weights, as in the paper.
weighted_std <- function(x, w) {
  x <- as.matrix(x)
  m <- colSums(w * x) / sum(w)
  xc <- sweep(x, 2, m)
  s <- sqrt(colSums(w * xc^2) / sum(w))
  sweep(xc, 2, s, "/")
}

#' PCC statistic: variance across eigenvectors of the adjusted R2
#'
#' @param comm sites x species (columns in the row order of `pev`).
#' @param pev species x (k - 1) eigenvector matrix.
#' @param env sites x p matrix (or vector) of predictors.
#' @return named vector: `var_r2` (test statistic) and `r2_pos` (mean of the
#'   positive adjusted R2, the value the authors recommend reporting).
pcc_stat <- function(comm, pev, env) {
  env <- as.matrix(env)
  n <- nrow(comm)
  p <- ncol(env)
  w_spp <- colSums(comm)
  w_site <- rowSums(comm)

  pev_std <- weighted_std(pev, w_spp)
  pcc <- (comm %*% pev_std) / w_site
  e_std <- weighted_std(env, w_site)

  # weighted regression of each PCC column on the standardized predictors
  xtwx <- crossprod(e_std, w_site * e_std)
  beta <- solve(xtwx, crossprod(e_std, w_site * pcc))
  fitted <- e_std %*% beta

  # tr(PL_i' PL_i) = sum of species weights (weighted variance 1), equal for all i
  r2 <- colSums(w_site * fitted^2) / sum(w_spp)
  r2_adj <- 1 - (n - 1) / (n - p - 1) * (1 - r2)

  c(var_r2 = stats::var(r2_adj), r2_pos = mean(r2_adj[r2_adj > 0]))
}

#' Row (site) and column (species) permutation tests of the PCC statistic
#'
#' With `site_env`, the permuted gradients replace the row permutations.
test_pcc <- function(comm, tree, env, nperm = 99, squared = FALSE, site_env = NULL) {
  comm <- comm[, tree$tip.label, drop = FALSE]
  pev <- phylo_eigenvectors(tree, squared = squared)
  obs <- pcc_stat(comm, pev, env)
  env <- as.matrix(env)

  rows <- if (is.null(site_env)) {
    replicate(nperm, {
      i <- sample.int(nrow(comm))
      pcc_stat(comm[i, , drop = FALSE], pev, env)[["var_r2"]]
    })
  } else {
    # restricted site permutations given as permuted gradients ([test_dp_all()])
    vapply(seq_len(nrow(site_env)), \(b) pcc_stat(comm, pev, as.matrix(site_env[b, ]))[["var_r2"]], 0)
  }
  cols <- replicate(nperm, {
    j <- sample.int(ncol(comm))
    pcc_stat(comm[, j, drop = FALSE], pev, env)[["var_r2"]]
  })

  p_row <- perm_p(obs[["var_r2"]], rows)
  p_col <- perm_p(obs[["var_r2"]], cols)
  tibble::tibble(
    test = "pcc", stat = obs[["r2_pos"]],
    p_site = p_row, p_taxa = p_col, p_max = max(p_row, p_col)
  )
}
