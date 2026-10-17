# Phylogenetic fuzzy weighting (Pillar & Duarte 2010; Duarte et al. 2016)
# and permutation tests on D_P (Mantel, distance-based regression / ADONIS).

#' Matrix Q of phylogenetic fuzzy sets
#'
#' Equivalent to `SYNCSA::matrix.p()`: similarities s = 1 - d / max(d), each row
#' standardized to unit sum (S is symmetric, so rows = columns).
#' @param tree ultrametric `phylo`.
#' @return species x species matrix Q, in `tree$tip.label` order.
phylo_q <- function(tree) {
  d <- ape::cophenetic.phylo(tree)
  d <- d[tree$tip.label, tree$tip.label]
  s <- 1 - d / max(d)
  s / rowSums(s)
}

#' Matrix P (sites x species) of phylogeny-weighted composition
#' @param comm sites x species abundance or incidence matrix (columns in Q order).
#' @param q matrix Q from [phylo_q()].
matrix_p <- function(comm, q) {
  w <- comm / rowSums(comm)
  w %*% q
}

#' Bray-Curtis dissimilarities among rows of P (square-rooted by default)
d_p <- function(p, sqrt = TRUE) {
  d <- as.matrix(vegan::vegdist(p, method = "bray"))
  if (sqrt) base::sqrt(d) else d
}

#' Gower-centred matrix G = (I - 11'/n) (-D^2 / 2) (I - 11'/n)
gower_centre <- function(d) {
  a <- -0.5 * d^2
  a <- sweep(a, 1, rowMeans(a))
  sweep(a, 2, colMeans(a))
}

#' Distance-based regression statistic (McArdle & Anderson 2001) for one predictor
#'
#' Same pseudo-F and R2 as `vegan::adonis2(d ~ x)`, computed by quadratic forms so
#' that permutations are cheap.
#' @param g Gower-centred matrix.
#' @param x numeric predictor (length n).
db_stat <- function(g, x) {
  n <- length(x)
  e <- x - mean(x)
  ss_reg <- drop(crossprod(e, g %*% e)) / sum(e^2)
  ss_tot <- sum(diag(g))
  c(F = ss_reg / ((ss_tot - ss_reg) / (n - 2)), R2 = ss_reg / ss_tot)
}

#' Mantel statistic (Pearson) between two distance matrices
mantel_stat <- function(d1, d2) {
  lt <- lower.tri(d1)
  stats::cor(d1[lt], d2[lt])
}

#' Permutation p-value (observed value included in the null distribution)
perm_p <- function(obs, null) (sum(null >= obs) + 1) / (length(null) + 1)

#' Site- and taxa-shuffle tests of D_P against an environmental variable
#'
#' Site shuffle permutes sites relative to `env`. Taxa shuffle permutes species
#' across the tips of the tree (rows and columns of Q), keeping W fixed, and
#' recomputes P and D_P (Pillar & Duarte 2010; Duarte et al. 2016).
#' @param sqrt_adonis square-root the Bray-Curtis matrix used in the ADONIS test.
#'   Duarte et al. 2016 state sqrt for both tests, but their published ADONIS R2
#'   is only reproduced with plain Bray-Curtis (docs/MEMORIA_PROJETO.md).
#'   The Mantel test always uses sqrt(Bray-Curtis).
#' @return tibble with statistic, observed value and p-values (site, taxa, max).
test_pfw <- function(comm, tree, env, nperm = 999, sqrt_adonis = TRUE) {
  q <- phylo_q(tree)
  comm <- comm[, tree$tip.label, drop = FALSE]
  bc <- d_p(matrix_p(comm, q), sqrt = FALSE)
  dp <- sqrt(bc)
  de <- as.matrix(stats::dist(env))
  g <- gower_centre(if (sqrt_adonis) dp else bc)

  obs_db <- db_stat(g, env)
  obs_mantel <- mantel_stat(dp, de)

  site <- replicate(nperm, {
    i <- sample.int(length(env))
    c(db_stat(g, env[i])[["F"]], mantel_stat(dp, de[i, i]))
  })

  k <- ncol(comm)
  taxa <- replicate(nperm, {
    j <- sample.int(k)
    bc_null <- d_p(matrix_p(comm, q[j, j]), sqrt = FALSE)
    dp_null <- sqrt(bc_null)
    g_null <- gower_centre(if (sqrt_adonis) dp_null else bc_null)
    c(db_stat(g_null, env)[["F"]], mantel_stat(dp_null, de))
  })

  p_site <- c(perm_p(obs_db[["F"]], site[1, ]), perm_p(obs_mantel, site[2, ]))
  p_taxa <- c(perm_p(obs_db[["F"]], taxa[1, ]), perm_p(obs_mantel, taxa[2, ]))

  tibble::tibble(
    test = c("adonis", "mantel"),
    stat = c(obs_db[["R2"]], obs_mantel),
    p_site = p_site,
    p_taxa = p_taxa,
    p_max = pmax(p_site, p_taxa)
  )
}
