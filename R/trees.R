# Trees with separately controlled branching-time shape (Pybus & Harvey's gamma)
# and topological imbalance (Colless). Gamma depends only on the internode
# intervals and Colless only on the topology, so drawing them from independent
# processes gives (near) independent control. Design: docs/PLANO_SIMULACAO.md,
# section 2.

#' Random topology under Aldous' beta-splitting model
#'
#' A clade of m tips splits into i and m - i tips with probability proportional
#' to Gamma(beta + i + 1) Gamma(beta + m - i + 1) / (Gamma(i + 1) Gamma(m - i + 1)).
#' beta = 0 is Yule, beta = -1.5 is PDA, beta -> -2 is caterpillar-like and
#' large beta is balanced.
#' @return `phylo` without branch lengths; tips t1..tn.
rtopology_beta <- function(n, beta = 0) {
  stopifnot(n >= 2, beta > -2)
  split_clade <- function(tips) {
    m <- length(tips)
    if (m == 1) return(tips)
    i <- seq_len(m - 1)
    lw <- lgamma(beta + i + 1) + lgamma(beta + m - i + 1) - lgamma(i + 1) - lgamma(m - i + 1)
    k <- sample(i, 1, prob = exp(lw - max(lw)))
    paste0("(", split_clade(tips[seq_len(k)]), ",", split_clade(tips[-seq_len(k)]), ")")
  }
  ape::read.tree(text = paste0(split_clade(paste0("t", seq_len(n))), ";"))
}

#' Internode intervals of a pure-birth process with lineage-dependent rates
#'
#' With k lineages each splits at rate k^c, so the interval with k lineages is
#' Exp(k^(1 + c)). c = 0 is Yule; c < 0 slows diversification as lineages
#' accumulate (diversity-dependence-like, gamma < 0); c > 0 accelerates it
#' (gamma > 0). The final interval (n lineages, up to the present) is drawn
#' from the same law.
#' @return intervals g_2, ..., g_n (length n - 1), scaled to total depth 1.
intervals_lineage_rate <- function(n, c = 0) {
  k <- 2:n
  g <- stats::rexp(n - 1, rate = k^(1 + c))
  g / sum(g)
}

#' Assign node ages to a topology from internode intervals
#'
#' Internal nodes are ranked by a uniform random linear extension of the tree
#' order (the ranking a Yule process implies given the topology): the next node
#' to split is drawn among the available ones with probability proportional to
#' the number of internal nodes in its subtree.
#' @param topo `phylo` topology (binary, rooted).
#' @param g intervals from [intervals_lineage_rate()] (length Ntip - 1).
assign_ages <- function(topo, g) {
  n <- ape::Ntip(topo)
  stopifnot(length(g) == n - 1, ape::is.binary(topo))
  ages_sorted <- rev(cumsum(rev(g)))      # ages of the n - 1 splits, oldest first
  parent <- topo$edge[, 1]
  child <- topo$edge[, 2]
  n_int <- vapply(ape::prop.part(topo), length, 0) - 1   # internal nodes per subtree
  age <- numeric(n + topo$Nnode)
  available <- n + 1
  for (r in seq_len(n - 1)) {
    v <- if (length(available) == 1) available else
      sample(available, 1, prob = n_int[available - n])
    age[v] <- ages_sorted[r]
    kids <- child[parent == v]
    available <- c(setdiff(available, v), kids[kids > n])
  }
  topo$edge.length <- age[parent] - age[child]
  topo
}

#' Ultrametric tree with root depth 1 and given shape parameters
#' @param beta Aldous' beta (imbalance; lower = more imbalanced).
#' @param c lineage-dependence of the split rate (gamma; lower = earlier splits).
sim_tree_shape <- function(n, beta = 0, c = 0) {
  assign_ages(rtopology_beta(n, beta), intervals_lineage_rate(n, c))
}

#' Shape summaries: Pybus & Harvey's gamma and Colless index normalized against
#' its Yule expectation (`treestats`), comparable across tree sizes.
tree_shape <- function(tree) {
  c(gamma = treestats::gamma_statistic(tree),
    colless = treestats::colless(tree, normalization = "yule"))
}

#' All `treestats` summary statistics for one tree, raw and normalized
#'
#' `treestats::calc_all_stats()`; normalized versions (suffix `_norm`) are given
#' only for the statistics that `treestats` normalizes (size-standardized).
#' @return one-row tibble.
tree_stats_all <- function(tree) {
  raw <- unlist(treestats::calc_all_stats(tree, normalize = FALSE))
  nrm <- unlist(treestats::calc_all_stats(tree, normalize = TRUE))
  changed <- names(nrm)[abs(nrm - raw) > 1e-12 * pmax(1, abs(raw))]
  tibble::as_tibble(as.list(c(raw, stats::setNames(nrm[changed], paste0(changed, "_norm")))))
}
