#' Simulate a tree with separately controlled imbalance and branching times
#'
#' @description Simulates an ultrametric tree whose topology and branching
#' times come from two independent processes: the topology from Aldous'
#' beta-splitting model and the internode intervals from a pure-birth process
#' in which the splitting rate depends on the number of lineages. Because the
#' gamma statistic depends only on the internode intervals and the Colless
#' index only on the topology, `beta` mainly moves the imbalance of the tree and
#' `c` mainly moves its gamma statistic, so that the two can be varied
#' independently, for instance to build a grid of tree shapes for simulation
#' studies.
#'
#' @details
#' \strong{Topology.} A clade of \eqn{m} tips splits into clades of \eqn{i} and
#' \eqn{m - i} tips with probability proportional to
#' \deqn{\frac{\Gamma(\beta + i + 1)\, \Gamma(\beta + m - i + 1)}{\Gamma(i + 1)\, \Gamma(m - i + 1)}}{
#' Gamma(beta + i + 1) Gamma(beta + m - i + 1) / (Gamma(i + 1) Gamma(m - i + 1))}
#' (Aldous 1996). `beta = 0` gives the Yule (ERM) distribution of topologies,
#' `beta = -1.5` the PDA (uniform) distribution, values approaching -2
#' increasingly caterpillar-like trees and large values increasingly balanced
#' trees.
#'
#' \strong{Branching times.} With \eqn{k} lineages, each lineage splits at rate
#' \eqn{k^c}, so the interval during which there are \eqn{k} lineages is
#' exponential with rate \eqn{k^{1 + c}}{k^(1 + c)}. `c = 0` gives the internode
#' intervals of a Yule process; `c < 0` slows diversification as lineages
#' accumulate (as under diversity dependence), which gives negative values of
#' the gamma statistic (Pybus & Harvey 2000); `c > 0` gives positive values.
#' The intervals are rescaled so that the root depth is 1.
#'
#' \strong{Combining both.} The \eqn{n - 1} branching times are assigned to the
#' internal nodes in an order drawn uniformly among the orders compatible with
#' the topology (a uniform random linear extension of the tree order, which is
#' the order a Yule process implies given its topology): at each step, the next
#' node to split is drawn among the nodes whose parent has already split, with
#' probability proportional to the number of internal nodes in its subtree.
#'
#' \strong{Reaching target values.} `beta` and `c` set the expected shape, but
#' the realized gamma statistic and Colless index vary among trees and depend on
#' the number of tips. To obtain trees near given target values, calibrate the
#' mean of each statistic against its parameter for the tree size of interest
#' (for instance with [gamma_statistic()] and [colless()] on a few hundred trees
#' per parameter value), invert the curves, and keep only the simulated trees
#' whose statistics fall within a tolerance of the targets.
#'
#' @param n number of tips (at least 2).
#' @param beta parameter of the beta-splitting model (> -2); lower values give
#'   more imbalanced topologies.
#' @param c dependence of the splitting rate on the number of lineages; lower
#'   values give earlier branching times (more negative gamma).
#' @return an ultrametric object of class `phylo` with root depth 1 and tips
#'   labelled `t1`, ..., `tn`.
#' @references
#' Aldous, D. (1996). Probability distributions on cladograms. In D. Aldous &
#' R. Pemantle (Eds.), Random discrete structures (pp. 1-18). Springer.
#' \doi{10.1007/978-1-4612-0719-1_1}
#'
#' Pybus, O. G., & Harvey, P. H. (2000). Testing macro-evolutionary models
#' using incomplete molecular phylogenies. Proceedings of the Royal Society B,
#' 267, 2267-2272. \doi{10.1098/rspb.2000.1278}
#' @export
#' @examples
#' set.seed(1)
#' yule_like <- treestats::sim_tree_shape(n = 100)
#' early_imbalanced <- treestats::sim_tree_shape(n = 100, beta = -1.5, c = -1)
#' # gamma moves with c, Colless with beta
#' treestats::gamma_statistic(yule_like)
#' treestats::gamma_statistic(early_imbalanced)
#' treestats::colless(yule_like, normalization = "yule")
#' treestats::colless(early_imbalanced, normalization = "yule")
sim_tree_shape <- function(n, beta = 0, c = 0) {
  if (length(n) != 1 || !is.numeric(n) || n < 2 || n != round(n)) {
    stop("n has to be a single integer of at least 2")
  }
  if (length(beta) != 1 || !is.numeric(beta) || beta <= -2) {
    stop("beta has to be a single number larger than -2")
  }
  if (length(c) != 1 || !is.numeric(c) || !is.finite(c)) {
    stop("c has to be a single finite number")
  }
  topology <- rtopology_beta(n, beta)
  intervals <- intervals_lineage_rate(n, c)
  return(assign_ages(topology, intervals))
}

#' Random topology under the beta-splitting model
#' @param n number of tips
#' @param beta beta parameter
#' @return phylo object without branch lengths
#' @keywords internal
rtopology_beta <- function(n, beta) {
  split_clade <- function(tips) {
    m <- length(tips)
    if (m == 1) return(tips)
    i <- seq_len(m - 1)
    lw <- lgamma(beta + i + 1) + lgamma(beta + m - i + 1) -
          lgamma(i + 1) - lgamma(m - i + 1)
    k <- sample(i, 1, prob = exp(lw - max(lw)))
    paste0("(", split_clade(tips[seq_len(k)]), ",",
           split_clade(tips[-seq_len(k)]), ")")
  }
  return(ape::read.tree(text = paste0(split_clade(paste0("t", seq_len(n))), ";")))
}

#' Internode intervals of a pure-birth process with lineage-dependent rate
#' @param n number of tips
#' @param c dependence of the splitting rate on the number of lineages
#' @return vector of n - 1 intervals summing to 1, from 2 to n lineages
#' @keywords internal
intervals_lineage_rate <- function(n, c) {
  k <- 2:n
  g <- stats::rexp(n - 1, rate = k^(1 + c))
  return(g / sum(g))
}

#' Assign branching times to the internal nodes of a topology
#' @param topology binary rooted phylo object
#' @param intervals output of intervals_lineage_rate()
#' @return phylo object with branch lengths
#' @keywords internal
assign_ages <- function(topology, intervals) {
  n <- ape::Ntip(topology)
  ages_sorted <- rev(cumsum(rev(intervals)))  # oldest split first
  parent <- topology$edge[, 1]
  child <- topology$edge[, 2]
  n_internal <- vapply(ape::prop.part(topology), length, 0) - 1
  age <- numeric(n + topology$Nnode)
  available <- n + 1
  for (r in seq_len(n - 1)) {
    v <- if (length(available) == 1) available else
      sample(available, 1, prob = n_internal[available - n])
    age[v] <- ages_sorted[r]
    kids <- child[parent == v]
    available <- c(setdiff(available, v), kids[kids > n])
  }
  topology$edge.length <- age[parent] - age[child]
  return(topology)
}
