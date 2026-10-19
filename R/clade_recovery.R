# Objective 2: recovery of a clade associated with the gradient. One clade is
# planted: the optima of its species are shifted along E, and the other species
# have optima without phylogenetic signal. Each method gives a score per node,
# and the recovery metrics compare the best-scoring node with the planted one.

#' Draw the clade to plant
#'
#' Uniform among the internal nodes (root excluded) whose clade holds between
#' `min_frac` and `max_frac` of the tips.
#' @return list with `node`, `tips` (labels of the clade) and `rel_age` (node
#'   height over root height: near 1 for basal clades, near 0 for recent ones).
plant_clade <- function(tree, min_frac = 0.1, max_frac = 0.3) {
  k <- ape::Ntip(tree)
  nodes <- seq(k + 2, k + tree$Nnode)
  sizes <- vapply(nodes, \(nd) length(ape::extract.clade(tree, nd)$tip.label), integer(1))
  eligible <- nodes[sizes >= min_frac * k & sizes <= max_frac * k]
  stopifnot(length(eligible) > 0)
  node <- eligible[sample.int(length(eligible), 1)]
  depth <- ape::node.depth.edgelength(tree)
  list(node = node, tips = ape::extract.clade(tree, node)$tip.label,
       rel_age = 1 - depth[node] / max(depth))
}

#' One metacommunity with a planted clade
#'
#' Optima: iid N(0, 1) plus `shift` for the species of the planted clade,
#' rescaled to [-1, 101]; niche breadth and maximum abundance without signal,
#' as in [sim_species_traits()]. The gradient analysed is the one that
#' generated the data.
#' @return [finish_metacommunity()] output plus `planted` ([plant_clade()]).
sim_metacommunity_clade <- function(tree, shift = 2, incidence = TRUE, n_sites = 100,
                                    mu_tol = 10, min_frac = 0.1, max_frac = 0.3) {
  planted <- plant_clade(tree, min_frac, max_frac)
  tr <- sim_species_traits(tree, 0, FALSE, FALSE, mu_tol)
  in_clade <- tree$tip.label %in% planted$tips
  tr$u <- rescale_to(stats::rnorm(ape::Ntip(tree)) + shift * in_clade)
  env <- stats::runif(n_sites, 0, 100)
  mc <- finish_metacommunity(sim_comm_traits(tr, env, tree$tip.label, incidence), tree,
                             env, TRUE, tr)
  c(mc, list(planted = planted))
}

#' Score of every internal node for a vector of species scores
#'
#' |correlation| between membership of the clade (0/1) and the species scores;
#' the root (all species) has no score.
#' @param sp_score named vector, one value per tip of `tree`.
#' @return tibble with `node` and `score`.
node_scores <- function(tree, sp_score) {
  k <- ape::Ntip(tree)
  sp_score <- sp_score[tree$tip.label]
  nodes <- seq(k + 2, k + tree$Nnode)
  score <- vapply(nodes, \(nd) {
    member <- tree$tip.label %in% ape::extract.clade(tree, nd)$tip.label
    abs(stats::cor(member, sp_score))
  }, numeric(1))
  tibble::tibble(node = nodes, score = score)
}

#' Recovery of the planted clade
#'
#' The planted clade is identified in `tree` (which may have lost species) as
#' the most recent common ancestor of its species that are present. A node and
#' its complement define the same split, so the match uses the Jaccard index
#' between the best-scoring clade and either the planted clade or its
#' complement.
#' @param scores output of [node_scores()].
#' @param planted_tips labels of the planted clade.
#' @return one-row tibble: `hit` (best node is the planted one), `rank` of the
#'   planted node among all nodes (1 = best) and `jaccard`.
recovery_metrics <- function(tree, scores, planted_tips) {
  present <- intersect(planted_tips, tree$tip.label)
  stopifnot(length(present) >= 2)
  true_node <- ape::getMRCA(tree, present)
  best <- scores$node[which.max(scores$score)]
  best_tips <- ape::extract.clade(tree, best)$tip.label
  jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  j <- max(jaccard(best_tips, present),
           jaccard(best_tips, setdiff(tree$tip.label, present)))
  rank <- if (true_node %in% scores$node) {
    sum(scores$score > scores$score[scores$node == true_node]) + 1L
  } else NA_integer_
  tibble::tibble(hit = j == 1, rank = rank, jaccard = j)
}
