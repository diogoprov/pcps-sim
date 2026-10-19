# Objective 2: recovery of a clade associated with the gradient. One clade is
# planted: the optima of its species are shifted along E, and the other species
# have optima without phylogenetic signal. Each method gives a score per node,
# and the recovery metrics compare the best-scoring node with the planted one.

#' Tip labels of the clade of every internal node except the root
#'
#' @return named list, one element per node (names = node numbers).
clade_tips <- function(tree) {
  parts <- ape::prop.part(tree)
  k <- ape::Ntip(tree)
  out <- lapply(parts[-1], \(p) tree$tip.label[p])
  stats::setNames(out, seq(k + 2, k + tree$Nnode))
}

#' Draw the clade to plant
#'
#' Uniform among the internal nodes (root excluded) whose clade holds between
#' `min_frac` and `max_frac` of the tips. With `depth = "old"` or `"young"`,
#' only among the eligible nodes at or above, or below, their median age, so
#' that basal and recent clades of the same size can be compared (H3); with
#' fixed clade sizes most eligible nodes are recent, so an unrestricted draw
#' would rarely plant a basal clade. Ages assume an ultrametric tree.
#' @return list with `node`, `tips` (labels of the clade) and `rel_age` (node
#'   height over root height: near 1 for basal clades, near 0 for recent ones).
plant_clade <- function(tree, min_frac = 0.1, max_frac = 0.3,
                        depth = c("any", "old", "young")) {
  depth <- match.arg(depth)
  k <- ape::Ntip(tree)
  tips <- clade_tips(tree)
  sizes <- lengths(tips)
  height <- ape::node.depth.edgelength(tree)
  rel_age <- 1 - height / max(height)
  eligible <- as.integer(names(tips))[sizes >= min_frac * k & sizes <= max_frac * k]
  if (length(eligible) == 0)
    stop(sprintf("no clade with %g-%g of the %d tips", min_frac, max_frac, k))
  if (depth != "any") {
    old <- rel_age[eligible] >= stats::median(rel_age[eligible])
    eligible <- eligible[if (depth == "old") old else !old]
  }
  node <- eligible[sample.int(length(eligible), 1)]
  list(node = node, tips = tips[[as.character(node)]], rel_age = rel_age[node])
}

#' One metacommunity with a planted clade
#'
#' Optima: u = -1 + 102 * Phi(z + shift * [species in the planted clade]), with
#' z iid N(0, 1), so the other species have no phylogenetic signal and spread
#' over the gradient, and the shift has the same meaning for every tree size
#' (a rescaling by the sample range would make it depend on N). Niche breadth
#' and maximum abundance without signal, as in [sim_species_traits()]. The
#' gradient analysed is the one that generated the data.
#' @param depth passed to [plant_clade()].
#' @return [finish_metacommunity()] output plus `planted` ([plant_clade()]).
sim_metacommunity_clade <- function(tree, shift = 2, depth = "any", incidence = TRUE,
                                    n_sites = 100, mu_tol = 10, min_frac = 0.1,
                                    max_frac = 0.3) {
  planted <- plant_clade(tree, min_frac, max_frac, depth)
  tr <- sim_species_traits(tree, lambda = 0, sigma_signal = FALSE, h_signal = FALSE,
                           mu_tol = mu_tol)
  in_clade <- tree$tip.label %in% planted$tips
  tr$u <- -1 + 102 * stats::pnorm(stats::rnorm(ape::Ntip(tree)) + shift * in_clade)
  env <- stats::runif(n_sites, 0, 100)
  mc <- finish_metacommunity(sim_comm_traits(tr, env, tree$tip.label, incidence), tree,
                             env, env_original = TRUE, traits = tr)
  c(mc, list(planted = planted))
}

#' Score of every internal node for a vector of species scores
#'
#' |correlation| between membership of the clade (0/1) and the species scores;
#' the root (all species) has no score. The two children of the root define the
#' same split and get the same score.
#' @param sp_score named vector with one value for every tip of `tree`; NA
#'   (a species constant across sites) counts as 0, no association.
#' @return tibble with `node` and `score` (NA if the species scores are
#'   constant).
node_scores <- function(tree, sp_score) {
  stopifnot(all(tree$tip.label %in% names(sp_score)))
  sp_score <- sp_score[tree$tip.label]
  sp_score[is.na(sp_score)] <- 0
  tips <- clade_tips(tree)
  score <- if (stats::sd(sp_score) > 0) {
    vapply(tips, \(t) abs(stats::cor(tree$tip.label %in% t, sp_score)), numeric(1))
  } else rep(NA_real_, length(tips))
  tibble::tibble(node = as.integer(names(tips)), score = unname(score))
}

#' Recovery of the planted clade
#'
#' The planted clade is identified in `tree` (which may have lost species) as
#' the most recent common ancestor of its species that are present. A node and
#' its complement define the same split, so the match uses the Jaccard index
#' between the best-scoring clade and either the planted clade or its
#' complement; ties in the rank count in favour of the planted node.
#' @param scores output of [node_scores()].
#' @param planted_tips labels of the planted clade.
#' @return one-row tibble: `hit` (best clade is the planted split), `rank` of
#'   the planted node among all nodes (1 = best) and `jaccard`; all NA when
#'   fewer than two planted species remain, the planted clade became the root
#'   or the scores are NA.
recovery_metrics <- function(tree, scores, planted_tips) {
  present <- intersect(planted_tips, tree$tip.label)
  if (length(present) < 2 || all(is.na(scores$score)))
    return(tibble::tibble(hit = NA, rank = NA_integer_, jaccard = NA_real_))
  true_node <- ape::getMRCA(tree, present)
  # the planted split is lost if pruning left only planted species
  if (!true_node %in% scores$node)
    return(tibble::tibble(hit = NA, rank = NA_integer_, jaccard = NA_real_))
  best <- scores$node[which.max(scores$score)]
  best_tips <- ape::extract.clade(tree, best)$tip.label
  jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  j <- max(jaccard(best_tips, present),
           jaccard(best_tips, setdiff(tree$tip.label, present)))
  rank <- sum(scores$score > scores$score[scores$node == true_node], na.rm = TRUE) + 1L
  tibble::tibble(hit = j == 1, rank = rank, jaccard = j)
}
