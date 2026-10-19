set.seed(31)
tree <- sim_tree_shape(40)
mc <- sim_metacommunity_clade(tree, shift = 3)

test_that("the best phylofactor split matches phylofactor::gpf", {
  skip_if_not_installed("phylofactor")
  comm <- mc$comm
  rownames(comm) <- paste0("s", seq_len(nrow(comm)))
  df <- data.table::data.table(Species = rep(colnames(comm), each = nrow(comm)),
                               Sample = rep(rownames(comm), ncol(comm)),
                               Successes = as.vector(comm),
                               env = rep(mc$env, ncol(comm)))
  df$Failures <- 1 - df$Successes
  pf <- suppressWarnings(phylofactor::gpf(df, mc$tree,
                                          frmla.phylo = cbind(Successes, Failures) ~ phylo * env,
                                          nfactors = 1, family = stats::binomial(),
                                          algorithm = "phylo"))
  # gpf warns about fitted probabilities of 0 or 1 for splits that separate
  # presences perfectly, as phylofactor_node_scores() does internally
  ref <- lapply(pf$groups[[1]], \(i) mc$tree$tip.label[i])
  s <- phylofactor_node_scores(mc$comm, mc$tree, mc$env)
  best <- ape::extract.clade(mc$tree, s$node[which.max(s$score)])$tip.label
  expect_true(setequal(best, ref[[1]]) || setequal(best, ref[[2]]))
})

test_that("species scores are named by tip and follow the shifted clade", {
  in_clade <- mc$tree$tip.label %in% mc$planted$tips
  for (sp in list(pcps_species_scores(mc$comm, mc$tree, mc$env),
                  pcc_species_scores(mc$comm, mc$tree, mc$env),
                  evopca_species_scores(mc$comm, mc$tree, mc$env))) {
    expect_named(sp, mc$tree$tip.label)
    expect_gt(abs(mean(sp[in_clade]) - mean(sp[!in_clade])), stats::sd(sp) / 2)
  }
})

test_that("a strong planted clade is recovered by evoPCA and phylofactor", {
  out <- recover_clade_all(mc)
  expect_setequal(out$method, c("pcps", "pcc", "evopca", "phylofactor"))
  expect_false(anyNA(out$jaccard))
  # PCPS spreads the signal over related species through the fuzzy weights of
  # P, and the first PCC eigenvectors follow the deep splits, so both pick a
  # larger clade here and are not required to hit this small recent clade
  expect_true(all(out$hit[out$method %in% c("evopca", "phylofactor")]))
})

test_that("PCPS species scores are the correlations of PCPS::pcps", {
  skip_if_not_installed("PCPS")
  comm <- mc$comm
  rownames(comm) <- paste0("s", seq_len(nrow(comm)))
  ref <- PCPS::pcps(comm, ape::cophenetic.phylo(mc$tree))
  sp <- pcps_species_scores(mc$comm, mc$tree, mc$env)
  # the axis is chosen by its correlation with E; signs of axes are arbitrary
  a <- which.max(abs(stats::cor(ref$vectors[, 1:4], mc$env)))
  expect_equal(abs(unname(sp)), abs(unname(ref$correlations[, a])), tolerance = 1e-8)
})

test_that("a species constant across sites does not break the scores", {
  tree5 <- ape::rcoal(12)
  sp <- stats::setNames(c(NA, rnorm(11)), tree5$tip.label)
  expect_false(anyNA(node_scores(tree5, sp)$score))
})
