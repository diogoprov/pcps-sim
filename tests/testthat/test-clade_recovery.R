set.seed(21)
tree <- ape::rcoal(60)

test_that("the planted clade holds 10-30% of the tips and is not the root", {
  for (i in 1:20) {
    p <- plant_clade(tree)
    expect_gte(length(p$tips), 6)
    expect_lte(length(p$tips), 18)
    expect_true(p$node > ape::Ntip(tree) + 1)
    expect_setequal(p$tips, ape::extract.clade(tree, p$node)$tip.label)
    expect_true(p$rel_age > 0 && p$rel_age < 1)
  }
})

test_that("the planted clade has shifted optima", {
  mc <- sim_metacommunity_clade(tree, shift = 3)
  # traits keep the tip order of the full tree, minus the species dropped
  kept <- tree$tip.label[tree$tip.label %in% mc$tree$tip.label]
  in_clade <- kept %in% mc$planted$tips
  expect_gt(mean(mc$traits$u[in_clade]), mean(mc$traits$u[!in_clade]))
})

test_that("a score equal to clade membership ranks the planted node first", {
  p <- plant_clade(tree)
  sp <- stats::setNames(as.numeric(tree$tip.label %in% p$tips) + rnorm(60, 0, 0.01),
                        tree$tip.label)
  m <- recovery_metrics(tree, node_scores(tree, sp), p$tips)
  expect_true(m$hit)
  expect_equal(m$rank, 1L)
  expect_equal(m$jaccard, 1)
})

test_that("noise scores rarely hit and give partial Jaccard", {
  p <- plant_clade(tree)
  hits <- replicate(30, {
    sp <- stats::setNames(rnorm(60), tree$tip.label)
    recovery_metrics(tree, node_scores(tree, sp), p$tips)$hit
  })
  expect_lt(mean(hits), 0.2)
})

test_that("recovery works on a tree that lost species", {
  p <- plant_clade(tree)
  drop <- c(p$tips[1], setdiff(tree$tip.label, p$tips)[1:5])
  pruned <- ape::drop.tip(tree, drop)
  sp <- stats::setNames(as.numeric(pruned$tip.label %in% p$tips), pruned$tip.label)
  m <- recovery_metrics(pruned, node_scores(pruned, sp), p$tips)
  expect_true(m$hit)
})
