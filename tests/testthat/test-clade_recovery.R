set.seed(21)
tree <- ape::rcoal(60)

test_that("clade_tips matches extract.clade for every node", {
  ct <- clade_tips(tree)
  for (nd in names(ct)) {
    expect_setequal(ct[[nd]], ape::extract.clade(tree, as.integer(nd))$tip.label)
  }
})

test_that("the planted clade holds 10-30% of the tips and is not the root", {
  set.seed(1)
  for (i in 1:20) {
    p <- plant_clade(tree)
    expect_gte(length(p$tips), 6)
    expect_lte(length(p$tips), 18)
    expect_true(p$node > ape::Ntip(tree) + 1)
    expect_setequal(p$tips, ape::extract.clade(tree, p$node)$tip.label)
    expect_true(p$rel_age > 0 && p$rel_age < 1)
  }
})

test_that("old clades are older than young clades of the same size range", {
  set.seed(2)
  old <- replicate(20, plant_clade(tree, depth = "old")$rel_age)
  young <- replicate(20, plant_clade(tree, depth = "young")$rel_age)
  expect_gt(min(old), max(young))
})

test_that("no eligible clade gives an informative error", {
  expect_error(plant_clade(ape::rcoal(5)), "no clade")
})

test_that("the planted clade has shifted optima", {
  set.seed(3)
  mc <- sim_metacommunity_clade(tree, shift = 3)
  # traits keep the tip order of the full tree, minus the species dropped
  kept <- tree$tip.label[tree$tip.label %in% mc$tree$tip.label]
  in_clade <- kept %in% mc$planted$tips
  expect_gt(mean(mc$traits$u[in_clade]), mean(mc$traits$u[!in_clade]))
  expect_true(all(mc$traits$u >= -1 & mc$traits$u <= 101))
})

test_that("a score equal to clade membership ranks the planted node first", {
  set.seed(4)
  p <- plant_clade(tree)
  sp <- stats::setNames(as.numeric(tree$tip.label %in% p$tips) + rnorm(60, 0, 0.01),
                        tree$tip.label)
  m <- recovery_metrics(tree, node_scores(tree, sp), p$tips)
  expect_true(m$hit)
  expect_equal(m$rank, 1L)
  expect_equal(m$jaccard, 1)
})

test_that("the complement split of a root child counts as a hit", {
  root_kids <- tree$edge[tree$edge[, 1] == ape::Ntip(tree) + 1, 2]
  kid <- root_kids[root_kids > ape::Ntip(tree)][1]
  planted <- ape::extract.clade(tree, kid)$tip.label
  sp <- stats::setNames(as.numeric(tree$tip.label %in% planted), tree$tip.label)
  m <- recovery_metrics(tree, node_scores(tree, sp), planted)
  expect_true(m$hit)
  expect_equal(m$rank, 1L)
})

test_that("noise scores rarely hit", {
  set.seed(5)
  p <- plant_clade(tree)
  hits <- replicate(30, {
    sp <- stats::setNames(rnorm(60), tree$tip.label)
    recovery_metrics(tree, node_scores(tree, sp), p$tips)$hit
  })
  expect_lt(mean(hits), 0.2)
})

test_that("constant scores and lost clades give NA metrics", {
  set.seed(6)
  p <- plant_clade(tree)
  flat <- stats::setNames(rep(1, 60), tree$tip.label)
  expect_true(is.na(recovery_metrics(tree, node_scores(tree, flat), p$tips)$hit))
  pruned <- ape::keep.tip(tree, p$tips[1:3])
  sp <- stats::setNames(rnorm(3), pruned$tip.label)
  expect_true(is.na(recovery_metrics(pruned, node_scores(pruned, sp), p$tips)$hit))
})

test_that("recovery works on a tree that lost species", {
  set.seed(7)
  p <- plant_clade(tree)
  drop <- c(p$tips[1], setdiff(tree$tip.label, p$tips)[1:5])
  pruned <- ape::drop.tip(tree, drop)
  sp <- stats::setNames(as.numeric(pruned$tip.label %in% p$tips), pruned$tip.label)
  m <- recovery_metrics(pruned, node_scores(pruned, sp), p$tips)
  expect_true(m$hit)
})
