test_that("a species alone converges to N* = (r - 1) / alpha", {
  set.seed(2)
  n <- mc_dynamics(matrix(5, 50, 1), matrix(0.05), dispersal = 0, n_steps = 300)
  expect_equal(mean(n), 80, tolerance = 0.05)
})

test_that("dispersal moves individuals without creating or losing them", {
  # one step without seeding: the same seed gives the same N before dispersal
  r <- matrix(1, 20, 3)
  a <- matrix(0, 3, 3)
  set.seed(9)
  n1 <- mc_dynamics(r, a, dispersal = 0, n_steps = 1, seed_steps = integer(0),
                    seed_mean = 3)
  set.seed(9)
  n2 <- mc_dynamics(r, a, dispersal = 0.5, n_steps = 1, seed_steps = integer(0),
                    seed_mean = 3)
  expect_equal(colSums(n1), colSums(n2))
  expect_false(isTRUE(all.equal(n1, n2)))
})

test_that("competition matrices: symmetric, unit diagonal, phylo ordered by relatedness", {
  set.seed(1)
  tr <- sim_tree_shape(30)
  for (type in c("none", "random", "phylo")) {
    a <- competition_matrix(tr, type)
    expect_true(isSymmetric(unname(a)))
    expect_equal(unname(diag(a)), rep(0.05, 30))
  }
  set.seed(3)
  ar <- competition_matrix(tr, "random")
  set.seed(3)
  ap <- competition_matrix(tr, "phylo")
  up <- upper.tri(ar)
  expect_equal(sort(ar[up]), sort(ap[up]))
  s <- (1 - ape::cophenetic.phylo(tr) / max(ape::cophenetic.phylo(tr)))[tr$tip.label, tr$tip.label]
  # pairs that share the same MRCA tie on s; across distinct s values the
  # coefficients are monotone
  rng <- tapply(ap[up], s[up], range)
  lo <- vapply(rng, `[`, 0, 1)
  hi <- vapply(rng, `[`, 0, 2)
  expect_true(all(lo[-1] >= hi[-length(hi)]))
})

test_that("sim_metacommunity_B returns the layer-A structure", {
  set.seed(6)
  tr <- sim_tree_shape(40)
  mc <- sim_metacommunity_B(tr, 1, TRUE, TRUE, TRUE, dispersal = 0.01,
                            competition = "phylo", n_steps = 200)
  expect_named(mc, c("comm", "tree", "env", "n_dropped_sites", "n_dropped_spp", "traits"))
  expect_identical(colnames(mc$comm), mc$tree$tip.label)
  expect_equal(nrow(mc$comm), length(mc$env))
  expect_true(all(mc$comm %in% c(0, 1)))
  expect_equal(ncol(mc$comm) + mc$n_dropped_spp, 40)
  expect_true(all(rowSums(mc$comm) > 0))
})

test_that("calibrate_pool_size returns the last pool size it simulated", {
  set.seed(8)
  arch <- list(dispersal = 0.01, competition = "phylo")
  sc <- list(lambda = 1, sigma_signal = FALSE, h_signal = FALSE)
  # tol = 0 cannot be met: the pool size must stay at the one simulated
  cal <- calibrate_pool_size(20, arch, sc, reps = 1, tol = 0, max_iter = 1)
  expect_equal(cal$n_pool, 20)
  expect_false(cal$converged)
  expect_equal(nrow(archetypes_layerB()), 4)
})
