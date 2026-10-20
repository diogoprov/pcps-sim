set.seed(21)
tree <- sim_tree_shape(40, beta = -0.5, c = -0.8)

test_that("rank_to_marginal keeps the target values and the order of x", {
  x <- rnorm(40)
  target <- runif(40, 0, 30)
  y <- rank_to_marginal(x, target)
  expect_equal(sort(unname(y)), sort(target))
  expect_equal(order(y), order(x))
})

test_that("Pagel lambda = 0 gives a diagonal covariance (no signal)", {
  x <- replicate(2000, sim_trait_lambda(tree, 0)[1:2])
  expect_lt(abs(cor(x[1, ], x[2, ])), 0.08)
})

test_that("sim_metacommunity_A returns aligned comm, tree and env", {
  for (sc in seq_len(nrow(scenarios_layerA()))) {
    s <- scenarios_layerA()[sc, ]
    mc <- sim_metacommunity_A(tree, s$lambda, s$sigma_signal, s$h_signal,
                              s$env_original, incidence = sc %% 2 == 0)
    expect_identical(colnames(mc$comm), mc$tree$tip.label)
    expect_equal(nrow(mc$comm), length(mc$env))
    expect_true(all(rowSums(mc$comm) > 0), all(colSums(mc$comm) > 0))
    expect_equal(ape::Ntip(mc$tree) + mc$n_dropped_spp, 40)
  }
})

test_that("uniform optima share one marginal distribution and keep the BM order", {
  set.seed(3)
  tr <- ape::rcoal(40)
  set.seed(11)
  a <- sim_species_traits(tr, 1, FALSE, FALSE, optima = "uniform")
  expect_equal(sort(a$u), seq(-1, 101, length.out = 40))
  set.seed(11)
  b <- sim_species_traits(tr, 1, FALSE, FALSE)
  expect_equal(rank(a$u), rank(b$u))
  expect_equal(range(b$u), c(-1, 101))
})
