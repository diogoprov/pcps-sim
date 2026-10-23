context("sim_tree_shape")

test_that("usage", {
  set.seed(42)
  phy <- treestats::sim_tree_shape(n = 50)
  testthat::expect_true(inherits(phy, "phylo"))
  testthat::expect_equal(ape::Ntip(phy), 50)
  testthat::expect_true(ape::is.binary(phy))
  testthat::expect_true(ape::is.ultrametric(phy))
  testthat::expect_equal(max(ape::node.depth.edgelength(phy)), 1)
  testthat::expect_true(all(phy$edge.length >= 0))

  testthat::expect_error(treestats::sim_tree_shape(n = 1))
  testthat::expect_error(treestats::sim_tree_shape(n = 10, beta = -2))
  testthat::expect_error(treestats::sim_tree_shape(n = 10, c = Inf))
})

test_that("c moves gamma and beta moves Colless", {
  set.seed(1)
  shape <- function(beta, c) {
    trees <- replicate(100, treestats::sim_tree_shape(80, beta, c),
                       simplify = FALSE)
    c(gamma = mean(vapply(trees, treestats::gamma_statistic, 0)),
      colless = mean(vapply(trees, treestats::colless, 0,
                            normalization = "yule")))
  }
  yule <- shape(0, 0)
  low_c <- shape(0, -1)
  low_beta <- shape(-1.5, 0)
  # Yule process: gamma and the Yule-normalized Colless index near zero
  testthat::expect_lt(abs(yule[["gamma"]]), 0.3)
  testthat::expect_lt(abs(yule[["colless"]]), 0.3)
  # c changes gamma, not Colless
  testthat::expect_lt(low_c[["gamma"]], yule[["gamma"]] - 1)
  testthat::expect_lt(abs(low_c[["colless"]] - yule[["colless"]]), 0.3)
  # beta changes Colless, not gamma
  testthat::expect_gt(low_beta[["colless"]], yule[["colless"]] + 0.5)
  testthat::expect_lt(abs(low_beta[["gamma"]] - yule[["gamma"]]), 0.3)
})
