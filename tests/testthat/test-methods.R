set.seed(1)
tree <- ape::rphylo(30, birth = 1, death = 0)
env <- runif(40, 0, 100)
comm <- sim_comm_static(sim_optima_pnk(tree, signal = TRUE), env)
ok <- drop_empty_sites(comm, env)
comm <- ok$comm
env <- ok$env

test_that("matrix_p reproduces SYNCSA::matrix.p", {
  ref <- SYNCSA::matrix.p(comm, ape::cophenetic.phylo(tree)[colnames(comm), colnames(comm)],
                          notification = FALSE)$matrix.P
  ours <- matrix_p(comm[, tree$tip.label], phylo_q(tree))
  expect_equal(unname(ours), unname(ref[, tree$tip.label]), tolerance = 1e-10)
})

test_that("db_stat reproduces vegan::adonis2 F and R2", {
  dp <- d_p(matrix_p(comm[, tree$tip.label], phylo_q(tree)))
  ref <- vegan::adonis2(stats::as.dist(dp) ~ env, permutations = 0)
  ours <- db_stat(gower_centre(dp), env)
  expect_equal(ours[["F"]], ref$F[1], tolerance = 1e-8)
  expect_equal(ours[["R2"]], ref$R2[1], tolerance = 1e-8)
})

test_that("mantel_stat reproduces vegan::mantel", {
  dp <- d_p(matrix_p(comm[, tree$tip.label], phylo_q(tree)))
  de <- as.matrix(dist(env))
  ref <- vegan::mantel(stats::as.dist(dp), stats::as.dist(de), permutations = 0)
  expect_equal(mantel_stat(dp, de), ref$statistic, tolerance = 1e-10)
})

test_that("weighted_std gives weighted mean 0 and variance 1", {
  x <- matrix(rnorm(60), 30)
  w <- runif(30)
  z <- weighted_std(x, w)
  expect_equal(colSums(w * z) / sum(w), c(0, 0), tolerance = 1e-10)
  expect_equal(colSums(w * z^2) / sum(w), c(1, 1), tolerance = 1e-10)
})

test_that("PCC per-eigenvector R2 equals squared weighted correlation (one predictor)", {
  c2 <- comm[, tree$tip.label]
  pev <- phylo_eigenvectors(tree)
  r <- rowSums(c2)
  pcc <- (c2 %*% weighted_std(pev, colSums(c2))) / r
  e <- weighted_std(env, r)
  r_w <- colSums(r * drop(e) * pcc) / sum(r)  # pcc columns have weighted mean 0
  n <- nrow(c2)
  r2_adj <- 1 - (n - 1) / (n - 2) * (1 - r_w^2)
  expect_equal(pcc_stat(c2, pev, env)[["var_r2"]], var(r2_adj), tolerance = 1e-10)
})
