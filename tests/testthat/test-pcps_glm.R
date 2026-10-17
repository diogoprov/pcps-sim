set.seed(11)
tree <- ape::rphylo(30, birth = 1, death = 0)
env <- runif(40, 0, 100)
comm <- sim_comm_static(sim_optima_pnk(tree, signal = TRUE), env)
ok <- drop_empty_sites(comm, env)
comm <- ok$comm[, tree$tip.label]
# PCPS drops species absent from all sites (and prunes them from Q); do the same
tree <- ape::keep.tip(tree, colnames(comm)[colSums(comm) > 0])
comm <- comm[, tree$tip.label]
env <- ok$env
rownames(comm) <- paste0("s", seq_len(nrow(comm)))

test_that("observed GLM F on PCPS1 matches PCPS::pcps.sig with FUN.GLM", {
  ref <- PCPS::pcps.sig(comm, ape::cophenetic.phylo(tree), data.frame(env = env, row.names = rownames(comm)),
                        FUN = PCPS::FUN.GLM, formula = pcps.1 ~ env,
                        choices = 1, runs = 9)
  ours <- f_axes(pcps_axes(comm, phylo_q(tree), 1), env)
  expect_equal(unname(ours), ref$obs.statistic, tolerance = 1e-8)
})

test_that("symmetric Procrustes on one axis leaves the taxa-shuffle F unchanged", {
  q <- phylo_q(tree)
  obs_ax <- pcps_axes(comm, q, 1)
  j <- sample.int(ncol(comm))
  null_ax <- pcps_axes(comm, q[j, j], 1)
  fitted_ax <- stats::fitted(vegan::procrustes(obs_ax, null_ax, symmetric = TRUE))
  expect_equal(f_axes(fitted_ax, env), f_axes(null_ax, env), tolerance = 1e-8)
})

test_that("test_dp_all gives the same observed statistics as the separate tests", {
  joint <- test_dp_all(comm, tree, env, nperm = 19)
  sep <- dplyr::bind_rows(test_pfw(comm, tree, env, nperm = 19, sqrt_adonis = TRUE),
                          test_pcps_glm(comm, tree, env, nperm = 19))
  expect_equal(joint$test, sep$test)
  expect_equal(joint$stat, sep$stat, tolerance = 1e-8)
  expect_true(all(joint$p_max >= joint$p_site & joint$p_max >= joint$p_taxa))
})

test_that("extra axis tests: Procrustes and RDA", {
  set.seed(3)
  base <- test_dp_all(comm, tree, env, nperm = 19)
  set.seed(3)
  ext <- test_dp_all(comm, tree, env, nperm = 19, extra_axes = TRUE)
  # the extra tests use no random numbers, so the original tests are unchanged
  expect_equal(ext[seq_len(nrow(base)), ], base)
  expect_equal(ext$test[-seq_len(nrow(base))],
               c("glm_pcps1_proc", "glm_pcps2_proc", "rda_pcps"))
  # observed Procrustes F equals the plain F, and so do the site-shuffle p
  expect_equal(ext$stat[5:6], ext$stat[3:4], tolerance = 1e-8)
  expect_equal(ext$p_site[5:6], ext$p_site[3:4])
  # RDA F matches vegan::rda on the same axes
  y <- pcps_axes(comm, phylo_q(tree), 1:2)
  a <- anova(vegan::rda(y ~ env), permutations = 0)
  expect_equal(ext$stat[7], a$F[1], tolerance = 1e-8)
})
