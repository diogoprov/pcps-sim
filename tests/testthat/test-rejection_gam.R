set.seed(3)
n <- 4000
d <- tibble::tibble(gamma = runif(n, -6, 1), colless = runif(n, 0, 3),
                    archetype = sample(c("a", "b"), n, replace = TRUE),
                    n_spp = sample(c(50, 100), n, replace = TRUE))
# rejection rises with gamma and does not depend on colless
d$p_max <- ifelse(runif(n) < plogis(-2 + 0.6 * d$gamma + 1.5), 0.01, 0.5)

test_that("the GAM recovers a gamma effect and no colless effect", {
  fit <- fit_rejection_gam(d)
  s <- summarise_rejection_gam(fit)
  expect_lt(s$p_value[s$term == "s(gamma)"], 1e-6)
  expect_gt(s$p_value[s$term == "s(colless)"], 0.01)
})

test_that("the marginal rate follows the simulated curve", {
  fit <- fit_rejection_gam(d)
  m <- marginal_rejection_gamma(fit, d, n_grid = 5)
  expect_equal(m$rate, plogis(-0.5 + 0.6 * m$gamma), tolerance = 0.05)
  expect_true(all(diff(m$rate) > 0))
})

test_that("a single archetype or N drops that factor", {
  fit <- fit_rejection_gam(d[d$archetype == "a" & d$n_spp == 50, ])
  expect_false(any(c("archetype", "n_spp") %in% all.vars(formula(fit))))
})
