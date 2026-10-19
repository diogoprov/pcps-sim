set.seed(3)
n <- 4000
d <- tibble::tibble(gamma = runif(n, -6, 1), colless = runif(n, 0, 3))
# rejection rises with gamma and does not depend on colless
d$p_max <- ifelse(runif(n) < plogis(-0.5 + 0.6 * d$gamma), 0.01, 0.5)
fit <- fit_rejection_gam(d)

test_that("the GAM recovers a gamma effect and no colless effect", {
  s <- summarise_rejection_gam(fit)
  expect_lt(s$p_value[s$term == "s(gamma)"], 1e-6)
  expect_gt(s$p_value[s$term == "s(colless)"], 0.01)
})

test_that("the curve follows the simulated rate and brackets it", {
  m <- rejection_curve_gamma(fit, n_grid = 5)
  truth <- plogis(-0.5 + 0.6 * m$gamma)
  expect_equal(m$rate, truth, tolerance = 0.05)
  expect_true(all(m$lower < m$rate & m$rate < m$upper))
  expect_true(all(diff(m$rate) > 0))
})
