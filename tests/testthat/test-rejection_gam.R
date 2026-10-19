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

test_that("grid points away from the observed gamma get no rate", {
  gap <- d[d$gamma < -4 | d$gamma > -1, ]
  m <- rejection_curve_gamma(fit_rejection_gam(gap), n_grid = 25)
  inside <- m$gamma > -3.7 & m$gamma < -1.3
  expect_true(all(is.na(m$rate[inside])))
  expect_false(anyNA(m$rate[!inside & (m$gamma < -4.3 | m$gamma > -0.7)]))
})

test_that("the within-cell slope recovers gamma effects inside cells only", {
  # two cells far apart in gamma, with gamma varying within each
  set.seed(4)
  w <- tibble::tibble(cell = rep(c("low", "high"), each = 3000),
                      gamma = c(rnorm(3000, -6, 0.6), rnorm(3000, 0, 0.6)),
                      colless = runif(6000, 0, 3))
  # rejection jumps between cells but is flat in gamma within them
  w$p_max <- ifelse(runif(6000) < ifelse(w$cell == "low", 0.2, 0.8), 0.01, 0.5)
  flat <- within_cell_gamma_slope(w)
  expect_true(flat$lower < 0 && flat$upper > 0)
  w$p_max <- ifelse(runif(6000) < plogis(1 * (w$gamma - ave(w$gamma, w$cell))), 0.01, 0.5)
  steep <- within_cell_gamma_slope(w)
  expect_equal(steep$slope, 1, tolerance = 0.15)
})
