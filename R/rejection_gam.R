#' Binomial GAM of rejection on the realized shape of the analysed tree
#'
#' Layer B cannot fix the shape of the analysed tree: the cell sets the pool
#' tree, and the dynamics prune it. The rejection is therefore modelled as a
#' smooth function of the realized gamma and Colless. Fit within one stratum
#' (test x scenario x archetype x N): the realized gamma scales with the number
#' of tips and differs among archetypes, so pooling strata would mix rates and
#' extrapolate gamma beyond the support of each stratum.
#'
#' @param d Rows of one stratum, with `p_max`, `gamma` and `colless`.
#' @param alpha Nominal level of the test.
#' @return A `mgcv::gam` fit.
fit_rejection_gam <- function(d, alpha = 0.05) {
  d$reject <- as.integer(d$p_max <= alpha)
  mgcv::gam(reject ~ s(gamma, k = 6) + s(colless, k = 6), family = stats::binomial(),
            data = d, method = "REML")
}

#' Rejection rate along the realized gamma, with a 95% interval
#'
#' Conditional on the median Colless of the stratum. The gamma grid spans the
#' 2.5% to 97.5% quantiles of the gamma the model was fit to. The realized
#' gamma clumps around the pool cells, so grid points farther than `max_dist`
#' from any observed gamma get `NA` instead of an interpolated rate.
#'
#' @param fit A fit from [fit_rejection_gam()].
#' @param n_grid Number of gamma values.
#' @param max_dist Largest distance to the nearest observed gamma.
#' @return A tibble with `gamma`, the rejection rate `rate` and its interval
#'   `lower`, `upper`.
rejection_curve_gamma <- function(fit, n_grid = 25, max_dist = 0.25) {
  m <- fit$model
  grid <- tibble::tibble(
    gamma = seq(stats::quantile(m$gamma, 0.025, names = FALSE),
                stats::quantile(m$gamma, 0.975, names = FALSE),
                length.out = n_grid),
    colless = stats::median(m$colless)
  )
  p <- lapply(stats::predict(fit, grid, se.fit = TRUE), as.numeric)
  observed <- sort(m$gamma)
  supported <- vapply(grid$gamma, \(g) min(abs(observed - g)) <= max_dist, logical(1))
  p$fit[!supported] <- NA
  z <- stats::qnorm(0.975)
  dplyr::mutate(grid, rate = stats::plogis(p$fit),
                lower = stats::plogis(p$fit - z * p$se.fit),
                upper = stats::plogis(p$fit + z * p$se.fit))
}

#' Effect of the realized gamma within the pool cells
#'
#' Logistic regression of rejection on the pool cell plus the realized gamma
#' and Colless. With the cell as a factor, the gamma slope only uses the
#' spread of the realized gamma around each pool tree, which comes from the
#' pruning, so it separates the analysed tree from the pool.
#'
#' @param d Rows of one stratum, with `p_max`, `cell`, `gamma` and `colless`.
#' @param alpha Nominal level of the test.
#' @return A one-row tibble with the slope of `gamma` on the logit scale, its
#'   95% Wald interval `lower`, `upper`, and `p_value`.
within_cell_gamma_slope <- function(d, alpha = 0.05) {
  d$reject <- as.integer(d$p_max <= alpha)
  fit <- stats::glm(reject ~ cell + gamma + colless, family = stats::binomial(), data = d)
  b <- summary(fit)$coefficients["gamma", ]
  z <- stats::qnorm(0.975)
  tibble::tibble(slope = b[["Estimate"]], lower = b[["Estimate"]] - z * b[["Std. Error"]],
                 upper = b[["Estimate"]] + z * b[["Std. Error"]], p_value = b[["Pr(>|z|)"]])
}

#' Summary of the smooth terms of a rejection GAM
#'
#' @param fit A fit from [fit_rejection_gam()].
#' @return A tibble with one row per smooth: `term`, `edf`, `chi_sq`, `p_value`,
#'   plus the deviance explained of the whole model.
summarise_rejection_gam <- function(fit) {
  s <- summary(fit)
  tibble::tibble(term = rownames(s$s.table), edf = s$s.table[, "edf"],
                 chi_sq = s$s.table[, "Chi.sq"], p_value = s$s.table[, "p-value"],
                 dev_expl = s$dev.expl)
}
