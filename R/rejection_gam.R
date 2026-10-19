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
  d$reject <- rejected(d, alpha)
  mgcv::gam(reject ~ s(gamma, k = 6) + s(colless, k = 6), family = stats::binomial(),
            data = d, method = "REML")
}

#' Rejection rate along the realized gamma, with a 95% interval
#'
#' Conditional on the median Colless of the stratum. The gamma grid spans the
#' 2.5% to 97.5% quantiles of the gamma the model was fit to. The realized
#' gamma clumps around the pool cells, so grid points farther than `max_dist`
#' from any observed gamma get `NA` instead of an interpolated rate. The
#' default is half the tolerance of the pool tree around its cell target.
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
  supported <- vapply(grid$gamma, \(g) min(abs(m$gamma - g)) <= max_dist, logical(1))
  p$fit[!supported] <- NA
  ci <- wald_ci(p$fit, p$se.fit)
  dplyr::mutate(grid, rate = stats::plogis(p$fit),
                lower = stats::plogis(ci$lower), upper = stats::plogis(ci$upper))
}

#' Effect of the pruning on gamma, given the pool tree
#'
#' Logistic regression of rejection on the pool cell, the gamma of the pool
#' tree, the change of gamma caused by the pruning (`gamma - gamma_pool`) and
#' the realized Colless. Each replicate draws its own pool tree within the
#' tolerance of the cell, so the cell alone does not fix the pool; the slope of
#' the shift is the effect of the analysed tree beyond its pool.
#'
#' @param d Rows of one stratum, with `p_max`, `cell`, `gamma_pool`, `gamma`
#'   and `colless`.
#' @param alpha Nominal level of the test.
#' @return A one-row tibble with the slope of the shift on the logit scale, its
#'   95% Wald interval `lower`, `upper`, `p_value`, the SD of the shift and
#'   whether the fit converged.
pruning_gamma_slope <- function(d, alpha = 0.05) {
  d$reject <- rejected(d, alpha)
  d$shift <- d$gamma - d$gamma_pool
  fit <- suppressWarnings(stats::glm(reject ~ cell + gamma_pool + shift + colless,
                                     family = stats::binomial(), data = d))
  b <- stats::coef(summary(fit))
  # an aliased shift (no spread) drops out of the coefficient table
  est <- if ("shift" %in% rownames(b)) b["shift", ] else c(NA, NA, NA, NA)
  ci <- wald_ci(est[[1]], est[[2]])
  tibble::tibble(slope = est[[1]], lower = ci$lower, upper = ci$upper, p_value = est[[4]],
                 shift_sd = stats::sd(d$shift), converged = fit$converged)
}

rejected <- function(d, alpha) as.integer(d$p_max <= alpha)

wald_ci <- function(est, se, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  list(lower = est - z * se, upper = est + z * se)
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
