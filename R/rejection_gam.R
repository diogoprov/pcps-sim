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
#' @return A `mgcv::bam` fit.
fit_rejection_gam <- function(d, alpha = 0.05) {
  d$reject <- as.integer(d$p_max <= alpha)
  mgcv::bam(reject ~ s(gamma, k = 6) + s(colless, k = 6), family = stats::binomial(),
            data = d, discrete = TRUE)
}

#' Rejection rate along the realized gamma, with a 95% interval
#'
#' Conditional on the median Colless of the stratum. The gamma grid spans the
#' 2.5% to 97.5% quantiles of the gamma the model was fit to.
#'
#' @param fit A fit from [fit_rejection_gam()].
#' @param n_grid Number of gamma values.
#' @return A tibble with `gamma`, the rejection rate `rate` and its interval
#'   `lower`, `upper`.
rejection_curve_gamma <- function(fit, n_grid = 25) {
  m <- fit$model
  grid <- tibble::tibble(
    gamma = seq(stats::quantile(m$gamma, 0.025, names = FALSE),
                stats::quantile(m$gamma, 0.975, names = FALSE),
                length.out = n_grid),
    colless = stats::median(m$colless)
  )
  p <- lapply(stats::predict(fit, grid, se.fit = TRUE), as.numeric)
  dplyr::mutate(grid, rate = stats::plogis(p$fit),
                lower = stats::plogis(p$fit - 1.96 * p$se.fit),
                upper = stats::plogis(p$fit + 1.96 * p$se.fit))
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
