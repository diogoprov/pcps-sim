#' Binomial GAM of rejection on the realized shape of the analysed tree
#'
#' Layer B cannot fix the shape of the analysed tree: the cell sets the pool
#' tree, and the dynamics prune it. The rejection is therefore modelled as a
#' smooth function of the realized gamma and Colless, with the archetype and N
#' as fixed effects. Fit separately for each test x scenario.
#'
#' @param d Rows of one test x scenario, with `p_max`, `gamma`, `colless`,
#'   `archetype` and `n_spp`.
#' @param alpha Nominal level of the test.
#' @param k Basis dimension of each smooth.
#' @return A `mgcv::bam` fit.
fit_rejection_gam <- function(d, alpha = 0.05, k = 6) {
  d <- prepare_gam_data(d, alpha)
  # a single archetype or N cannot enter as a factor
  terms <- c(sprintf("s(gamma, k = %d)", k), sprintf("s(colless, k = %d)", k),
             if (nlevels(d$archetype) > 1) "archetype",
             if (nlevels(d$n_spp) > 1) "n_spp")
  mgcv::bam(stats::reformulate(terms, "reject"), family = stats::binomial(),
            data = d, discrete = TRUE)
}

prepare_gam_data <- function(d, alpha = 0.05) {
  d |>
    dplyr::mutate(reject = as.integer(p_max <= alpha),
                  archetype = droplevels(factor(archetype)),
                  n_spp = droplevels(factor(n_spp)))
}

#' Marginal rejection rate along the realized gamma
#'
#' At each gamma on the grid, the gamma of every replicate is set to that
#' value and the predicted probabilities are averaged, so the other terms keep
#' their observed distribution (marginal standardization).
#'
#' @param fit A fit from [fit_rejection_gam()].
#' @param d The rows the model was fit to.
#' @param n_grid Number of gamma values, spread between the 2.5% and 97.5%
#'   quantiles of the observed gamma.
#' @param max_rows Averages over a random subset of at most this many rows, to
#'   keep the prediction cheap.
#' @return A tibble with `gamma` and the marginal rejection rate `rate`.
marginal_rejection_gamma <- function(fit, d, n_grid = 25, max_rows = 5000) {
  d <- prepare_gam_data(d)
  if (nrow(d) > max_rows) d <- d[sample.int(nrow(d), max_rows), ]
  grid <- seq(stats::quantile(d$gamma, 0.025), stats::quantile(d$gamma, 0.975),
              length.out = n_grid)
  rate <- vapply(grid, \(g) {
    mean(stats::predict(fit, transform(d, gamma = g), type = "response"))
  }, numeric(1))
  tibble::tibble(gamma = grid, rate = rate)
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
