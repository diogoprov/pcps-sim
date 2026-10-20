# Formal models of the simulation results (docs/PLANO_SIMULACAO.md, section 6):
# binomial models of the rejection or recovery rate, aggregated by design
# condition (1000 replicates each), with the gamma and Colless levels of the
# grid cell as factors (analysis/18_formal_models.R explains why).

#' Successes and trials by design condition
#'
#' @param d one row per replicate and test (or method), with a logical or 0/1
#'   column `success`.
#' @param by grouping columns (the design condition).
#' @return tibble with `by`, `k` (successes) and `n` (trials).
aggregate_binomial <- function(d, by) {
  d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(k = sum(success), n = dplyr::n(), .groups = "drop")
}

#' Priors on the logit scale (declared explicitly, docs/PLANO_SIMULACAO.md)
#'
#' Beta-binomial: phi is the precision of the replicate-level probabilities.
#' With 1000 trials per condition, a quasibinomial dispersion D means
#' phi ~ 999 / (D - 1) - 1; the fits give D = 3 to 14, so phi ~ 75 to 500, and
#' gamma(2, 0.01) (mean 200, sd 141) is weakly informative on that scale.
#' @param intercept_mean prior mean of the intercept (logit(0.05) ~ -2.9 for
#'   type I error, 0 for power and recovery).
model_priors <- function(intercept_mean = 0) {
  c(brms::set_prior(sprintf("normal(%g, 1.5)", intercept_mean), class = "Intercept"),
    brms::set_prior("normal(0, 1.5)", class = "b"),
    brms::set_prior("gamma(2, 0.01)", class = "phi"))
}

#' Convergence diagnostics of a brms fit
#'
#' @param fit a `brmsfit`.
#' @param max_treedepth the `max_treedepth` the model was sampled with.
#' @return one-row tibble: maximum R-hat, minimum bulk and tail ESS, number of
#'   divergent transitions and of transitions at the maximum tree depth.
fit_diagnostics <- function(fit, max_treedepth = 10) {
  draws <- posterior::summarise_draws(brms::as_draws_df(fit), "rhat", "ess_bulk",
                                      "ess_tail")
  np <- brms::nuts_params(fit)
  tibble::tibble(rhat_max = max(draws$rhat, na.rm = TRUE),
                 ess_bulk_min = min(draws$ess_bulk, na.rm = TRUE),
                 ess_tail_min = min(draws$ess_tail, na.rm = TRUE),
                 divergent = sum(np$Value[np$Parameter == "divergent__"]),
                 max_treedepth = sum(np$Value[np$Parameter == "treedepth__"] >=
                                     max_treedepth))
}

#' Posterior predictive check: coverage of the observed successes
#'
#' In-sample (the conditions the model was fit to), so the coverage is
#' optimistic; with `groups`, the successes are summed within each group first,
#' which checks the structure shared by the conditions of a group.
#' @param k observed successes, one per condition.
#' @param yrep draws x conditions matrix from `brms::posterior_predict()`.
#' @param groups optional factor, one value per condition.
#' @return one-row tibble with the fraction of conditions (or groups) whose
#'   observed successes fall in the 50% and 95% posterior predictive intervals.
ppc_coverage <- function(k, yrep, groups = NULL) {
  if (!is.null(groups)) {
    k <- tapply(k, groups, sum)
    yrep <- t(apply(yrep, 1, \(y) tapply(y, groups, sum)))
  }
  inside <- function(level) {
    q <- apply(yrep, 2, stats::quantile, probs = c((1 - level) / 2, (1 + level) / 2))
    mean(k >= q[1, ] & k <= q[2, ])
  }
  tibble::tibble(cover_50 = inside(0.5), cover_95 = inside(0.95))
}

#' Observed rate and posterior predictive rate (median, 95%) of every condition
#'
#' In-sample: the same conditions the model was fit to, so the agreement is
#' optimistic.
#' @param k,n successes and trials, one per condition.
#' @param yrep draws x conditions matrix from `brms::posterior_predict()`.
#' @return tibble with `observed`, `pred_lo`, `pred_mid` and `pred_hi`.
ppc_table <- function(k, n, yrep) {
  q <- apply(sweep(yrep, 2, n, "/"), 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
  tibble::tibble(observed = k / n, pred_lo = q[1, ], pred_mid = q[2, ], pred_hi = q[3, ])
}
