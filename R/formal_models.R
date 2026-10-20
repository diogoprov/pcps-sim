# Formal models of the simulation results (docs/PLANO_SIMULACAO.md, section 6):
# binomial models of the rejection or recovery rate, aggregated by design
# condition (1000 replicates each), with the tree shape of the cell as
# continuous covariates. mgcv to explore, brms for the final model.

#' Shape targets of every grid cell and N, standardized within N
#'
#' The gamma target scales with N, so it is standardized within each N; the
#' Colless target (normalized against Yule) is comparable across N.
#' @param grid tree grid with `n`, `cell`, `gamma_target`, `colless_target`.
cell_shapes <- function(grid) {
  grid |>
    dplyr::transmute(n_spp = n, cell, gamma_target = unname(gamma_target),
                     colless_target = unname(colless_target)) |>
    dplyr::group_by(n_spp) |>
    dplyr::mutate(gamma_z = as.numeric(scale(gamma_target))) |>
    dplyr::ungroup() |>
    dplyr::mutate(colless_z = as.numeric(scale(colless_target)))
}

#' Successes and trials by design condition
#'
#' @param d one row per replicate and test (or method), with a logical or 0/1
#'   column `success`.
#' @param by grouping columns (the design condition).
#' @return tibble with `by`, `k` (successes), `n` (trials) and `obs` (one id
#'   per row, for an observation-level random effect).
aggregate_binomial <- function(d, by) {
  d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(k = sum(success), n = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(obs = factor(dplyr::row_number()))
}

#' Priors on the logit scale (declared explicitly, docs/PLANO_SIMULACAO.md)
#'
#' @param intercept_mean prior mean of the intercept (logit(0.05) ~ -2.9 for
#'   type I error, 0 for power and recovery).
model_priors <- function(intercept_mean = 0) {
  c(brms::set_prior(sprintf("normal(%g, 1.5)", intercept_mean), class = "Intercept"),
    brms::set_prior("normal(0, 1.5)", class = "b"),
    brms::set_prior("exponential(1)", class = "sd"))
}

#' Convergence diagnostics of a brms fit
#'
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
#' The predictions draw new levels of the observation-level random effect, so
#' they integrate over it; with the fitted levels every condition would be
#' reproduced almost exactly and the check would be trivially passed.
#' @return one-row tibble with the fraction of conditions whose observed `k`
#'   falls in the 50% and 95% posterior predictive intervals.
ppc_coverage <- function(fit) {
  new <- fit$data
  new$obs <- paste0("new_", new$obs)
  yrep <- brms::posterior_predict(fit, newdata = new, allow_new_levels = TRUE,
                                  sample_new_levels = "gaussian")
  k <- fit$data$k
  inside <- function(level) {
    q <- apply(yrep, 2, stats::quantile, probs = c((1 - level) / 2, (1 + level) / 2))
    mean(k >= q[1, ] & k <= q[2, ])
  }
  tibble::tibble(cover_50 = inside(0.5), cover_95 = inside(0.95))
}
