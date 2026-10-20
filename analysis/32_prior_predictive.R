# Prior predictive checks of the formal models: each model refitted with
# sample_prior = "only" (same formula, data and priors), and the rates its
# prior predicts for the observed conditions summarized: quantiles over all
# draws and conditions, and the share of predicted rates below 0.01 or above
# 0.99. Writes data/derived/models/prior_predictive.csv.
# Usage: Rscript analysis/32_prior_predictive.R [models]

args <- commandArgs(trailingOnly = TRUE)
models <- if (length(args)) args else c("type_I", "power", "B_type_I", "B_power", "clade")
out <- purrr::map(models, \(m) {
  fit <- readRDS(file.path("data", "derived", "models", m, "brms_fit.rds"))
  t0 <- Sys.time()
  pr <- stats::update(fit, sample_prior = "only", chains = 4, cores = 4, iter = 1000,
                      seed = 20261020, refresh = 0, backend = "cmdstanr")
  rate <- brms::posterior_epred(pr, ndraws = 1000) / matrix(fit$data$n, 1000, nrow(fit$data), byrow = TRUE)
  q <- stats::quantile(rate, c(0.05, 0.25, 0.5, 0.75, 0.95))
  message(m, ": ", format(Sys.time() - t0))
  tibble::tibble(model = m, q05 = q[1], q25 = q[2], median = q[3], q75 = q[4], q95 = q[5],
                 below_0.01 = mean(rate < 0.01), above_0.99 = mean(rate > 0.99))
}) |>
  dplyr::bind_rows()
readr::write_csv(out, file.path("data", "derived", "models", "prior_predictive.csv"))
print(out)
