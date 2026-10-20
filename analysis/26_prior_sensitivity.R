# Prior sensitivity of the formal models by power-scaling (priorsense), without
# refitting. Coefficients of these models are many and weakly identified one
# by one, so the quantities that matter are the predicted rates of every
# condition (the rates and contrasts reported); both are assessed.
# Writes data/derived/models/prior_sensitivity.csv (coefficients) and
# data/derived/models/prior_sensitivity_rates.csv (one row per condition).
# Usage: Rscript analysis/26_prior_sensitivity.R

models <- c("type_I", "power", "B_type_I", "B_power", "clade")

rate_draws <- function(fit, newdata) {
  nch <- posterior::nchains(fit)
  nit <- posterior::niterations(fit)
  function(x, ...) {
    # posterior_epred returns the draws chain by chain, so they fill the array
    # iteration first
    e <- brms::posterior_epred(x, newdata = newdata) / newdata$n
    e <- t(t(e))
    a <- array(e, dim = c(nit, nch, ncol(e)),
               dimnames = list(NULL, NULL, paste0("rate[", seq_len(ncol(e)), "]")))
    posterior::as_draws_df(posterior::as_draws_array(a))
  }
}

coefs <- list()
rates <- list()
for (m in models) {
  fit <- readRDS(file.path("data", "derived", "models", m, "brms_fit.rds"))
  v <- grep("^(b_|phi$)", brms::variables(fit), value = TRUE)
  coefs[[m]] <- tibble::as_tibble(as.data.frame(priorsense::powerscale_sensitivity(fit, variable = v))) |>
    dplyr::mutate(model = m, .before = 1)
  conditions <- unique(fit$data[, setdiff(names(fit$data), "k"), drop = FALSE])
  ps <- priorsense::powerscale_sensitivity(fit, prediction = rate_draws(fit, conditions))
  d <- tibble::as_tibble(as.data.frame(ps))
  d <- d[grepl("^rate\\[", d$variable), ]
  rates[[m]] <- dplyr::bind_cols(tibble::tibble(model = m),
                                 conditions[as.integer(gsub("\\D", "", d$variable)), , drop = FALSE],
                                 dplyr::select(d, prior, likelihood, diagnosis))
}
coefs <- dplyr::bind_rows(coefs)
rates <- dplyr::bind_rows(rates)
readr::write_csv(coefs, file.path("data", "derived", "models", "prior_sensitivity.csv"))
readr::write_csv(rates, file.path("data", "derived", "models", "prior_sensitivity_rates.csv"))
summ <- function(d) dplyr::summarise(dplyr::group_by(d, model), n = dplyr::n(),
                                     flagged = sum(diagnosis != "-"),
                                     median_prior = stats::median(prior), max_prior = max(prior))
print(summ(coefs))
print(summ(rates))
