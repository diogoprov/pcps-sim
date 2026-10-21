# Lack of fit of the mean structure of the formal models (reviewer check,
# 2026-10-20). For each model of analysis/18_formal_models.R and
# analysis/20_formal_models_layerB.R:
#   1. Pearson lack-of-fit test of the fixed structure against the binomial
#      (each aggregated row is one condition, so there is no replicate-level
#      information on extra-binomial variation: the dispersion of 2-9 is lack
#      of fit of the means);
#   2. (part "re") the same fixed structure with a binomial likelihood and a
#      random effect per condition, (1 | condition), in brms; the predicted
#      rate of every condition is compared with that of the beta-binomial
#      model (median and maximum absolute difference, and correlation).
# Writes data/derived/models/mean_structure.csv (part 1) and
# data/derived/models/<model>_re/ (part 2: fit, diagnostics, comparison).
# Usage: Rscript analysis/35_mean_structure.R <glm|re> [models]

args <- commandArgs(trailingOnly = TRUE)
part <- match.arg(args[[1]], c("glm", "re"))
models <- if (length(args) >= 2) args[-1] else c("type_I", "power", "B_type_I", "B_power", "clade")
for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
root <- file.path("data", "derived", "models")
iter <- as.integer(Sys.getenv("PCPS_RE_ITER", "2000"))  # a small value for a quick check
condition_of <- function(agg) {
  keys <- setdiff(names(agg), c("k", "n", "g", "cl", "n_f", "shift_f"))
  interaction(agg[keys], drop = TRUE)
}

if (part == "glm") {
  out <- purrr::map(models, \(m) {
    fit <- readRDS(file.path(root, m, "brms_fit.rds"))
    agg <- fit$data
    main <- stats::glm(stats::update(stats::formula(fit)$formula, cbind(k, n - k) ~ .),
                       family = stats::binomial(), data = agg)
    pearson <- sum(stats::residuals(main, type = "pearson")^2)
    tibble::tibble(model = m, conditions = nrow(agg), df = main$df.residual,
                   pearson_chi2 = pearson, dispersion = pearson / main$df.residual,
                   p_lack_of_fit = stats::pchisq(pearson, main$df.residual, lower.tail = FALSE),
                   share_resid_above_2 = mean(abs(stats::residuals(main, type = "pearson")) > 2))
  }) |>
    dplyr::bind_rows()
  readr::write_csv(out, file.path(root, "mean_structure.csv"))
  print(out)
} else {
  for (m in models) {
    fit <- readRDS(file.path(root, m, "brms_fit.rds"))
    agg <- fit$data
    raw <- readr::read_csv(file.path(root, m, "aggregated.csv"), show_col_types = FALSE)
    stopifnot(nrow(raw) == nrow(agg), all(raw$k == agg$k))  # same rows, same order
    agg$condition <- condition_of(raw)
    fixed <- as.character(stats::formula(fit)$formula)[3]
    out_dir <- file.path(root, paste0(m, "_re"))
    dir.create(out_dir, showWarnings = FALSE)
    t0 <- Sys.time()
    re <- brms::brm(
      stats::as.formula(paste("k | trials(n) ~", fixed, "+ (1 | condition)")),
      data = agg, family = stats::binomial(),
      prior = c(subset(model_priors(if (grepl("type_I", m)) -2.9 else 0), class != "phi"),
                brms::prior(exponential(2), class = sd)),
      chains = 4, cores = 4, iter = iter, warmup = iter %/% 2, seed = 20261020,
      backend = "cmdstanr", control = list(adapt_delta = 0.95, max_treedepth = 12),
      file = file.path(out_dir, "brms_fit"), file_refit = "always", refresh = 0
    )
    secs <- as.numeric(Sys.time() - t0, units = "secs")
    diag <- dplyr::mutate(fit_diagnostics(re, 12), seconds = secs)
    readr::write_csv(diag, file.path(out_dir, "diagnostics.csv"))
    p_main <- apply(brms::posterior_epred(fit), 2, stats::median) / agg$n
    p_re <- apply(brms::posterior_epred(re), 2, stats::median) / agg$n
    cmp <- tibble::tibble(model = m, median_abs_diff = stats::median(abs(p_re - p_main)),
                          max_abs_diff = max(abs(p_re - p_main)), cor = stats::cor(p_re, p_main),
                          sd_condition = stats::median(brms::as_draws_df(re)$sd_condition__Intercept))
    readr::write_csv(cmp, file.path(out_dir, "comparison.csv"))
    print(dplyr::bind_cols(cmp, diag))
  }
}
