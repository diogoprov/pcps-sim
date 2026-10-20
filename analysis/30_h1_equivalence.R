# Equivalence test of H1 on the type I model of Layer A: for every test, null
# scenario, N and data type, the posterior probability that the type I error at
# high imbalance differs from that at Yule imbalance by less than 0.025 (half
# the width of Bradley's liberal interval, [0.025, 0.075]), at Yule and at low
# gamma. Writes data/derived/models/h1_equivalence.csv.
# Usage: Rscript analysis/30_h1_equivalence.R

margin <- 0.025
fit <- readRDS(file.path("data", "derived", "models", "type_I", "brms_fit.rds"))
levs <- list(test = levels(fit$data$test), scenario = levels(fit$data$scenario))
present <- dplyr::distinct(fit$data, n_f, data_type) |>
  dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
grid <- tidyr::crossing(!!!levs, present, n = 1)
rate_at <- function(g, cl) brms::posterior_epred(fit, newdata = dplyr::mutate(grid, g = g, cl = cl))
out <- purrr::map(c("g_yule", "g_low"), \(g) {
  d <- rate_at(g, "c_high") - rate_at(g, "c_yule")
  dplyr::mutate(dplyr::select(grid, test, scenario, n_f, data_type), gamma = g,
                median = apply(d, 2, stats::median), p_equivalent = colMeans(abs(d) < margin))
}) |>
  dplyr::bind_rows()
readr::write_csv(out, file.path("data", "derived", "models", "h1_equivalence.csv"))
print(summary(out$p_equivalent))
print(dplyr::summarise(out, conditions = dplyr::n(), above_0.95 = sum(p_equivalent > 0.95),
                       min = min(p_equivalent)))
print(dplyr::arrange(out, p_equivalent) |> utils::head(5))
