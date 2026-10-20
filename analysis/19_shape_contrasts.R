# Posterior contrasts from the brms models of analysis/18_formal_models.R, for
# every N and data type (clade: shift 2; abundance only at N = 100, where it
# was simulated). Median and 95% credible interval of:
#   shape    - rate at low gamma minus Yule (at Yule Colless), and at high
#              Colless minus Yule (at Yule gamma), by test x scenario or
#              method x depth (H3, H4)
#   nominal  - |type I - 0.05| at high Colless minus at Yule Colless (H1: does
#              imbalance move the type I error away from the nominal level)
#   scenario - type I in S2sigma and in S2h minus S2 (the same gradient without
#              signal in sigma or h), at Yule and at low gamma (H2)
# Usage: Rscript analysis/19_shape_contrasts.R

for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
root <- file.path("data", "derived", "models")
alpha <- 0.05

summarise_draws_by <- function(d, rows, ...) {
  q <- apply(d, 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
  dplyr::mutate(rows, ..., median = q[2, ], lower = q[1, ], upper = q[3, ])
}

# grid of every level of `grp` x N x data type present in the data (clade:
# shift 2)
base_grid <- function(fit, grp, ...) {
  levs <- lapply(stats::setNames(grp, grp), \(v) levels(fit$data[[v]]))
  present <- dplyr::distinct(fit$data, n_f, data_type) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  out <- tidyr::crossing(!!!levs, present, n = 1, ...)
  if ("shift_f" %in% names(fit$data)) out$shift_f <- "2"
  out
}
rate_at <- function(fit, grid, g, cl) {
  brms::posterior_epred(fit, newdata = dplyr::mutate(grid, g = g, cl = cl))
}

shape <- purrr::map(c("type_I", "power", "clade"), \(m) {
  fit <- readRDS(file.path(root, m, "brms_fit.rds"))
  grid <- base_grid(fit, if (m == "clade") c("method", "depth") else c("test", "scenario"))
  rows <- dplyr::select(grid, -n, -dplyr::any_of("shift_f"))
  yule <- rate_at(fit, grid, "g_yule", "c_yule")
  dplyr::bind_rows(
    summarise_draws_by(rate_at(fit, grid, "g_low", "c_yule") - yule, rows, model = m,
                       contrast = "gamma low - yule"),
    summarise_draws_by(rate_at(fit, grid, "g_yule", "c_high") - yule, rows, model = m,
                       contrast = "colless high - yule")
  )
}) |>
  dplyr::bind_rows()

fit <- readRDS(file.path(root, "type_I", "brms_fit.rds"))
grid <- base_grid(fit, c("test", "scenario"))
nominal <- summarise_draws_by(
  abs(rate_at(fit, grid, "g_yule", "c_high") - alpha) -
    abs(rate_at(fit, grid, "g_yule", "c_yule") - alpha),
  dplyr::select(grid, test, scenario, n_f, data_type), model = "type_I",
  contrast = "|colless high - 0.05| - |yule - 0.05|")

scenario <- purrr::map(c("g_yule", "g_low"), \(g) {
  at <- \(sc) rate_at(fit, base_grid(fit, "test", scenario = sc), g, "c_yule")
  rows <- dplyr::select(base_grid(fit, "test", scenario = "S2"), test, n_f, data_type)
  s2 <- at("S2")
  dplyr::bind_rows(
    summarise_draws_by(at("S2sigma") - s2, rows, model = "type_I",
                       contrast = paste("S2sigma - S2 at", g)),
    summarise_draws_by(at("S2h") - s2, rows, model = "type_I",
                       contrast = paste("S2h - S2 at", g))
  )
}) |>
  dplyr::bind_rows()

readr::write_csv(dplyr::bind_rows(shape, nominal, scenario), file.path(root, "contrasts.csv"))
