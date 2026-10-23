# Raw rates as the primary evidence (point 7 of the Codex review, 2026-10-22).
#   1. raw_contrasts.csv: every contrast of analysis/19_shape_contrasts.R and
#      analysis/20_formal_models_layerB.R recomputed from the raw rates, with a
#      95% Monte Carlo interval that respects the design: cells of the grid and
#      archetypes use different trees or communities (unpaired), whereas S2sigma
#      and S2h share the tree of S2 (paired, from paired_contrasts.csv). The
#      contrast on the distance from 0.05 uses a parametric bootstrap. Standard
#      errors use the Jeffreys estimate (k + 0.5) / (n + 1), so that rates of 0
#      or 1 do not give intervals of zero width; the bootstrap interval of the
#      absolute distance is approximate near 0.05, where it is not smooth.
#      `agree`: the raw and model intervals both exclude zero with the same
#      sign, or both include it.
#   2. model_discrepancy.csv: rate of every condition from the beta-binomial
#      model, from the model with a random effect per condition
#      (analysis/35_mean_structure.R) and raw, with its 95% Jeffreys interval,
#      for the conditions where the two models differ by more than 0.02;
#      `bb_in_raw_ci` tells whether the beta-binomial rate lies in that interval.
# Usage: Rscript analysis/38_raw_contrasts.R

for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
root <- file.path("data", "derived", "models")
alpha <- 0.05
z <- stats::qnorm(0.975)
n_boot <- 10000
gap_threshold <- 0.02
read_model_csv <- \(...) readr::read_csv(file.path(root, ...), show_col_types = FALSE)

# raw rate of each condition, aggregated as in the models
raw_rates <- function(m) {
  read_model_csv(m, "aggregated.csv") |>
    dplyr::mutate(p = k / n, n_f = as.character(n_f), dplyr::across(dplyr::any_of("shift"), as.character))
}
# Monte Carlo variance of a rate, with the Jeffreys estimate of p
p_jeffreys <- \(p, n) (p * n + 0.5) / (n + 1)
var_jeffreys <- \(p, n) p_jeffreys(p, n) * (1 - p_jeffreys(p, n)) / n
# unpaired difference a - b with its normal Monte Carlo interval
unpaired <- function(a, b) {
  d <- dplyr::inner_join(a, b, by = setdiff(intersect(names(a), names(b)), c("p", "n")),
                         suffix = c("_a", "_b"))
  se <- sqrt(var_jeffreys(d$p_a, d$n_a) + var_jeffreys(d$p_b, d$n_b))
  dplyr::mutate(d, raw = p_a - p_b, raw_lower = raw - z * se, raw_upper = raw + z * se) |>
    dplyr::select(-p_a, -p_b, -n_a, -n_b)
}
at_cell <- \(d, cell_value) d |>
  dplyr::filter(cell == cell_value) |>
  dplyr::select(-cell, -dplyr::any_of(c("k", "g", "cl", "n_spp", "shift_f")))

# 1a. static grid: low gamma and high Colless against Yule, by test or method
shape_raw <- purrr::map(c("type_I", "power", "clade"), \(m) {
  r <- raw_rates(m)
  if (m == "clade") r <- dplyr::filter(r, shift == "2")
  yule <- at_cell(r, "g_yule__c_yule")
  dplyr::bind_rows(
    dplyr::mutate(unpaired(at_cell(r, "g_low__c_yule"), yule), contrast = "gamma low - yule"),
    dplyr::mutate(unpaired(at_cell(r, "g_yule__c_high"), yule), contrast = "colless high - yule")
  ) |>
    dplyr::mutate(model = m) |>
    dplyr::select(-dplyr::any_of("shift"))
}) |>
  dplyr::bind_rows()

# 1b. H1: distance from the nominal level, high Colless minus Yule (bootstrap)
set.seed(20261022)
rates_type_I <- raw_rates("type_I")
p_boot <- \(p, n) stats::rbinom(n_boot, n, p_jeffreys(p, n)) / n
nominal_raw <- dplyr::inner_join(at_cell(rates_type_I, "g_yule__c_high"),
                                 at_cell(rates_type_I, "g_yule__c_yule"),
                                 by = c("test", "scenario", "n_f", "data_type"),
                                 suffix = c("_h", "_y")) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    raw = abs(p_h - alpha) - abs(p_y - alpha),
    q = list(stats::quantile(abs(p_boot(p_h, n_h) - alpha) - abs(p_boot(p_y, n_y) - alpha),
                             c(0.025, 0.975))),
    raw_lower = q[[1]], raw_upper = q[[2]]) |>
  dplyr::ungroup() |>
  dplyr::transmute(test, scenario, n_f, data_type, raw, raw_lower, raw_upper, model = "type_I",
                   contrast = "|colless high - 0.05| - |yule - 0.05|")

# 1c. S2sigma and S2h minus S2 (paired by tree), at Yule Colless
scenario_raw <- read_model_csv("paired_contrasts.csv") |>
  dplyr::filter(cell %in% c("g_yule__c_yule", "g_low__c_yule")) |>
  dplyr::transmute(test, n_f = as.character(n_spp), data_type, model = "type_I",
                   contrast = paste(sub(" - S2", "", contrast), "- S2 at", sub("__c_yule", "", cell)),
                   raw = diff, raw_lower = diff - z * se_paired, raw_upper = diff + z * se_paired)

static <- dplyr::bind_rows(shape_raw, nominal_raw, scenario_raw) |>
  dplyr::inner_join(read_model_csv("contrasts.csv") |> dplyr::mutate(n_f = as.character(n_f)),
                    by = c("model", "contrast", "test", "scenario", "method", "depth", "n_f",
                           "data_type"))

# 1d. process-based: archetypes against species sorting, N = 100 (unpaired)
process <- purrr::map(c("B_type_I", "B_power"), \(m) {
  r <- dplyr::filter(raw_rates(m), n_f == "100")
  base <- dplyr::filter(r, archetype == "sorting__none") |> dplyr::select(-archetype, -k, -n_spp)
  dplyr::filter(r, archetype != "sorting__none") |>
    dplyr::select(-k, -n_spp) |>
    unpaired(base) |>
    dplyr::mutate(model = m, contrast = paste(archetype, "- sorting__none"), archetype = NULL) |>
    dplyr::inner_join(read_model_csv(m, "archetype_contrasts.csv"),
                      by = c("test", "scenario", "cell", "contrast"))
}) |>
  dplyr::bind_rows()

# -1 or 1 when the interval excludes zero on that side, 0 otherwise
side <- \(lo, hi) (lo > 0) - (hi < 0)
out <- dplyr::bind_rows(static, process) |>
  dplyr::mutate(agree = side(raw_lower, raw_upper) == side(lower, upper),
                abs_gap = abs(raw - median))
readr::write_csv(out, file.path(root, "raw_contrasts.csv"))

# 2. conditions where the beta-binomial and the random-effect model differ
discrepancy <- purrr::map(c("type_I", "power", "clade", "B_type_I", "B_power"), \(m) {
  fit <- readRDS(file.path(root, m, "brms_fit.rds"))
  re <- readRDS(file.path(root, paste0(m, "_re"), "brms_fit.rds"))
  raw <- raw_rates(m)
  stopifnot(nrow(raw) == nrow(fit$data), all(raw$k == fit$data$k))
  dplyr::mutate(raw, model = m,
                p_bb = apply(brms::posterior_epred(fit), 2, stats::median) / n,
                p_re = apply(brms::posterior_epred(re), 2, stats::median) / n)
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(gap = p_re - p_bb,
                p_lower = stats::qbeta(0.025, k + 0.5, n - k + 0.5),
                p_upper = stats::qbeta(0.975, k + 0.5, n - k + 0.5),
                bb_in_raw_ci = p_bb >= p_lower & p_bb <= p_upper) |>
  dplyr::filter(abs(gap) > gap_threshold) |>
  dplyr::arrange(dplyr::desc(abs(gap)))
readr::write_csv(discrepancy, file.path(root, "model_discrepancy.csv"))

print(dplyr::count(out, model, contrast, agree))
print(dplyr::count(discrepancy, model, bb_in_raw_ci))
