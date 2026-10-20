# Tables of the Supporting Information (manuscript/si/):
#   S1 - targets of the tree-shape grid (gamma, Colless normalized against
#        Yule) by N, with the calibrated (beta, c) and acceptance rates;
#   S2 - every rate (rejection rate or clade-recovery hit rate) by condition,
#        with its Monte Carlo standard error and 95% Monte Carlo interval;
#   S3 - convergence and posterior predictive diagnostics of the formal models.
# Usage: Rscript analysis/24_si_tables.R

models <- file.path("data", "derived", "models")
out_dir <- file.path("manuscript", "si")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
level_labels <- c(yule = "Yule", mid = "median", low = "10% quantile", high = "90% quantile")

grid <- readRDS(file.path("data", "derived", "trees", "tree_grid_k200.rds"))
s1 <- grid |>
  dplyr::transmute(N = n,
                   gamma_level = level_labels[sub("g_(.*)__.*", "\\1", cell)],
                   colless_level = level_labels[sub(".*__c_", "", cell)],
                   gamma_target = round(unname(gamma_target), 2),
                   colless_target = round(unname(colless_target), 2),
                   c = round(c, 3), beta = round(beta, 3),
                   acceptance = round(acceptance, 3))
readr::write_csv(s1, file.path(out_dir, "table_S1_grid.csv"))

# rates with their Monte Carlo standard error, sqrt(p (1 - p) / K)
with_mcse <- function(d, model) {
  d |>
    dplyr::mutate(model = model, rate = k / n, mcse = sqrt(rate * (1 - rate) / n),
                  mc_lower = pmax(rate - 1.96 * mcse, 0), mc_upper = pmin(rate + 1.96 * mcse, 1),
                  .before = 1)
}
s2 <- dplyr::bind_rows(
  with_mcse(readr::read_csv(file.path(models, "type_I", "aggregated.csv"), show_col_types = FALSE),
            "Layer A, type I error"),
  with_mcse(readr::read_csv(file.path(models, "power", "aggregated.csv"), show_col_types = FALSE),
            "Layer A, power"),
  with_mcse(readr::read_csv(file.path(models, "B_type_I", "aggregated.csv"), show_col_types = FALSE),
            "Layer B, type I error"),
  with_mcse(readr::read_csv(file.path(models, "B_power", "aggregated.csv"), show_col_types = FALSE),
            "Layer B, power"),
  with_mcse(readr::read_csv(file.path(models, "clade", "aggregated.csv"), show_col_types = FALSE),
            "Clade recovery")
) |>
  dplyr::select(model, dplyr::any_of(c("test", "method", "scenario", "archetype", "cell",
                                       "n_spp", "data_type", "depth", "shift")),
                k, n, rate, mcse, mc_lower, mc_upper)
readr::write_csv(s2, file.path(out_dir, "table_S2_rates_mcse.csv"))

s3 <- purrr::map(c("type_I", "power", "B_type_I", "B_power", "clade"), \(m) {
  readr::read_csv(file.path(models, m, "diagnostics.csv"), show_col_types = FALSE) |>
    dplyr::mutate(model = m, .before = 1)
}) |>
  dplyr::bind_rows()
readr::write_csv(s3, file.path(out_dir, "table_S3_diagnostics.csv"))
message("S1: ", nrow(s1), " rows; S2: ", nrow(s2), " rows; S3: ", nrow(s3), " rows")
