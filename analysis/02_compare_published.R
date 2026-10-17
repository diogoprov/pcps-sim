# Compare layer-A validation runs with the published tables.
# Flags cells whose published rate falls outside the 95% binomial interval of ours.
# Usage: Rscript analysis/02_compare_published.R <n_reps>

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[[1]]) else 1000L
for (f in list.files("R", full.names = TRUE)) source(f)
val <- file.path("data", "derived", "validation")

binom_flag <- function(ours, pub, n) {
  ci <- vapply(seq_along(ours), \(i)
    stats::binom.test(round(ours[i] * n[i]), n[i])$conf.int, numeric(2))
  pub >= ci[1, ] & pub <= ci[2, ]
}

compare <- function(ours, pub, keys) {
  ours |>
    tidyr::pivot_longer(c(p_site, p_taxa, p_max), names_to = "rate", values_to = "ours") |>
    dplyr::inner_join(
      tidyr::pivot_longer(pub, c(p_site, p_taxa, p_max), names_to = "rate",
                          values_to = "published"),
      by = c(keys, "rate")
    ) |>
    dplyr::mutate(in_ci = binom_flag(ours, published, n))
}

res_d <- readr::read_csv(file.path(val, sprintf("duarte2016_reps%d.csv", n_reps)),
                         show_col_types = FALSE)
pub_d <- readr::read_csv("data/published/duarte2016_table1.csv", show_col_types = FALSE)
cmp_d <- compare(rejection_rates(res_d, test, scenario, n_spp),
                 dplyr::select(pub_d, -stat), c("test", "scenario", "n_spp"))
stat_d <- res_d |>
  dplyr::summarise(ours = mean(stat), .by = c(test, scenario, n_spp)) |>
  dplyr::inner_join(dplyr::select(pub_d, test, scenario, n_spp, published = stat),
                    by = c("test", "scenario", "n_spp"))

res_p <- readr::read_csv(file.path(val, sprintf("pnk2015_reps%d.csv", n_reps)),
                         show_col_types = FALSE)
pub_p <- readr::read_csv("data/published/pnk2015_table1.csv", show_col_types = FALSE)
cmp_p <- compare(rejection_rates(res_p, squared, scenario, mu_tol), pub_p,
                 c("scenario", "mu_tol"))

readr::write_csv(cmp_d, file.path(val, "compare_duarte2016.csv"))
readr::write_csv(cmp_p, file.path(val, "compare_pnk2015.csv"))
print(cmp_d, n = Inf)
print(stat_d, n = Inf)
print(cmp_p, n = Inf)
