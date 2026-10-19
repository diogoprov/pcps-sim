# Layer B: rejection rate per pool cell, placed at the realized gamma of the
# analysed tree. The realized gamma clumps around the pool cells, so one point
# per cell (Wilson 95% interval on the rate; horizontal bar = 2.5-97.5%
# quantiles of the realized gamma) shows the data without interpolating
# between the clumps, as the GAM curves of analysis/14_gam_layerB.R would.
# Usage: Rscript analysis/15_fig_layerB_cells.R

alpha <- 0.05
tests <- c("adonis", "mantel", "glm_pcps1", "glm_pcps2", "rda_pcps", "pcc")
in_file <- file.path("data", "derived", "layerB", "results.csv")
out_dir <- file.path("data", "derived", "layerB_gam")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

res <- data.table::fread(in_file, select = c("scenario", "archetype", "n_spp", "cell",
                                             "test", "p_max", "gamma")) |>
  tibble::as_tibble() |>
  dplyr::filter(test %in% tests)

cells <- res |>
  dplyr::group_by(scenario, archetype, n_spp, cell, test) |>
  dplyr::summarise(k = sum(p_max <= alpha), n = dplyr::n(),
                   gamma_mid = stats::median(gamma),
                   gamma_lo = stats::quantile(gamma, 0.025, names = FALSE),
                   gamma_hi = stats::quantile(gamma, 0.975, names = FALSE),
                   .groups = "drop") |>
  dplyr::mutate(rate = k / n,
                ci = purrr::map2(k, n, \(k, n) stats::prop.test(k, n, correct = FALSE)$conf.int),
                lower = purrr::map_dbl(ci, 1), upper = purrr::map_dbl(ci, 2),
                colless = sub(".*__c_", "", cell),
                test = factor(test, tests),
                # side by side within a cell, so the tests do not overlap
                x = gamma_mid + (as.integer(test) - (length(tests) + 1) / 2) * 0.12) |>
  dplyr::select(-ci)
readr::write_csv(cells, file.path(out_dir, "cell_rates.csv"))

for (n_tips in sort(unique(cells$n_spp))) {
  p <- cells |>
    dplyr::filter(n_spp == n_tips) |>
    ggplot2::ggplot(ggplot2::aes(x, rate, colour = test, shape = colless)) +
    ggplot2::geom_hline(yintercept = alpha, linetype = 2, colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = gamma_lo, xmax = gamma_hi), width = 0,
                           orientation = "y", alpha = 0.2) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), size = 0.25) +
    ggplot2::facet_grid(archetype ~ scenario, scales = "free_y") +
    ggplot2::labs(x = expression("Realized " * gamma * " of the analysed tree (median, 95% range)"),
                  y = sprintf("Rejection rate (N = %d)", n_tips), colour = "Test",
                  shape = "Colless of the pool cell") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(file.path(out_dir, sprintf("cell_rates_n%d.png", n_tips)), p,
                  width = 16, height = 9, dpi = 150)
}
