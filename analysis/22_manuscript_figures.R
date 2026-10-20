# Figures of the manuscript (manuscript/figures/), from the formal models
# (analysis/18_formal_models.R, analysis/20_formal_models_layerB.R) and the
# regional pools (analysis/21_regional_pools.R). Rates are posterior medians at
# N = 100 and incidence data (clade recovery: shift 2).
# Usage: Rscript analysis/22_manuscript_figures.R

models <- file.path("data", "derived", "models")
out_dir <- file.path("manuscript", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
test_labels <- c(adonis = "PERMANOVA", mantel = "Mantel", glm_pcps1 = "GLM PCPS1",
                 glm_pcps2 = "GLM PCPS2", rda_pcps = "RDA PCPS", pcc = "PCC")
method_labels <- c(pcps = "PCPS", pcc = "PCC", evopca = "evoPCA",
                   phylofactor = "Phylofactorization")
level_labels <- c(g_yule = "Yule", g_mid = "mid", g_low = "low",
                  c_yule = "Yule", c_mid = "mid", c_high = "high")
scenario_label <- function(x) sub("sigma", "\u03c3", x)
theme_ms <- ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 strip.background = ggplot2::element_rect(fill = "grey95"))

# rate over the gamma x Colless grid, one panel per group x row
surface_plot <- function(d, group, row, labels, fill_limits, fill_name) {
  d <- d |>
    dplyr::mutate(g = factor(level_labels[g], c("low", "mid", "Yule")),
                  cl = factor(level_labels[cl], c("Yule", "mid", "high")),
                  group = factor(labels[.data[[group]]], labels),
                  row = scenario_label(.data[[row]]))
  ggplot2::ggplot(d, ggplot2::aes(cl, g, fill = rate)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rate),
                                    colour = rate > fill_limits[1] + 0.5 * diff(fill_limits)),
                       size = 2.2, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white")) +
    ggplot2::scale_fill_distiller(palette = "YlGnBu", direction = 1, limits = fill_limits,
                                  oob = scales::squish, name = fill_name) +
    ggplot2::facet_grid(row ~ group) +
    ggplot2::labs(x = "Colless level", y = expression(gamma ~ "level")) +
    theme_ms
}

type1 <- readr::read_csv(file.path(models, "type_I", "shape_surface.csv"), show_col_types = FALSE)
p <- surface_plot(dplyr::filter(type1, test %in% names(test_labels)), "test", "scenario",
                  test_labels, c(0, 0.17), "Type I\nerror")
ggplot2::ggsave(file.path(out_dir, "fig2_type1.png"), p, width = 17, height = 14, units = "cm",
                dpi = 300)

power <- readr::read_csv(file.path(models, "power", "shape_surface.csv"), show_col_types = FALSE)
p <- surface_plot(dplyr::filter(power, test %in% names(test_labels)), "test", "scenario",
                  test_labels, c(0, 1), "Power")
ggplot2::ggsave(file.path(out_dir, "fig3_power.png"), p, width = 17, height = 9, units = "cm",
                dpi = 300)

clade <- readr::read_csv(file.path(models, "clade", "shape_surface.csv"), show_col_types = FALSE)
p <- surface_plot(clade, "method", "depth", method_labels, c(0, 1), "Recovery")
ggplot2::ggsave(file.path(out_dir, "fig5_clade.png"), p, width = 15, height = 8, units = "cm",
                dpi = 300)

# layer B: power by archetype and pool cell
archetype_labels <- c(sorting__none = "Sorting", sorting__phylo = "Sorting + phylo. competition",
                      mass_effects__none = "Mass effects",
                      mass_effects__phylo = "Mass effects + phylo. competition")
cell_labels <- c(g_yule__c_yule = "Yule/Yule", g_low__c_yule = "low/Yule",
                 g_yule__c_high = "Yule/high", g_low__c_high = "low/high",
                 g_mid__c_mid = "mid/mid")
b_power <- readr::read_csv(file.path(models, "B_power", "rates.csv"), show_col_types = FALSE) |>
  dplyr::mutate(test = factor(test_labels[test], test_labels),
                archetype = factor(archetype_labels[archetype], archetype_labels),
                cell = factor(cell_labels[cell], cell_labels),
                scenario = scenario_label(scenario))
p <- ggplot2::ggplot(b_power, ggplot2::aes(cell, median, colour = archetype)) +
  ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), size = 0.15,
                           position = ggplot2::position_dodge(width = 0.6)) +
  ggplot2::facet_grid(scenario ~ test) +
  ggplot2::scale_colour_manual(values = c("#0072B2", "#56B4E9", "#D55E00", "#CC79A7"),
                               name = NULL) +
  ggplot2::labs(x = expression("Pool tree cell (" * gamma * "/Colless level)"),
                y = "Power (N = 100)") +
  theme_ms +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "fig4_layerB_power.png"), p, width = 17, height = 11,
                units = "cm", dpi = 300)

# Box 1: regional pools, empirical anuran clades and the grid targets
grid <- readRDS(file.path("data", "derived", "trees", "tree_grid_k200.rds")) |>
  dplyr::transmute(n_spp = factor(n), gamma = unname(gamma_target),
                   colless = unname(colless_target))
clades <- readr::read_csv(file.path("data", "derived", "trees", "empirical_shape_portik2023.csv"),
                          show_col_types = FALSE) |>
  dplyr::filter(source == "clade")
pools <- readr::read_csv(file.path("data", "derived", "pools", "pool_shapes.csv"),
                         show_col_types = FALSE)
p <- ggplot2::ggplot() +
  ggplot2::geom_point(data = clades, ggplot2::aes(colless, gamma), colour = "grey75", size = 0.6) +
  ggplot2::geom_point(data = grid, ggplot2::aes(colless, gamma, shape = n_spp), size = 1.8) +
  ggplot2::geom_errorbar(data = pools, ggplot2::aes(x = colless_boot_median,
                                                    ymin = gamma_boot_lo, ymax = gamma_boot_hi,
                                                    colour = pool), width = 0) +
  ggplot2::geom_errorbar(data = pools, ggplot2::aes(y = gamma_boot_median,
                                                    xmin = colless_boot_lo, xmax = colless_boot_hi,
                                                    colour = pool), width = 0,
                         orientation = "y") +
  ggplot2::geom_point(data = pools, ggplot2::aes(colless_boot_median, gamma_boot_median,
                                                 colour = pool), size = 2.4) +
  ggplot2::geom_point(data = pools, ggplot2::aes(colless, gamma, colour = pool), shape = 5,
                      size = 2.2, stroke = 0.8, show.legend = FALSE) +
  ggplot2::scale_shape_manual(values = c(1, 2, 0), name = "Grid targets, N") +
  ggplot2::scale_colour_manual(values = c("#1b9e77", "#7570b3", "#d95f02"), name = "Regional pool") +
  ggplot2::labs(x = "Colless index (normalized against Yule)", y = expression(gamma)) +
  theme_ms
ggplot2::ggsave(file.path(out_dir, "box1_pools.png"), p, width = 15, height = 10, units = "cm",
                dpi = 300)
