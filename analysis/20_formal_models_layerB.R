# Formal models of layer B (process-based metacommunities): beta-binomial
# models of the rejection rate aggregated by condition (1000 replicates each),
# for type I error (scenarios under H0) and power (under H1). The question is
# how the tests fare under different assembly processes (archetypes: species
# sorting or mass effects, with or without phylogenetic competition) and pool
# tree shapes. Layer B has five pool cells (the four corners of the grid and
# the centre), so the cell enters as a 5-level factor; its effect may differ
# among tests and archetypes, and the effects of archetype, scenario and cell
# may differ with N (mass effects with phylogenetic competition has no N = 200,
# so that combination rests on the other interactions). The cell is the shape
# of the pool tree: under phylogenetic competition the pruning changes the
# analysed tree (gamma more negative), so the contrast of competition mixes the
# process with a change of tree shape. Same priors, sampler settings and checks
# as analysis/18_formal_models.R (fit_formal_model() in R/formal_models.R).
# Outputs in <out_root>/B_<type_I|power>/.
# Usage: Rscript analysis/20_formal_models_layerB.R <type_I|power> [chains] [iter] [out_root]

args <- commandArgs(trailingOnly = TRUE)
model <- match.arg(args[[1]], c("type_I", "power"))
chains <- if (length(args) >= 2) as.integer(args[[2]]) else 4L
iter <- if (length(args) >= 3) as.integer(args[[3]]) else 2000L
out_root <- if (length(args) >= 4) args[[4]] else file.path("data", "derived", "models")

for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
out_dir <- file.path(out_root, paste0("B_", model))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
alpha <- 0.05
tests <- c("adonis", "mantel", "glm_pcps1", "glm_pcps2", "rda_pcps", "pcc")
scenarios <- if (model == "type_I") c("S0", "S1", "S2", "S2sigma", "S2h") else
  c("S3a", "S3b", "S3sigma")
cells <- c("g_yule__c_yule", "g_low__c_yule", "g_yule__c_high", "g_low__c_high",
           "g_mid__c_mid")

agg <- data.table::fread(file.path("data", "derived", "layerB", "results.csv"),
                         select = c("scenario", "archetype", "cell", "n_spp", "test",
                                    "p_max")) |>
  tibble::as_tibble() |>
  dplyr::filter(scenario %in% scenarios, test %in% tests) |>
  dplyr::mutate(success = p_max <= alpha)
stopifnot(!anyNA(agg$success))
agg <- aggregate_binomial(agg, c("test", "scenario", "archetype", "cell", "n_spp")) |>
  dplyr::mutate(test = factor(test, tests), scenario = factor(scenario, scenarios),
                archetype = factor(archetype), cell = factor(cell, cells),
                n_f = factor(n_spp))
stopifnot(!anyNA(agg$cell), length(unique(agg$n)) == 1)
readr::write_csv(agg, file.path(out_dir, "aggregated.csv"))
message("B_", model, ": ", nrow(agg), " conditions")

fixed <- "test * archetype * (scenario + cell) + test * n_f * (archetype + scenario + cell)"
fit <- fit_formal_model(fixed, agg,
                        groups = interaction(agg$test, agg$archetype, agg$n_f, agg$cell,
                                             drop = TRUE),
                        out_dir = out_dir, intercept_mean = if (model == "type_I") -2.9 else 0,
                        chains = chains, iter = iter,
                        # with the N interactions the type I model hit tree depth 10
                        # in ~1% of the transitions
                        max_treedepth = 12)

# Rates by test x archetype x scenario x cell at N = 100 (median, 95% CI), and
# the contrast of each archetype with species sorting without competition
grid <- tidyr::crossing(test = levels(agg$test), archetype = levels(agg$archetype),
                        scenario = levels(agg$scenario), cell = levels(agg$cell),
                        n_f = "100", n = 1)
ep <- brms::posterior_epred(fit, newdata = grid)
q <- apply(ep, 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
grid <- dplyr::mutate(grid, median = q[2, ], lower = q[1, ], upper = q[3, ])
readr::write_csv(grid, file.path(out_dir, "rates.csv"))
ref <- grid$archetype == "sorting__none"
key <- c("test", "scenario", "cell")
contr <- purrr::map(setdiff(levels(agg$archetype), "sorting__none"), \(a) {
  this <- grid$archetype == a
  i_ref <- match(interaction(grid[this, key]), interaction(grid[ref, key]))
  d <- ep[, this, drop = FALSE] - ep[, which(ref)[i_ref], drop = FALSE]
  qq <- apply(d, 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
  dplyr::mutate(grid[this, key], contrast = paste(a, "- sorting__none"),
                median = qq[2, ], lower = qq[1, ], upper = qq[3, ])
}) |>
  dplyr::bind_rows()
readr::write_csv(contr, file.path(out_dir, "archetype_contrasts.csv"))

p <- ggplot2::ggplot(grid, ggplot2::aes(cell, median, colour = archetype)) +
  ggplot2::geom_hline(yintercept = alpha, linetype = 2, colour = "grey50") +
  ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), size = 0.2,
                           position = ggplot2::position_dodge(width = 0.6)) +
  ggplot2::facet_grid(scenario ~ test, scales = "free_y") +
  ggplot2::labs(x = "Pool tree cell", y = "Rejection rate (N = 100)", colour = "Archetype") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggplot2::ggsave(file.path(out_dir, "rates.png"), p, width = 14, height = 8, dpi = 150)
