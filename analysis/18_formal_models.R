# Formal models of the simulation results (R/formal_models.R): binomial
# models aggregated by design condition (1000 replicates each), for
#   type_I - layer A, scenarios under H0 (rejection rate of each test)
#   power  - layer A, scenarios under H1
#   clade  - objective 2, recovery of the planted clade (hit rate of each method)
# Tree shape enters as the gamma and Colless levels of the grid cell (3 x 3
# factors, with their interaction: the gamma x imbalance surface). The plan asks
# for continuous shape covariates, but every cell has a single target per axis
# and layer A did not keep the shape of each tree, so a "continuous" term would
# only interpolate three points (mostly Yule vs the rest); the factors keep the
# surface without assuming a form. Each model: an exploratory quasibinomial GLM
# (for the dispersion) and the final brms model (beta-binomial for the
# overdispersion, explicit priors, 4 chains), with R-hat, ESS, divergences and a
# posterior predictive check. An observation-level random effect was tried
# first: with 1000 trials per condition its ~1000 latent values made almost
# every transition hit the maximum tree depth. Outputs in <out_root>/<model>/.
# Usage: Rscript analysis/18_formal_models.R <type_I|power|clade> [chains] [iter] [out_root]

args <- commandArgs(trailingOnly = TRUE)
model <- match.arg(args[[1]], c("type_I", "power", "clade"))
chains <- if (length(args) >= 2) as.integer(args[[2]]) else 4L
iter <- if (length(args) >= 3) as.integer(args[[3]]) else 2000L
out_root <- if (length(args) >= 4) args[[4]] else file.path("data", "derived", "models")

for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
out_dir <- file.path(out_root, model)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
alpha <- 0.05

# H2 needs the shape effects to differ among scenarios (scenario x shape); H4
# needs them to differ among tests or methods; H3 needs gamma x depth. The
# shape effects may also differ with N and data type (shape x N, shape x data
# type), so the contrasts at one N do not stand for the others by construction
layer_a <- paste("test * scenario * (g + cl) + test * (g:cl + n_f + data_type) +",
                 "test:(g + cl):(n_f + data_type)")
spec <- list(
  type_I = list(file = c("full", "results.csv"), group = "test",
                scenarios = c("S0", "S1", "S2", "S2sigma", "S2h"), fixed = layer_a,
                intercept_mean = -2.9),
  power = list(file = c("full", "results.csv"), group = "test",
               scenarios = c("S3a", "S3b", "S3sigma"), fixed = layer_a,
               intercept_mean = 0),
  clade = list(file = c("clade", "results.csv"), group = "method",
               fixed = paste("method * (g * depth + cl + g:cl + shift_f + n_f + data_type) +",
                             "method:(g + cl):(n_f + data_type + shift_f)"),
               intercept_mean = 0)
)[[model]]

cols <- if (model == "clade") c("method", "cell", "n_spp", "data_type", "depth", "shift", "hit") else
  c("test", "scenario", "cell", "n_spp", "data_type", "p_max")
raw <- data.table::fread(do.call(file.path, as.list(c("data", "derived", spec$file))),
                         select = cols) |>
  tibble::as_tibble()
if (model == "clade") {
  raw <- dplyr::mutate(raw, success = hit)
  by <- c("method", "cell", "n_spp", "data_type", "depth", "shift")
} else {
  raw <- raw |>
    dplyr::filter(scenario %in% spec$scenarios) |>
    dplyr::mutate(success = p_max <= alpha)
  by <- c("test", "scenario", "cell", "n_spp", "data_type")
}
stopifnot(!anyNA(raw$success))
agg <- aggregate_binomial(raw, by) |>
  dplyr::mutate(g = factor(sub("__.*", "", cell), c("g_yule", "g_mid", "g_low")),
                cl = factor(sub(".*__", "", cell), c("c_yule", "c_mid", "c_high")),
                n_f = factor(n_spp),
                dplyr::across(dplyr::any_of(c(spec$group, "scenario", "data_type", "depth")),
                              factor))
if (model == "clade") agg$shift_f <- factor(agg$shift)
stopifnot(!anyNA(agg$g), !anyNA(agg$cl), length(unique(agg$n)) == 1)
readr::write_csv(agg, file.path(out_dir, "aggregated.csv"))
message(model, ": ", nrow(agg), " conditions")
rm(raw)

# Fit and checks (R/formal_models.R). Grouped check: sums over the conditions of
# each test (or method) x N x cell, where the fixed structure must reproduce the
# shape x N pattern
fit <- fit_formal_model(spec$fixed, agg,
                        groups = interaction(agg[[spec$group]], agg$n_f, agg$cell, drop = TRUE),
                        out_dir = out_dir, intercept_mean = spec$intercept_mean,
                        chains = chains, iter = iter)

# Key figure: rate over the gamma x Colless grid by test or method, at N = 100 and
# incidence (clade: shift 2)
grid <- tidyr::crossing(g = levels(agg$g), cl = levels(agg$cl),
                        grp = levels(agg[[spec$group]]), n_f = "100",
                        data_type = "incidence", n = 1)
names(grid)[names(grid) == "grp"] <- spec$group
grid <- if (model == "clade") {
  tidyr::crossing(grid, depth = levels(agg$depth), shift_f = "2")
} else {
  tidyr::crossing(grid, scenario = levels(agg$scenario))
}
ep <- brms::posterior_epred(fit, newdata = grid)
grid$rate <- apply(ep, 2, stats::median)
readr::write_csv(grid, file.path(out_dir, "shape_surface.csv"))
facet_row <- if (model == "clade") "depth" else "scenario"
p <- ggplot2::ggplot(grid, ggplot2::aes(cl, g, fill = rate)) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rate)), size = 2.5) +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::facet_grid(stats::reformulate(spec$group, facet_row)) +
  ggplot2::labs(x = "Colless level", y = "Gamma level", fill = "Rate") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(out_dir, "shape_surface.png"), p, width = 12, height = 8, dpi = 150)
