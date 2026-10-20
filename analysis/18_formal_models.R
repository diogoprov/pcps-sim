# Formal models of the simulation results (R/formal_models.R): binomial
# models aggregated by design condition, with the gamma and Colless targets of
# the grid cell as continuous covariates (standardized), for
#   type_I - layer A, scenarios under H0 (rejection rate of each test)
#   power  - layer A, scenarios under H1
#   clade  - objective 2, recovery of the planted clade (hit rate of each method)
# Each model: an exploratory GAM (mgcv, quasibinomial) and the final brms model
# (binomial, observation-level random effect for the lack of fit, explicit
# priors, 4 chains), with R-hat, ESS, divergences and a posterior predictive
# check. Outputs in data/derived/models/<model>/.
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
shapes <- cell_shapes(readRDS(file.path("data", "derived", "trees", "tree_grid_k200.rds")))

if (model %in% c("type_I", "power")) {
  scen <- if (model == "type_I") c("S0", "S1", "S2", "S2sigma", "S2h") else
    c("S3a", "S3b", "S3sigma")
  agg <- data.table::fread(file.path("data", "derived", "full", "results.csv"),
                           select = c("scenario", "cell", "n_spp", "data_type", "test",
                                      "p_max")) |>
    tibble::as_tibble() |>
    dplyr::filter(scenario %in% scen) |>
    dplyr::mutate(success = p_max <= alpha) |>
    aggregate_binomial(c("test", "scenario", "cell", "n_spp", "data_type"))
  group <- "test"
  fixed <- "test * (scenario + gamma_z + colless_z + n_f + data_type)"
  smooth <- "test * (scenario + n_f + data_type) + s(gamma_z, by = test, k = 3) +
             s(colless_z, by = test, k = 3)"
  intercept_mean <- if (model == "type_I") -2.9 else 0
} else {
  agg <- data.table::fread(file.path("data", "derived", "clade", "results.csv"),
                           select = c("cell", "n_spp", "data_type", "depth", "shift",
                                      "method", "hit")) |>
    tibble::as_tibble() |>
    dplyr::mutate(success = hit) |>
    aggregate_binomial(c("method", "cell", "n_spp", "data_type", "depth", "shift"))
  group <- "method"
  fixed <- "method * (gamma_z * depth + colless_z + shift_f + n_f + data_type)"
  smooth <- "method * (depth + shift_f + n_f + data_type) +
             s(gamma_z, by = interaction(method, depth), k = 3) +
             s(colless_z, by = method, k = 3)"
  intercept_mean <- 0
}
agg <- agg |>
  dplyr::left_join(shapes, by = c("n_spp", "cell")) |>
  dplyr::mutate(n_f = factor(n_spp), dplyr::across(dplyr::all_of(group), factor))
if ("shift" %in% names(agg)) agg$shift_f <- factor(agg$shift)
stopifnot(!anyNA(agg$gamma_z), length(unique(agg$n)) == 1)
readr::write_csv(agg, file.path(out_dir, "aggregated.csv"))
message(model, ": ", nrow(agg), " conditions")

# GAM (exploratory) ---------------------------------------------------------------
gam_fit <- mgcv::gam(stats::as.formula(paste("cbind(k, n - k) ~", smooth)),
                     family = stats::quasibinomial(), data = agg, method = "REML")
s <- summary(gam_fit)
readr::write_csv(tibble::tibble(term = rownames(s$s.table), edf = s$s.table[, "edf"],
                                F = s$s.table[, "F"], p_value = s$s.table[, "p-value"],
                                dev_expl = s$dev.expl, scale = gam_fit$scale),
                 file.path(out_dir, "gam_smooths.csv"))
message("GAM: deviance explained ", round(s$dev.expl, 3), ", scale ", round(gam_fit$scale, 2))

# brms (final) ---------------------------------------------------------------------
t0 <- Sys.time()
fit <- brms::brm(
  stats::as.formula(paste("k | trials(n) ~", fixed, "+ (1 | obs)")),
  data = agg, family = stats::binomial(), prior = model_priors(intercept_mean),
  chains = chains, cores = chains, iter = iter, warmup = iter %/% 2, seed = 20261019,
  backend = "cmdstanr", control = list(adapt_delta = 0.95),
  file = file.path(out_dir, "brms_fit"), refresh = 0
)
message("brms: ", format(Sys.time() - t0))
diag <- dplyr::bind_cols(fit_diagnostics(fit), ppc_coverage(fit))
readr::write_csv(diag, file.path(out_dir, "diagnostics.csv"))
print(diag)
fe <- brms::fixef(fit)
readr::write_csv(tibble::as_tibble(fe, rownames = "term"), file.path(out_dir, "fixef.csv"))

# posterior predictive check plot and gamma effects by test or method
grDevices::png(file.path(out_dir, "ppc.png"), width = 1400, height = 900, res = 150)
print(brms::pp_check(fit, type = "scatter_avg", ndraws = min(200, brms::ndraws(fit))))
grDevices::dev.off()
ce <- brms::conditional_effects(fit, effects = paste0("gamma_z:", group),
                                conditions = data.frame(n = 1))
ggplot2::ggsave(file.path(out_dir, "gamma_effect.png"), plot(ce, plot = FALSE)[[1]],
                width = 9, height = 6, dpi = 150)
