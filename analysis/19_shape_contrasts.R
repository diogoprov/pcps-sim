# Posterior contrasts of the tree-shape effects from the brms models of
# analysis/18_formal_models.R: rate at low gamma minus Yule, and at high
# Colless minus Yule, for every test x scenario (or method x depth), at N = 100,
# incidence (clade: shift 2). Median and 95% credible interval.
# Usage: Rscript analysis/19_shape_contrasts.R

for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
root <- file.path("data", "derived", "models")

contrasts <- purrr::map(c("type_I", "power", "clade"), \(m) {
  fit <- readRDS(file.path(root, m, "brms_fit.rds"))
  grp <- if (m == "clade") c("method", "depth") else c("test", "scenario")
  base <- tidyr::crossing(!!!lapply(stats::setNames(grp, grp), \(v) levels(fit$data[[v]])),
                          n_f = "100", data_type = "incidence", n = 1)
  if (m == "clade") base$shift_f <- "2"
  at <- \(g, cl) brms::posterior_epred(fit, newdata = dplyr::mutate(base, g = g, cl = cl))
  yule <- at("g_yule", "c_yule")
  summ <- \(d, what) {
    q <- apply(d, 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
    dplyr::mutate(base[grp], model = m, contrast = what, median = q[2, ], lower = q[1, ],
                  upper = q[3, ])
  }
  dplyr::bind_rows(summ(at("g_low", "c_yule") - yule, "gamma low - yule"),
                   summ(at("g_yule", "c_high") - yule, "colless high - yule"))
}) |>
  dplyr::bind_rows()
readr::write_csv(contrasts, file.path(root, "shape_contrasts.csv"))
