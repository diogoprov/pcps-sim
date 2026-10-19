# Layer B: binomial GAM of the rejection rate on the realized shape (gamma,
# Colless) of the analysed tree, within each test x scenario x archetype x N
# (R/rejection_gam.R). Asks whether the loss of power at low gamma of the pool
# also holds along the gamma of the tree that is analysed. The realized and the
# pool gamma are strongly correlated, so this does not separate the two.
# Usage: Rscript analysis/14_gam_layerB.R <workers>

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[[1]]) else 8L
alpha <- 0.05

for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
in_file <- file.path("data", "derived", "layerB", "results.csv")
out_dir <- file.path("data", "derived", "layerB_gam")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

res <- data.table::fread(in_file, select = c("scenario", "archetype", "n_spp", "test",
                                             "p_max", "gamma", "colless")) |>
  tibble::as_tibble()
stopifnot(!anyNA(res))
strata <- split(res, list(res$test, res$scenario, res$archetype, res$n_spp), drop = TRUE)

t0 <- Sys.time()
fits <- parallel::mclapply(strata, \(d) {
  fit <- fit_rejection_gam(d, alpha)
  key <- d[1, c("test", "scenario", "archetype", "n_spp")]
  list(summary = dplyr::bind_cols(key, summarise_rejection_gam(fit), n = nrow(d)),
       curve = dplyr::bind_cols(key, rejection_curve_gamma(fit)))
}, mc.cores = workers)
# a worker that errors returns a try-error, one that is killed returns NULL
stopifnot(length(fits) == length(strata),
          vapply(fits, \(x) is.list(x) && !inherits(x, "try-error"), logical(1)))
message(sprintf("%d models in %.1f min", length(fits),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

smooths <- dplyr::bind_rows(lapply(fits, `[[`, "summary"))
curves <- dplyr::bind_rows(lapply(fits, `[[`, "curve"))
readr::write_csv(smooths, file.path(out_dir, "gam_smooths.csv"))
readr::write_csv(curves, file.path(out_dir, "gam_gamma_curves.csv"))

# N = 100 is the only N with all four archetypes
p <- curves |>
  dplyr::filter(n_spp == 100) |>
  ggplot2::ggplot(ggplot2::aes(gamma, rate, colour = test, fill = test)) +
  ggplot2::geom_hline(yintercept = alpha, linetype = 2, colour = "grey50") +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  ggplot2::geom_line() +
  ggplot2::facet_grid(archetype ~ scenario, scales = "free") +
  ggplot2::labs(x = expression("Realized " * gamma * " of the analysed tree"),
                y = "Rejection rate (N = 100)", colour = "Test", fill = "Test") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(out_dir, "gam_gamma_curves_n100.png"), p, width = 16, height = 9, dpi = 150)
