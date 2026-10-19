# Layer B: binomial GAM of the rejection rate on the realized shape (gamma,
# Colless) of the analysed tree, with archetype and N as fixed effects, for
# each test x scenario (R/rejection_gam.R). Asks whether the loss of power at
# low gamma of the pool also holds along the gamma of the tree that is analysed.
# Usage: Rscript analysis/14_gam_layerB.R <workers> [max rows per model]

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[[1]]) else 8L
max_rows <- if (length(args) >= 2) as.integer(args[[2]]) else Inf

for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
in_file <- file.path("data", "derived", "layerB", "results.csv")
out_dir <- file.path("data", "derived", "layerB_gam")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

res <- data.table::fread(in_file, select = c("scenario", "archetype", "n_spp", "test",
                                             "p_max", "gamma", "colless")) |>
  tibble::as_tibble()
stopifnot(!anyNA(res))
groups <- split(res, list(res$test, res$scenario), drop = TRUE)

set.seed(20261019L)
groups <- lapply(groups, \(d) if (nrow(d) > max_rows) d[sample.int(nrow(d), max_rows), ] else d)

t0 <- Sys.time()
fits <- parallel::mclapply(groups, \(d) {
  fit <- fit_rejection_gam(d)
  key <- d[1, c("test", "scenario")]
  list(summary = dplyr::bind_cols(key, summarise_rejection_gam(fit)),
       curve = dplyr::bind_cols(key, marginal_rejection_gamma(fit, d)))
}, mc.cores = workers, mc.set.seed = TRUE)
stopifnot(!vapply(fits, inherits, logical(1), "try-error"))
message(sprintf("%d models in %.1f min", length(fits),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

smooths <- dplyr::bind_rows(lapply(fits, `[[`, "summary"))
curves <- dplyr::bind_rows(lapply(fits, `[[`, "curve"))
readr::write_csv(smooths, file.path(out_dir, "gam_smooths.csv"))
readr::write_csv(curves, file.path(out_dir, "gam_gamma_curves.csv"))

p <- ggplot2::ggplot(curves, ggplot2::aes(gamma, rate, colour = test)) +
  ggplot2::geom_hline(yintercept = 0.05, linetype = 2, colour = "grey50") +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~scenario, scales = "free_y") +
  ggplot2::labs(x = expression("Realized " * gamma * " of the analysed tree"),
                y = "Marginal rejection rate", colour = "Test") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(out_dir, "gam_gamma_curves.png"), p, width = 10, height = 7, dpi = 150)
