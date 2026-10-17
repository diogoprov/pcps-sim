# Validation of layer A against the published tables:
#   Duarte et al. 2016 (MEE) Table 1 -- PFW, Mantel and ADONIS, mu_tol = 10
#   Peres-Neto & Kembel 2015 (Plant Ecol) Table 1 -- PCC, one gradient
# Usage: Rscript analysis/01_validate_layerA.R <n_reps> <workers>

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[[1]]) else 200L
workers <- if (length(args) >= 2) as.integer(args[[2]]) else 16L

for (f in list.files("R", full.names = TRUE)) source(f)
future::plan(future::multisession, workers = workers)
opts <- furrr::furrr_options(seed = 20261016L, globals = TRUE,
                             packages = c("dplyr", "tibble"))
out_dir <- file.path("data", "derived", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Duarte et al. 2016 -------------------------------------------------------
grid_d <- tidyr::crossing(scenarios_4(), n_spp = c(50, 100, 200),
                          rep = seq_len(n_reps))
t0 <- Sys.time()
res_d <- furrr::future_pmap(
  grid_d,
  \(scenario, phylo, env, n_spp, rep) rep_duarte2016(n_spp, phylo, env) |>
    dplyr::mutate(scenario = scenario, n_spp = n_spp, rep = rep),
  .options = opts
) |> dplyr::bind_rows()
message("Duarte 2016: ", format(Sys.time() - t0))
readr::write_csv(res_d, file.path(out_dir, sprintf("duarte2016_reps%d.csv", n_reps)))

# Peres-Neto & Kembel 2015 ---------------------------------------------------
set.seed(2015)
tree_pnk <- ape::rphylo(50, birth = 1, death = 0)
grid_p <- tidyr::crossing(scenarios_4(), mu_tol = c(5, 10),
                          squared = c(FALSE, TRUE), rep = seq_len(n_reps))
t0 <- Sys.time()
res_p <- furrr::future_pmap(
  grid_p,
  \(scenario, phylo, env, mu_tol, squared, rep)
    rep_pnk2015(tree_pnk, phylo, env, mu_tol = mu_tol, squared = squared) |>
    dplyr::mutate(scenario = scenario, mu_tol = mu_tol, squared = squared, rep = rep),
  .options = opts
) |> dplyr::bind_rows()
message("PNK 2015: ", format(Sys.time() - t0))
readr::write_csv(res_p, file.path(out_dir, sprintf("pnk2015_reps%d.csv", n_reps)))

print(rejection_rates(res_d, test, scenario, n_spp), n = Inf)
print(rejection_rates(res_p, squared, scenario, mu_tol), n = Inf)
