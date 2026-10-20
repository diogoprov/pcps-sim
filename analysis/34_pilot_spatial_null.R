# Pilot of the spatial version of layer B under the null hypothesis (S0: the
# analysed gradient is independent of the one that generated the data, but
# both are spatially autocorrelated): type I error of the site test with torus
# translations (as in the layerB_spatial profile) against the plain site
# shuffle, under species sorting and mass effects without competition, in Yule
# pool trees at N = 100, incidence data.
# Writes data/derived/pilot_spatial/results.csv.
# Usage: Rscript analysis/34_pilot_spatial_null.R <reps> <workers>

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 200L
workers <- if (length(args) >= 2) as.integer(args[[2]]) else 16L
for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
future::plan(future::multisession, workers = workers)
out_dir <- file.path("data", "derived", "pilot_spatial")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pools <- feasible_pools(layerB_pool_file())
jobs <- tidyr::crossing(archetype = c("sorting__none", "mass_effects__none"), rep = seq_len(reps)) |>
  dplyr::left_join(archetypes_layerB(), by = "archetype")
res <- furrr::future_pmap(jobs, \(archetype, rep, regime, dispersal, competition) {
  sp <- spatial_landscape()
  n_pool <- pools$n_pool[pools$archetype == archetype & pools$scenario == "S0" & pools$target == 100]
  tr <- sim_tree_shape(n_pool)
  mc <- sim_metacommunity_B(tr, 0, FALSE, FALSE, FALSE, dispersal = dispersal,
                            competition = competition, incidence = TRUE, space = sp)
  torus <- analyse_metacommunity(mc, c("pfw", "glm", "pcc"), 999, 999, extra_axes = TRUE,
                                 site_env = torus_site_env(mc, sp))
  plain <- analyse_metacommunity(mc, c("pfw", "glm", "pcc"), 999, 999, extra_axes = TRUE)
  dplyr::bind_rows(dplyr::mutate(torus, null = "torus"), dplyr::mutate(plain, null = "site shuffle")) |>
    dplyr::mutate(archetype = archetype, rep = rep)
}, .options = furrr::furrr_options(seed = 20261020L)) |>
  dplyr::bind_rows()
readr::write_csv(res, file.path(out_dir, "results.csv"))
s <- dplyr::summarise(dplyr::group_by(res, archetype, null, test),
                      site = mean(p_site <= 0.05), max = mean(p_max <= 0.05), .groups = "drop")
print(tidyr::pivot_wider(s, names_from = null, values_from = c(site, max)), n = Inf)
