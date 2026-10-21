# Results of the spatial version of layer B (profile layerB_spatial): raw
# rejection rates with Monte Carlo SE by test, scenario, archetype, pool cell
# and data type (max test, with torus translations as the site null), and the
# same rates of the non-spatial layer B (incidence, N = 100) for comparison;
# share of replicates with empty sites dropped (where the torus null is not
# exact).
# Writes data/derived/models/layerB_spatial_rates.csv.
# Usage: Rscript analysis/37_layerB_spatial.R

cols <- c("scenario", "archetype", "cell", "n_spp", "data_type", "test", "p_max", "n_dropped_sites")
rates <- function(file, run) {
  d <- data.table::fread(file, select = cols)[n_spp == 100]
  d[, .(k = sum(p_max <= 0.05), n = .N, dropped = mean(n_dropped_sites > 0)),
    by = .(scenario, archetype, cell, data_type, test)][
    , `:=`(rate = k / n, mcse = sqrt(k / n * (1 - k / n) / n), run = run)]
}
out <- data.table::rbindlist(list(
  rates(file.path("data", "derived", "layerB_spatial", "results.csv"), "spatial"),
  rates(file.path("data", "derived", "layerB", "results.csv"), "non-spatial")
))
readr::write_csv(out, file.path("data", "derived", "models", "layerB_spatial_rates.csv"))
message("replicates with dropped sites (spatial): ",
        round(stats::weighted.mean(out[run == "spatial"]$dropped, out[run == "spatial"]$n), 4))
s <- out[cell == "g_yule__c_yule" & scenario %in% c("S0", "S2", "S2sigma", "S3b")]
print(data.table::dcast(s, archetype + scenario + test ~ run + data_type, value.var = "rate"),
      digits = 2, nrows = 200)
