# GLM on PCPS axes with symmetric Procrustes of the null axes onto the observed
# ones (as in PCPS::pcps.sig), and RDA of PCPS axes 1-2 on the gradient (which
# needs no Procrustes), against the plain GLM used in the main text: raw
# rejection rates of the max test, with Monte Carlo standard errors, by layer,
# scenario, gamma level, N and data type (Layer B: archetype and pool cell).
# Layer A from the run with extra axes (data/derived/axes/results.csv),
# Layer B from data/derived/layerB/results.csv.
# Writes data/derived/models/procrustes_rates.csv.
# Usage: Rscript analysis/27_procrustes_variant.R

tests <- c("glm_pcps1", "glm_pcps1_proc", "glm_pcps2", "glm_pcps2_proc", "rda_pcps")
rate_table <- function(file, by) {
  d <- data.table::fread(file, select = c(by, "test", "p_max"))
  d <- d[test %in% tests]
  d[, .(k = sum(p_max <= 0.05), n = .N), by = c(by, "test")][
    , `:=`(rate = k / n, mcse = sqrt(k / n * (1 - k / n) / n))]
}
a <- rate_table(file.path("data", "derived", "axes", "results.csv"),
                c("scenario", "cell", "n_spp", "data_type"))
a[, `:=`(layer = "A", archetype = NA_character_)]
b <- rate_table(file.path("data", "derived", "layerB", "results.csv"),
                c("scenario", "archetype", "cell", "n_spp", "data_type"))
b[, layer := "B"]
out <- data.table::rbindlist(list(a, b), use.names = TRUE)
out[, gamma := sub("__.*", "", cell)]
readr::write_csv(out, file.path("data", "derived", "models", "procrustes_rates.csv"))

# summary: type I in S2sigma and power in S3b, by gamma level (mean over the
# Colless levels), at N = 100 and 200 with incidence data
s <- out[data_type == "incidence" & n_spp %in% c(100, 200) & scenario %in% c("S2", "S2sigma", "S3b"),
         .(rate = round(mean(rate), 3)), by = .(layer, scenario, test, n_spp, gamma)]
print(data.table::dcast(s, layer + scenario + n_spp + gamma ~ test, value.var = "rate"))
