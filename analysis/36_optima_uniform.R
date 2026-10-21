# Results of the run with optima of fixed marginal distribution (profile
# optima_uniform): raw power of every test by cell, N and data type, with
# Monte Carlo SE, and the loss with low gamma (low minus Yule gamma, at Yule
# Colless) next to the same loss in the main run (and, for RDA PCPS and the
# Procrustes GLM, in the run with the extra axes).
# Writes data/derived/models/optima_uniform_rates.csv and optima_uniform_loss.csv.
# Usage: Rscript analysis/36_optima_uniform.R

rates <- function(file) {
  d <- data.table::fread(file, select = c("scenario", "cell", "n_spp", "data_type", "test", "p_max"))
  d <- d[scenario %in% c("S3a", "S3b", "S3sigma")]
  d[, .(k = sum(p_max <= 0.05), n = .N), by = .(scenario, cell, n_spp, data_type, test)][
    , `:=`(rate = k / n, mcse = sqrt(k / n * (1 - k / n) / n))]
}
u <- rates(file.path("data", "derived", "optima_uniform", "results.csv"))
readr::write_csv(u, file.path("data", "derived", "models", "optima_uniform_rates.csv"))
main <- data.table::rbindlist(list(
  rates(file.path("data", "derived", "full", "results.csv")),
  rates(file.path("data", "derived", "axes", "results.csv"))[test %in% c("rda_pcps", "glm_pcps1_proc", "glm_pcps2_proc")]
))
loss <- function(d, run) {
  w <- data.table::dcast(d[cell %in% c("g_low__c_yule", "g_yule__c_yule")],
                         scenario + n_spp + data_type + test ~ cell, value.var = "rate")
  w[, .(scenario, n_spp, data_type, test, yule = g_yule__c_yule, low = g_low__c_yule,
        loss = g_low__c_yule - g_yule__c_yule, run = run)]
}
out <- data.table::rbindlist(list(loss(u, "uniform optima"), loss(main, "main")))
readr::write_csv(out, file.path("data", "derived", "models", "optima_uniform_loss.csv"))
print(data.table::dcast(out[scenario == "S3b" & data_type == "incidence"],
                        test + n_spp ~ run, value.var = "loss"), digits = 2)
