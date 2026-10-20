# Paired Monte Carlo contrasts between scenarios: replicate r of every scenario
# uses tree r of its cell, so the difference in rejection between S2sigma (or
# S2h) and S2 can be computed replicate by replicate. For each test, cell, N
# and data type: the raw difference, its paired Monte Carlo standard error
# (sd of the replicate differences / sqrt(K)) and the unpaired one, for Layer A.
# Writes data/derived/models/paired_contrasts.csv.
# Usage: Rscript analysis/31_paired_contrasts.R

d <- data.table::fread(file.path("data", "derived", "full", "results.csv"),
                       select = c("scenario", "cell", "n_spp", "data_type", "test", "rep", "p_max"))
d <- d[scenario %in% c("S2", "S2sigma", "S2h")]
d[, reject := as.integer(p_max <= 0.05)]
w <- data.table::dcast(d, cell + n_spp + data_type + test + rep ~ scenario, value.var = "reject")
out <- data.table::rbindlist(lapply(c("S2sigma", "S2h"), \(sc) {
  w[!is.na(get(sc)) & !is.na(S2),
    .(contrast = paste(sc, "- S2"), K = .N, diff = mean(get(sc) - S2),
      se_paired = stats::sd(get(sc) - S2) / sqrt(.N),
      se_unpaired = sqrt(mean(get(sc)) * (1 - mean(get(sc))) / .N + mean(S2) * (1 - mean(S2)) / .N)),
    by = .(cell, n_spp, data_type, test)]
}))
out[, `:=`(excl_paired = abs(diff) > 1.96 * se_paired, excl_unpaired = abs(diff) > 1.96 * se_unpaired)]
readr::write_csv(out, file.path("data", "derived", "models", "paired_contrasts.csv"))
print(out[, .(n = .N, se_ratio = round(stats::median(se_paired / se_unpaired), 3),
              excl_paired = sum(excl_paired), excl_unpaired = sum(excl_unpaired)),
          by = .(contrast, test)])
