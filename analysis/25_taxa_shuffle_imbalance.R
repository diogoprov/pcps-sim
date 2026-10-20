# H1 as stated in the plan concerns the taxa shuffle alone: its rejection rate
# in the null scenarios without phylogenetic structure in the optima (S0, S2),
# by Colless level (pooled over the gamma levels), and the difference between
# high imbalance and Yule with a 95% Monte Carlo interval.
# Writes data/derived/models/taxa_shuffle_imbalance.csv.
# Usage: Rscript analysis/25_taxa_shuffle_imbalance.R

d <- data.table::fread(file.path("data", "derived", "full", "results.csv"),
                       select = c("scenario", "cell", "n_spp", "data_type", "test", "p_taxa"))
d <- d[scenario %in% c("S0", "S2") & !is.na(p_taxa)]
d[, colless := sub(".*__c_", "", cell)]
rates <- d[, .(k = sum(p_taxa <= 0.05), n = .N),
           by = .(scenario, test, n_spp, data_type, colless)]
rates[, rate := k / n]
w <- data.table::dcast(rates, scenario + test + n_spp + data_type ~ colless,
                       value.var = c("rate", "n"))
w[, diff := rate_high - rate_yule]
w[, mc_se := sqrt(rate_high * (1 - rate_high) / n_high + rate_yule * (1 - rate_yule) / n_yule)]
w[, excludes_zero := abs(diff) > 1.96 * mc_se]
readr::write_csv(w, file.path("data", "derived", "models", "taxa_shuffle_imbalance.csv"))
message("combinations: ", nrow(w), "; interval excludes zero: ", sum(w$excludes_zero),
        "; Yule ", paste(round(range(w$rate_yule), 3), collapse = "-"),
        "; high ", paste(round(range(w$rate_high), 3), collapse = "-"))
