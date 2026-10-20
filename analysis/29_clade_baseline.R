# Chance level of the clade-recovery methods: hit rate, median rank of the
# planted node and mean Jaccard similarity with shift 0 (the planted clade has
# no effect on the optima), against shift 0.5 and 2, by method, N, depth and
# data type, averaged over the cells of the grid.
# Writes data/derived/models/clade_baseline.csv.
# Usage: Rscript analysis/29_clade_baseline.R

d <- data.table::fread(file.path("data", "derived", "clade", "results.csv"),
                       select = c("method", "cell", "n_spp", "data_type", "depth", "shift",
                                  "hit", "rank", "jaccard"))
s <- d[shift %in% c(0, 0.5, 2),
       .(hit = mean(hit), hit_max_cell = max(tapply(hit, cell, mean)),
         median_rank = stats::median(rank, na.rm = TRUE), jaccard = mean(jaccard, na.rm = TRUE)),
       by = .(method, n_spp, data_type, depth, shift)]
readr::write_csv(s, file.path("data", "derived", "models", "clade_baseline.csv"))
print(data.table::dcast(s[data_type == "incidence"], method + n_spp + depth ~ shift,
                        value.var = "hit"), digits = 3)
print(s[shift == 0, .(hit_range = paste(round(range(hit), 3), collapse = "-"),
                      max_cell = round(max(hit_max_cell), 3),
                      rank = paste(range(median_rank), collapse = "-")), by = method])
