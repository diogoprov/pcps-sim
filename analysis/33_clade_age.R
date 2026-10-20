# Recovery of the planted clade as a function of its relative age (age of its
# stem node relative to the root, recorded per replicate), to separate the
# effect of clade age from that of gamma, which also sets node ages. For shift
# 2 and incidence data: binomial GLM of hit on relative age (standardized),
# gamma level and their interaction, by method and N, and the hit rate by
# tertile of relative age within each gamma level.
# Writes data/derived/models/clade_age.csv and clade_age_glm.csv.
# Usage: Rscript analysis/33_clade_age.R

d <- data.table::fread(file.path("data", "derived", "clade", "results.csv"),
                       select = c("method", "cell", "n_spp", "data_type", "depth", "shift",
                                  "hit", "rel_age"))
d <- d[shift == 2 & data_type == "incidence"]
d[, gamma := factor(sub("__.*", "", cell), c("g_yule", "g_mid", "g_low"))]
d[, age_tertile := cut(rel_age, stats::quantile(rel_age, 0:3 / 3), include.lowest = TRUE,
                       labels = c("young", "middle", "old")), by = .(n_spp, gamma)]
tert <- d[, .(hit = mean(hit), rel_age = stats::median(rel_age), reps = .N),
          by = .(method, n_spp, gamma, age_tertile)]
readr::write_csv(tert, file.path("data", "derived", "models", "clade_age.csv"))
glms <- d[, {
  z <- as.numeric(scale(rel_age))
  f <- stats::glm(hit ~ z * gamma, family = stats::binomial())
  co <- summary(f)$coefficients
  data.table::data.table(term = rownames(co), estimate = co[, 1], se = co[, 2])
}, by = .(method, n_spp)]
readr::write_csv(glms, file.path("data", "derived", "models", "clade_age_glm.csv"))
print(data.table::dcast(tert[n_spp == 100], method + gamma ~ age_tertile, value.var = "hit"), digits = 2)
print(glms[n_spp == 100 & term %in% c("z", "z:gammag_low")], digits = 2)
