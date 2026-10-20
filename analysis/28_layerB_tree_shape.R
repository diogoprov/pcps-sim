# Shape of the pool tree and of the analysed tree in Layer B (one row per
# replicate: test PERMANOVA only, since the trees do not depend on the test),
# by archetype, pool cell and N: medians of gamma and of the normalized Colless
# index, and the number of species left after the dynamics.
# Writes data/derived/models/layerB_tree_shape.csv.
# Usage: Rscript analysis/28_layerB_tree_shape.R

d <- data.table::fread(file.path("data", "derived", "layerB", "results.csv"),
                       select = c("scenario", "archetype", "cell", "n_spp", "test", "n_realized",
                                  "gamma_pool", "colless_pool", "gamma", "colless"))
d <- d[test == "adonis"]
s <- d[, .(gamma_pool = stats::median(gamma_pool), gamma_analysed = stats::median(gamma),
           colless_pool = stats::median(colless_pool), colless_analysed = stats::median(colless),
           species_left = stats::median(n_realized)),
       by = .(archetype, cell, n_spp)]
readr::write_csv(s, file.path("data", "derived", "models", "layerB_tree_shape.csv"))
print(s[n_spp == 100][order(cell, archetype)], digits = 3)
