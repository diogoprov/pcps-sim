# All treestats statistics for the grid trees and the empirical anuran clades,
# to relate rejection rates to tree shape beyond gamma and Colless.
# Usage: Rscript analysis/08_tree_stats.R <workers>  (after analysis/05 and 06)

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[[1]]) else 16L

for (f in list.files("R", full.names = TRUE)) source(f)
future::plan(future::multisession, workers = workers)
opts <- furrr::furrr_options(seed = TRUE, globals = TRUE)
out_dir <- file.path("data", "derived", "trees")

grid <- readRDS(file.path(out_dir, "tree_grid_k200.rds"))
jobs <- grid |>
  dplyr::select(n, cell, trees) |>
  dplyr::mutate(tree_id = purrr::map(trees, seq_along)) |>
  tidyr::unnest(c(trees, tree_id))

t0 <- Sys.time()
stats_grid <- furrr::future_map(jobs$trees, tree_stats_all, .options = opts) |>
  dplyr::bind_rows()
message("Grid trees: ", format(Sys.time() - t0))
dplyr::bind_cols(dplyr::select(jobs, n, cell, tree_id), stats_grid) |>
  readr::write_csv(file.path(out_dir, "tree_stats_grid_k200.csv"))

# Empirical clades (same extraction and scaling as analysis/05)
tr <- ape::read.tree("data/raw/portik2023_timetree.tre")
genus <- sub("_.*", "", tr$tip.label)
anura <- ape::extract.clade(
  tr, ape::getMRCA(tr, tr$tip.label[genus %in% c("Ascaphus", "Rana")])
)
sizes <- vapply(ape::prop.part(anura), length, 0)
nodes <- ape::Ntip(anura) + which(sizes >= 40 & sizes <= 250)
clades <- purrr::map(nodes, \(v) {
  t <- ape::extract.clade(anura, v)
  t$edge.length <- t$edge.length / max(ape::node.depth.edgelength(t))
  phytools::force.ultrametric(t, method = "extend", message = FALSE)
})
stats_emp <- furrr::future_map(clades, tree_stats_all, .options = opts) |>
  dplyr::bind_rows()
dplyr::bind_cols(tibble::tibble(node = nodes, n = vapply(clades, ape::Ntip, 0)), stats_emp) |>
  readr::write_csv(file.path(out_dir, "tree_stats_empirical_portik2023.csv"))
message("Done: ", format(Sys.time() - t0))
