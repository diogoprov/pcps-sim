# Empirical reference for the tree-shape grid: gamma and Yule-normalized Colless
# of anuran subtrees from the Portik et al. 2023 timetree.
#   (a) clades with 40-250 tips (nested clades included, so not independent);
#   (b) random tip subsets of 50, 100 and 200 anuran species (incomplete sampling
#       pushes gamma down, so these bracket the low end).
# Usage: Rscript analysis/05_tree_shape_empirical.R

for (f in list.files("R", full.names = TRUE)) source(f)

# Subtrees inherit rounding error in the branch lengths (Myr, 6 decimals); scale
# to root depth 1 and extend tips to remove it before the ultrametric check.
shape_emp <- function(t) {
  t$edge.length <- t$edge.length / max(ape::node.depth.edgelength(t))
  as.list(tree_shape(phytools::force.ultrametric(t, method = "extend", message = FALSE)))
}
out_dir <- file.path("data", "derived", "trees")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tr <- ape::read.tree("data/raw/portik2023_timetree.tre")
genus <- sub("_.*", "", tr$tip.label)
anura <- ape::extract.clade(
  tr, ape::getMRCA(tr, tr$tip.label[genus %in% c("Ascaphus", "Rana")])
)
message("Anura: ", ape::Ntip(anura), " tips")

n_anura <- ape::Ntip(anura)
sizes <- vapply(ape::prop.part(anura), length, 0)
nodes <- n_anura + which(sizes >= 40 & sizes <= 250)
clades <- purrr::map(nodes, \(v) ape::extract.clade(anura, v))
emp_clades <- tibble::tibble(
  source = "clade", node = nodes, n = vapply(clades, ape::Ntip, 0)
) |>
  dplyr::bind_cols(dplyr::bind_rows(purrr::map(clades, shape_emp)))

set.seed(2023)
emp_random <- tidyr::crossing(n = c(50, 100, 200), rep = 1:200) |>
  dplyr::mutate(shape = purrr::map(n, \(k) {
    keep <- sample(anura$tip.label, k)
    shape_emp(ape::keep.tip(anura, keep))
  })) |>
  tidyr::unnest_wider(shape) |>
  dplyr::mutate(source = "random")

emp <- dplyr::bind_rows(emp_clades, emp_random)
readr::write_csv(emp, file.path(out_dir, "empirical_shape_portik2023.csv"))

emp |>
  dplyr::mutate(size = cut(n, c(0, 75, 150, Inf), labels = c("~50", "~100", "~200"))) |>
  dplyr::summarise(
    k = dplyr::n(),
    dplyr::across(c(gamma, colless), list(
      q10 = \(x) quantile(x, 0.1), q50 = median, q90 = \(x) quantile(x, 0.9))),
    .by = c(source, size)
  ) |>
  dplyr::arrange(source, size) |>
  print(width = Inf)
