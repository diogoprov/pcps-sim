# Exploratory: per-tree PCC rejection rates (analysis/03) against tree shape.
# Usage: Rscript analysis/04_pnk_tree_shape.R (after 03 with 20 trees x 200 reps)
for (f in list.files("R", full.names = TRUE)) source(f)
colless_raw <- function(t) { n <- ape::Ntip(t); nt <- vapply(ape::prop.part(t), length, 0); ch <- split(t$edge[, 2], t$edge[, 1]); sz <- function(x) if (x <= n) 1 else nt[x - n]; sum(vapply(ch, \(k) abs(sz(k[1]) - sz(k[2])), 0)) }
res <- readr::read_csv("data/derived/validation/pnk2015_trees20_reps200.csv", show_col_types = FALSE)
set.seed(2015)
trees <- replicate(20, ape::rphylo(50, birth = 1, death = 0), simplify = FALSE)
shape <- tibble::tibble(tree_id = 1:20,
  gamma = vapply(trees, ape::gammaStat, 0),
  colless = vapply(trees, \(t) colless_raw(t), 0))
pt <- rejection_rates(res, squared, scenario, mu_tol, tree_id) |>
  dplyr::filter(!squared) |>
  dplyr::select(scenario, mu_tol, tree_id, p_taxa) |>
  tidyr::pivot_wider(names_from = c(scenario, mu_tol), values_from = p_taxa) |>
  dplyr::left_join(shape, by = "tree_id")
print(pt, n = Inf)
print(round(cor(pt[-1], method = "spearman"), 2))
