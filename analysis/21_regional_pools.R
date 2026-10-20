# Regional species pools from the empirical tadpole metacommunities
# (project "NRI vs PCPS", 2012; local copy in data/raw/nri_pcps/, not in git:
# compiled from the literature and from unpublished data shared with the
# author), as an empirical reference for the tree-shape grid.
#
# Step 1, names:
#   - taxon labels of each metacommunity from its NEXUS file (2012 names); Both
#     from Both_tree.txt (its NEXUS has codes only); Mainara from a table of its
#     codes (below);
#   - other abbreviated labels ("Bokcir", "A.albosignatus", "Halbopunctatus")
#     resolved by genus prefix and species-epithet prefix (the epithet stem,
#     without the last two letters, so that gender changes still match), first
#     against the full 2012 names of the other datasets, then against the valid
#     and then all ASW names (AmphiNom) of species recorded in Brazil,
#     Argentina, Paraguay, Uruguay, Bolivia or Peru (AmphibiaWeb ranges);
#   - names harmonized with the Amphibian Species of the World (ASW synonymy,
#     then asw_sync) and matched to the tips of the Portik et al. 2023 tree,
#     harmonized the same way;
#   - the rules set on 2026-10-19: undetermined taxa (sp., aff., cf., gr.,
#     "Xsp" codes) and codes that stay ambiguous are excluded; a name the ASW
#     does not know is replaced by the closest ASW name within 2 edits, when
#     unique (spelling errors; the replacement may also apply ASW synonymy).
# Step 2: three nested regional pools pruned from the Portik tree and from its
# 100 time-calibrated bootstrap trees, with gamma and Colless (normalized
# against Yule, R/trees.R).
# Writes, in data/derived/pools/: names.csv (every label, its resolution and
# status), datasets.csv, pool_species.csv and pool_shapes.csv.
# Usage: Rscript analysis/21_regional_pools.R

for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
raw_dir <- file.path("data", "raw", "nri_pcps")
out_dir <- file.path("data", "derived", "pools")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# one file per metacommunity (Wild19961.nex duplicates Wild1996.nex), with its
# locality and region; `locality_source` says whether it was confirmed
# or it is inferred from the file names, the coordinates file and the species
datasets <- tibble::tribble(
  ~dataset,      ~file,                       ~locality,                  ~region,              ~locality_source,
  "bocaina",     "bocaina_tadpole.trees.nex", "Serra da Bocaina, SP/RJ",  "Atlantic Forest SE", "inferred",
  "both",        "Both_tree.txt",             "Cacapava do Sul, RS",      "South (Pampa)",      "inferred",
  "costa",       "Costa_phylo.nex",           "Maranhao?",                "Northeast?",         "inferred",
  "crivelari",   "Crivelari_phylo.nex",       "Parana",                   "Atlantic Forest S",  "author",
  "decio_jatai", "Decio_phylo_jatai.nex",     "Jatai, Luiz Antonio, SP",  "Interior SP",        "inferred",
  "decio_assis", "Decio_phylo.nex",           "Assis, SP",                "Interior SP",        "inferred",
  "dias",        "Dias2007.nex",              "Ubatuba?, SP",             "Atlantic Forest SE", "inferred",
  "gomes",       "Gomes.nex",                 "Campos do Jordao?, SP",    "Atlantic Forest SE", "inferred",
  "mainara",     "Mainara.nex",               "coastal SP/RJ",            "Atlantic Forest SE", "inferred",
  "motta",       "Motta_phylo.nex",           "NW Sao Paulo",             "Interior SP",        "author",
  "peltzer",     "Peltzer_phylo.nex",         "Santa Fe?, Argentina",     "South (Argentina)",  "inferred",
  "rogerio",     "Rogeio_phylo.nex",          "SW Goias",                 "Cerrado GO",         "inferred",
  "vasconcelos", "Vasconcelos2011.nex",       "NW Sao Paulo",             "Interior SP",        "author",
  "wild",        "Wild1996.nex",              "Peru",                     "Amazonia (Peru)",    "inferred"
)
read_labels <- function(path) {
  if (grepl("[.]nex$", path)) nexus_labels(path) else ape::read.tree(path)$tip.label
}
labels <- purrr::map2(datasets$dataset, file.path(raw_dir, datasets$file), \(d, f)
  tibble::tibble(dataset = d, label = read_labels(f))) |>
  dplyr::bind_rows() |>
  dplyr::left_join(dplyr::select(datasets, dataset, region), by = "dataset")

# Mainara: codes only; names from the clades of Mainara_phylo.txt and the
# coastal SP/RJ fauna (inferred), except where confirmed (2026-10-19), in a
# local table (data/raw/nri_pcps/mainara_codes.csv, not in git: unpublished)
mainara <- readr::read_csv(file.path(raw_dir, "mainara_codes.csv"), show_col_types = FALSE)

# names a code may stand for (R/taxonomy.R)
labels$is_full <- grepl("^[A-Z][a-z]+_[a-z]", labels$label)
full_names <- unique(gsub("_", " ", labels$label[labels$is_full]))
dictionary <- name_dictionary(full_names[grepl(binomial_re, full_names)])
labels$undetermined <- is_undetermined(labels$label)
labels <- labels |>
  dplyr::mutate(
    name_2012 = dplyr::case_when(
      undetermined ~ NA_character_,
      is_full ~ gsub("_", " ", label),
      dataset == "mainara" ~ mainara$name[match(label, mainara$label)],
      TRUE ~ vapply(label, resolve_code, character(1), dictionary = dictionary)
    ),
    resolution = dplyr::case_when(
      is_full ~ "full name",
      dataset == "mainara" ~ paste("Mainara table:", mainara$source[match(label, mainara$label)]),
      !is.na(name_2012) ~ "code resolved",
      TRUE ~ "code unresolved"
    )
  )

labels$asw_name <- unname(sync_asw(stats::na.omit(labels$name_2012))[labels$name_2012])

labels$corrected <- !labels$undetermined & !is.na(labels$name_2012) & is.na(labels$asw_name)
# names the ASW does not know (spelling errors of full names or of the 2012
# names a code resolved to): the closest ASW name within 2 edits, if unique
labels$asw_name[labels$corrected] <- vapply(labels$name_2012[labels$corrected], fix_spelling,
                                            character(1), dictionary = dictionary)

portik_tree <- ape::read.tree(file.path("data", "raw", "portik2023_timetree.tre"))
portik_asw <- sync_asw(gsub("_", " ", portik_tree$tip.label))
labels <- labels |>
  dplyr::mutate(
    portik_tip = portik_tree$tip.label[match(asw_name, portik_asw)],
    status = dplyr::case_when(
      undetermined ~ "excluded: undetermined",
      is.na(asw_name) ~ "excluded: unresolved",
      is.na(portik_tip) ~ "excluded: not in Portik tree",
      startsWith(resolution, "Mainara table: inferred") ~ "kept: inferred (Mainara)",
      corrected ~ "kept: name corrected",
      TRUE ~ "kept"
    )
  )
readr::write_csv(labels, file.path(out_dir, "names.csv"))
readr::write_csv(datasets, file.path(out_dir, "datasets.csv"))
print(dplyr::count(labels, status))
kept <- dplyr::filter(labels, startsWith(status, "kept"))
kept |>
  dplyr::group_by(region) |>
  dplyr::summarise(species = dplyr::n_distinct(asw_name)) |>
  print(n = Inf)
message(dplyr::n_distinct(kept$asw_name), " species kept")

# Step 2: nested pools. The bootstrap trees differ from the point-estimate tree
# in topology and root age, so their 2.5-97.5% quantiles describe the variation
# among bootstrap trees, not a confidence interval around the point estimate
pools <- list(
  "Atlantic Forest" = c("Atlantic Forest SE", "Atlantic Forest S"),
  "Atlantic Forest + interior SP" = c("Atlantic Forest SE", "Atlantic Forest S", "Interior SP"),
  "South America" = unique(datasets$region)
)
boot <- ape::read.tree(file.path("data", "raw", "portik2023_timetree_bootstraps.tre"))
pool_shape <- purrr::imap(pools, \(regions, pool) {
  tips <- unique(kept$portik_tip[kept$region %in% regions])
  stopifnot(!anyNA(tips))
  shape_of <- \(tr) tree_shape(ape::keep.tip(tr, tips))
  point <- shape_of(portik_tree)
  b <- vapply(boot, shape_of, numeric(2))
  summarise_boot <- \(v) stats::quantile(v, c(0.025, 0.5, 0.975), names = FALSE)
  g <- summarise_boot(b["gamma", ])
  cl <- summarise_boot(b["colless", ])
  tibble::tibble(pool = pool, n = length(tips),
                 gamma = point[["gamma"]], gamma_boot_median = g[2], gamma_boot_lo = g[1],
                 gamma_boot_hi = g[3], colless = point[["colless"]],
                 colless_boot_median = cl[2], colless_boot_lo = cl[1], colless_boot_hi = cl[3])
}) |>
  dplyr::bind_rows()
readr::write_csv(pool_shape, file.path(out_dir, "pool_shapes.csv"))
readr::write_csv(dplyr::distinct(dplyr::select(kept, dataset, region, asw_name, portik_tip)),
                 file.path(out_dir, "pool_species.csv"))
print(pool_shape, width = Inf)
