# Regional species pools from the author's empirical tadpole metacommunities
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
#   - the author's rules (2026-10-19): undetermined taxa (sp., aff., cf., gr.,
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
binomial_re <- "^[A-Z][a-z]+ [a-z-]+$"

# taxon labels of a NEXUS file (TAXLABELS block)
nexus_labels <- function(path) {
  x <- paste(readLines(path, warn = FALSE, encoding = "latin1"), collapse = " ")
  block <- regmatches(x, regexpr("(?i)TAXLABELS(.*?);", x, perl = TRUE))
  lab <- strsplit(sub("(?i)TAXLABELS", "", sub(";$", "", block), perl = TRUE),
                  "[[:space:]]+")[[1]]
  unique(gsub("'", "", lab[lab != ""]))
}

# one file per metacommunity (Wild19961.nex duplicates Wild1996.nex), with its
# locality and region; `locality_source` says whether the author confirmed it
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
  "mainara",     "Mainara.nex",               "Mainara (-23.36, -44.85)", "Atlantic Forest SE", "inferred",
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
# coastal SP/RJ fauna (inferred), except where the author confirmed (2026-10-19)
mainara <- tibble::tribble(
  ~label,  ~name,                         ~source,
  "Hfab",  "Hypsiboas faber",             "inferred",
  "Halb",  "Hypsiboas albomarginatus",    "inferred",
  "Aeug",  "Aplastodiscus eugenioi",      "inferred",
  "Bhyl",  "Bokermannohyla hylax",        "inferred",
  "Tmes",  "Trachycephalus mesophaeus",   "inferred",
  "Dber",  "Dendropsophus berthalutzae",  "inferred",
  "Dmin",  "Dendropsophus minutus",       "inferred",
  "Dele",  "Dendropsophus elegans",       "inferred",
  "Stra",  "Scinax trapicheiroi",         "author (Ololygon trapicheiroi)",
  "Slit",  "Scinax littoralis",           "inferred",
  "Sarg",  "Scinax argyreornatus",        "inferred",
  "Sperp", "Scinax perpusillus",          "inferred",
  "Sper",  "Scinax perereca",             "inferred",
  "Shay",  "Scinax hayii",                "inferred",
  "Salt",  "Scinax alter",                "inferred",
  "Rorn",  "Rhinella ornata",             "inferred",
  "Llat",  "Leptodactylus latrans",       "inferred",
  "Patl",  "Physalaemus atlanticus",      "inferred",
  "Eova",  "Elachistocleis ovalis",       "author",
  "Ccar",  "Chiasmocleis carvalhoi",      "inferred"
)

# dictionary of names a code may stand for: 2012 binomials, and current and
# synonym ASW names of species in the study countries (synonyms map to the
# current name); a name may appear in both sources
labels$is_full <- grepl("^[A-Z][a-z]+_[a-z]", labels$label)
full_names <- unique(gsub("_", " ", labels$label[labels$is_full]))
binomials_2012 <- full_names[grepl(binomial_re, full_names)]
study_countries <- c("BR", "AR", "PY", "UY", "BO", "PE")
amphibiaweb <- AmphiNom::amphweb
local_species <- paste(amphibiaweb$genus, amphibiaweb$species)[
  vapply(strsplit(amphibiaweb$isocc, ","), \(x) any(x %in% study_countries), logical(1))]
asw_synonyms <- AmphiNom::asw_synonyms
local_synonyms <- asw_synonyms[asw_synonyms$species %in% local_species, ]
dictionary <- dplyr::bind_rows(
  tibble::tibble(name = binomials_2012, current = binomials_2012, source = "2012 names"),
  tibble::tibble(name = local_synonyms$synonyms, current = local_synonyms$species,
                 source = "ASW")
) |>
  dplyr::filter(grepl(binomial_re, name)) |>
  dplyr::distinct(name, source, .keep_all = TRUE)

# possible splits of a code into genus and epithet prefixes: "A.albosignatus";
# "Bokcir" (3 + 3), "Rorn" (1 + 3), "Hpu" (1 + 2), "Halbopunctatus" (1 +
# epithet); a 6-letter code may be 3 + 3 ("Bokcir") or 1 + epithet ("Dnanus")
code_splits <- function(code) {
  if (grepl("[.]", code)) {
    p <- strsplit(code, "[.]")[[1]]
    return(list(c(p[1], p[length(p)])))
  }
  split_at <- function(k) c(substr(code, 1, k), substr(code, k + 1, nchar(code)))
  if (nchar(code) == 6) list(split_at(3), split_at(1)) else list(split_at(1))
}
# candidates (rows of `dict`) for a code; the epithet prefix is compared on its
# stem, so "albopunctatus" also matches "albopunctata"
code_candidates <- function(code, dict) {
  genus <- tolower(sub(" .*", "", dict$name))
  epithet <- sub(".* ", "", dict$name)
  hits <- lapply(code_splits(code), \(parts) {
    parts <- tolower(parts)
    stem <- if (nchar(parts[2]) > 5) substr(parts[2], 1, nchar(parts[2]) - 2) else parts[2]
    dict[startsWith(genus, parts[1]) & startsWith(epithet, stem), ]
  })
  dplyr::distinct(dplyr::bind_rows(hits))
}
# a code is resolved when its candidates in the first dictionary that has any
# point to a single species: the 2012 names, then the valid ASW names (a code
# like "P.boiei" also matches old synonyms such as Phyllomedusa boiei), then
# all ASW names
resolve_code <- function(code) {
  sets <- list(dictionary[dictionary$source == "2012 names", ],
               dictionary[dictionary$source == "ASW" & dictionary$name == dictionary$current, ],
               dictionary[dictionary$source == "ASW", ])
  for (dict in sets) {
    hit <- code_candidates(code, dict)
    if (nrow(hit) > 0)
      return(if (length(unique(hit$current)) == 1) hit$name[1] else NA_character_)
  }
  NA_character_
}

# undetermined taxa, including the codes of "Genus sp." ("Hsp", "Lsp", "Ssp")
labels$undetermined <- grepl("(^|[ _.])(sp|aff|cf|gr)([ _.0-9]|$)", labels$label) |
  grepl("^[A-Z]sp[0-9]*$", labels$label)
labels <- labels |>
  dplyr::mutate(
    name_2012 = dplyr::case_when(
      undetermined ~ NA_character_,
      is_full ~ gsub("_", " ", label),
      dataset == "mainara" ~ mainara$name[match(label, mainara$label)],
      TRUE ~ vapply(label, resolve_code, character(1))
    ),
    resolution = dplyr::case_when(
      is_full ~ "full name",
      dataset == "mainara" ~ paste("Mainara table:", mainara$source[match(label, mainara$label)]),
      !is.na(name_2012) ~ "code resolved",
      TRUE ~ "code unresolved"
    )
  )

# harmonize with ASW (AmphiNom): a name listed as a synonym maps to its current
# name; otherwise asw_sync (NA when the name is not found)
sync <- function(x) {
  x <- unique(x)
  out <- asw_synonyms$species[match(x, asw_synonyms$synonyms)]
  todo <- is.na(out)
  if (any(todo)) {
    s <- suppressMessages(AmphiNom::asw_sync(x[todo]))
    out[todo] <- ifelse(s$status == "name_not_found", NA_character_,
                        s$ASW_names)[match(x[todo], s$query)]
  }
  stats::setNames(out, x)
}
labels$asw_name <- unname(sync(stats::na.omit(labels$name_2012))[labels$name_2012])

# names the ASW does not know (spelling errors of full names or of the 2012
# names a code resolved to): the closest ASW name within 2 edits, if unique
asw_dictionary <- dictionary[dictionary$source == "ASW", ]
fix_spelling <- function(name) {
  d <- utils::adist(name, asw_dictionary$name)[1, ]
  hit <- unique(asw_dictionary$current[d == min(d)])
  if (min(d) <= 2 && length(hit) == 1) hit else NA_character_
}
labels$corrected <- !labels$undetermined & !is.na(labels$name_2012) & is.na(labels$asw_name)
labels$asw_name[labels$corrected] <- vapply(labels$name_2012[labels$corrected], fix_spelling,
                                            character(1))

portik_tree <- ape::read.tree(file.path("data", "raw", "portik2023_timetree.tre"))
portik_asw <- sync(gsub("_", " ", portik_tree$tip.label))
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
