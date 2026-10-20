# Regional species pools from the author's empirical tadpole metacommunities
# (project "NRI vs PCPS", 2012; local copy in data/raw/nri_pcps/, not in git:
# compiled from the literature and from unpublished data shared with the
# author). Step 1, names:
#   - taxon labels of each metacommunity from its NEXUS file (2012 names);
#   - abbreviated labels ("Bokcir", "Hpu", "A.albosignatus", "Halbopunctatus")
#     resolved by genus prefix and species epithet prefix (the epithet stem,
#     without the last two letters, so that gender changes still match), first
#     against the full 2012 names of the other datasets and then against the
#     current and synonym names of the ASW (AmphiNom) for species recorded in
#     Brazil, Argentina, Paraguay, Uruguay, Bolivia or Peru (AmphibiaWeb ranges);
#     ambiguous or unresolved codes go to the review table with candidates;
#   - names harmonized with the Amphibian Species of the World through AmphiNom
#     (asw_sync), and matched to the tips of the Portik et al. 2023 tree, also
#     harmonized with asw_sync.
#   - the author's rules (2026-10-19): undetermined taxa (sp., aff., cf., gr.)
#     and codes that stay ambiguous are excluded; typos are corrected when a
#     single ASW name is close (agrep); doubtful resolutions are kept and
#     flagged until the author decides.
# Writes data/derived/pools/names.csv (every label with its resolution and
# status) and datasets.csv.
# Usage: Rscript analysis/21_regional_pools.R

raw_dir <- file.path("data", "raw", "nri_pcps")
out_dir <- file.path("data", "derived", "pools")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# taxon labels of a NEXUS file (TAXLABELS block)
nexus_labels <- function(path) {
  x <- paste(readLines(path, warn = FALSE, encoding = "latin1"), collapse = " ")
  block <- regmatches(x, regexpr("(?i)TAXLABELS(.*?);", x, perl = TRUE))
  lab <- strsplit(sub("(?i)TAXLABELS", "", sub(";$", "", block), perl = TRUE),
                  "[[:space:]]+")[[1]]
  unique(gsub("'", "", lab[lab != ""]))
}

# one NEXUS per metacommunity (Wild19961.nex duplicates Wild1996.nex), with its
# locality and region; `locality_source` says whether the author confirmed it
# or it is inferred from the file names, the coordinates file and the species
datasets <- tibble::tribble(
  ~dataset,      ~file,                       ~locality,                ~region,              ~locality_source,
  "bocaina",     "bocaina_tadpole.trees.nex", "Serra da Bocaina, SP/RJ", "Atlantic Forest SE", "inferred",
  "both",        "Both_tree.nex",             "Cacapava do Sul, RS",    "South (Pampa)",      "inferred",
  "costa",       "Costa_phylo.nex",           "Maranhao?",              "Northeast?",         "inferred",
  "crivelari",   "Crivelari_phylo.nex",       "Parana",                 "Atlantic Forest S",  "author",
  "decio_jatai", "Decio_phylo_jatai.nex",     "Jatai, Luiz Antonio, SP", "Interior SP",       "inferred",
  "decio_assis", "Decio_phylo.nex",           "Assis, SP",              "Interior SP",        "inferred",
  "dias",        "Dias2007.nex",              "Ubatuba?, SP",           "Atlantic Forest SE", "inferred",
  "gomes",       "Gomes.nex",                 "Campos do Jordao?, SP",  "Atlantic Forest SE", "inferred",
  "mainara",     "Mainara.nex",               "Mainara (-23.36, -44.85)", "Atlantic Forest SE", "inferred",
  "motta",       "Motta_phylo.nex",           "NW Sao Paulo",           "Interior SP",        "author",
  "peltzer",     "Peltzer_phylo.nex",         "Santa Fe?, Argentina",   "South (Argentina)",  "inferred",
  "rogerio",     "Rogeio_phylo.nex",          "SW Goias",               "Cerrado GO",         "inferred",
  "vasconcelos", "Vasconcelos2011.nex",       "NW Sao Paulo",           "Interior SP",        "author",
  "wild",        "Wild1996.nex",              "Peru",                   "Amazonia (Peru)",    "inferred"
)
labels <- purrr::map2(datasets$dataset, file.path(raw_dir, datasets$file), \(d, f)
  tibble::tibble(dataset = d, label = nexus_labels(f))) |>
  dplyr::bind_rows() |>
  dplyr::left_join(dplyr::select(datasets, dataset, region), by = "dataset")

# full names: "Genus_species" (possibly with aff./cf./sp.)
is_full <- grepl("^[A-Z][a-z]+_[a-z]", labels$label)
full_names <- unique(gsub("_", " ", labels$label[is_full]))
binomials <- full_names[grepl("^[A-Z][a-z]+ [a-z-]+$", full_names)]

# dictionary of names a code may stand for: 2012 binomials, and current and
# synonym ASW names of species in the study countries (synonyms map to the
# current name)
countries <- c("BR", "AR", "PY", "UY", "BO", "PE")
aweb <- AmphiNom::amphweb
local_species <- paste(aweb$genus, aweb$species)[
  vapply(strsplit(aweb$isocc, ","), \(x) any(x %in% countries), logical(1))]
syn <- AmphiNom::asw_synonyms
syn <- syn[syn$species %in% local_species, ]
dictionary <- dplyr::bind_rows(
  tibble::tibble(name = binomials, current = binomials, source = "2012 names"),
  tibble::tibble(name = syn$synonyms, current = syn$species, source = "ASW")
) |>
  dplyr::filter(grepl("^[A-Z][a-z]+ [a-z-]+$", name)) |>
  dplyr::distinct(name, .keep_all = TRUE)

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
# candidates (names in `dict`) for a code; the epithet prefix is compared on its
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
# point to a single species: first the 2012 names, then the valid ASW names
# (a code like "P.boiei" also matches old synonyms such as Phyllomedusa boiei),
# then all ASW names
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
candidates_text <- function(code) {
  hit <- code_candidates(code, dictionary)
  paste(utils::head(unique(hit$current), 6), collapse = "; ")
}

labels <- labels |>
  dplyr::mutate(name_2012 = ifelse(is_full, gsub("_", " ", label),
                                   vapply(label, resolve_code, character(1))),
                resolution = dplyr::case_when(is_full ~ "full name",
                                              !is.na(name_2012) ~ "code resolved",
                                              TRUE ~ "code unresolved"))

# harmonize with ASW (AmphiNom): a name listed as a synonym maps to its current
# name; otherwise asw_sync (NA when the name is not found)
all_synonyms <- AmphiNom::asw_synonyms
sync <- function(x) {
  x <- unique(x)
  out <- all_synonyms$species[match(x, all_synonyms$synonyms)]
  todo <- is.na(out)
  if (any(todo)) {
    s <- suppressMessages(AmphiNom::asw_sync(x[todo]))
    out[todo] <- ifelse(s$status == "name_not_found", NA_character_,
                        s$ASW_names)[match(x[todo], s$query)]
  }
  stats::setNames(out, x)
}
query <- unique(stats::na.omit(labels$name_2012))
asw <- sync(query)
labels$asw_name <- unname(asw[labels$name_2012])

portik <- ape::read.tree(file.path("data", "raw", "portik2023_timetree.tre"))$tip.label
portik_names <- gsub("_", " ", portik)
portik_asw <- sync(portik_names)
portik_species <- unique(stats::na.omit(portik_asw))
labels <- labels |>
  dplyr::mutate(in_portik = !is.na(asw_name) & asw_name %in% portik_species,
                portik_tip = portik[match(asw_name, portik_asw[portik_names])])

# undetermined taxa, including the codes of "Genus sp." ("Hsp", "Lsp", "Ssp")
undetermined <- grepl("(^|[ _.])(sp|aff|cf|gr)([ _.0-9]|$)", labels$label) |
  grepl("^[A-Z]sp[0-9]*$", labels$label)
# typos in full names: the closest ASW name (the 2012 names include the typos
# themselves), if it is within 2 edits and points to a single current name
asw_dictionary <- dictionary[dictionary$source == "ASW", ]
fix_typo <- function(name) {
  d <- utils::adist(name, asw_dictionary$name)[1, ]
  hit <- unique(asw_dictionary$current[d == min(d)])
  if (min(d) <= 2 && length(hit) == 1) hit else NA_character_
}
typo <- !undetermined & labels$resolution == "full name" & is.na(labels$asw_name)
labels$asw_name[typo] <- vapply(labels$name_2012[typo], fix_typo, character(1))
# doubtful resolutions, kept until the author decides
doubtful <- c(both = "Lma", mainara = "Stra", mainara = "Eova", crivelari = "V.uranoscopa")
labels <- labels |>
  dplyr::mutate(
    in_portik = !is.na(asw_name) & asw_name %in% portik_species,
    portik_tip = portik[match(asw_name, portik_asw[portik_names])],
    status = dplyr::case_when(
      undetermined ~ "excluded: undetermined",
      is.na(asw_name) ~ "excluded: unresolved",
      !in_portik ~ "excluded: not in Portik tree",
      paste(dataset, label) %in% paste(names(doubtful), doubtful) ~ "kept: doubtful",
      typo ~ "kept: typo corrected",
      TRUE ~ "kept"
    )
  )
readr::write_csv(labels, file.path(out_dir, "names.csv"))
print(dplyr::count(labels, status))
readr::write_csv(datasets, file.path(out_dir, "datasets.csv"))
labels |>
  dplyr::filter(startsWith(status, "kept")) |>
  dplyr::group_by(region) |>
  dplyr::summarise(species = dplyr::n_distinct(asw_name)) |>
  print(n = Inf)
message(dplyr::n_distinct(labels$asw_name[startsWith(labels$status, "kept")]),
        " species kept")

# Step 2: three nested regional pools, pruned from the Portik tree and from its
# 100 time-calibrated bootstrap trees; shape (gamma, Colless normalized against
# Yule, as in R/trees.R) of each, against the grid targets of the nearest N
pools <- list(
  "Atlantic Forest" = c("Atlantic Forest SE", "Atlantic Forest S"),
  "Atlantic Forest + interior SP" = c("Atlantic Forest SE", "Atlantic Forest S", "Interior SP"),
  "South America" = unique(datasets$region)
)
kept <- dplyr::filter(labels, startsWith(status, "kept"))
for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
boot <- ape::read.tree(file.path("data", "raw", "portik2023_timetree_bootstraps.tre"))
pool_shape <- purrr::imap(pools, \(regions, pool) {
  tips <- unique(kept$portik_tip[kept$region %in% regions])
  shape_of <- \(tr) tree_shape(ape::keep.tip(tr, tips))
  b <- vapply(boot, shape_of, numeric(2))
  tibble::tibble(pool = pool, n = length(tips),
                 gamma = shape_of(ape::read.tree(file.path("data", "raw",
                                                           "portik2023_timetree.tre")))[["gamma"]],
                 gamma_lo = stats::quantile(b["gamma", ], 0.025, names = FALSE),
                 gamma_hi = stats::quantile(b["gamma", ], 0.975, names = FALSE),
                 colless = shape_of(ape::read.tree(file.path("data", "raw",
                                                             "portik2023_timetree.tre")))[["colless"]],
                 colless_lo = stats::quantile(b["colless", ], 0.025, names = FALSE),
                 colless_hi = stats::quantile(b["colless", ], 0.975, names = FALSE))
}) |>
  dplyr::bind_rows()
readr::write_csv(pool_shape, file.path(out_dir, "pool_shapes.csv"))
readr::write_csv(dplyr::transmute(kept, dataset, region, asw_name, portik_tip) |>
                   dplyr::distinct(), file.path(out_dir, "pool_species.csv"))
print(pool_shape)
