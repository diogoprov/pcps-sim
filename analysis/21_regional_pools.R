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
# Writes data/derived/pools/names.csv (every label with its resolution) and
# names_to_review.csv (what needs the author's taxonomic decision).
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

# one NEXUS per metacommunity (Wild19961.nex duplicates Wild1996.nex)
datasets <- tibble::tribble(
  ~dataset,      ~file,
  "bocaina",     "bocaina_tadpole.trees.nex",
  "both",        "Both_tree.nex",
  "costa",       "Costa_phylo.nex",
  "crivelari",   "Crivelari_phylo.nex",
  "decio_jatai", "Decio_phylo_jatai.nex",
  "decio_assis", "Decio_phylo.nex",
  "dias",        "Dias2007.nex",
  "gomes",       "Gomes.nex",
  "mainara",     "Mainara.nex",
  "motta",       "Motta_phylo.nex",
  "peltzer",     "Peltzer_phylo.nex",
  "rogerio",     "Rogeio_phylo.nex",
  "vasconcelos", "Vasconcelos2011.nex",
  "wild",        "Wild1996.nex"
)
labels <- purrr::map2(datasets$dataset, file.path(raw_dir, datasets$file), \(d, f)
  tibble::tibble(dataset = d, label = nexus_labels(f))) |>
  dplyr::bind_rows()

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

# split a code into genus and epithet prefixes
code_parts <- function(code) {
  if (grepl("[.]", code)) {
    p <- strsplit(code, "[.]")[[1]]
    return(c(p[1], p[length(p)]))
  }
  # "Bokcir" (3 + 3), "Rorn" (1 + 3), "Hpu" (1 + 2), "Halbopunctatus" (1 + epithet)
  k <- if (nchar(code) == 6) 3 else 1
  c(substr(code, 1, k), substr(code, k + 1, nchar(code)))
}
# candidates (names in `dict`) for a code; the epithet prefix is compared on its
# stem, so "albopunctatus" also matches "albopunctata"
code_candidates <- function(code, dict) {
  parts <- tolower(code_parts(code))
  stem <- if (nchar(parts[2]) > 5) substr(parts[2], 1, nchar(parts[2]) - 2) else parts[2]
  genus <- tolower(sub(" .*", "", dict$name))
  epithet <- sub(".* ", "", dict$name)
  dict[startsWith(genus, parts[1]) & startsWith(epithet, stem), ]
}
# a code is resolved when all its candidates in the first dictionary that has
# any point to the same species
resolve_code <- function(code) {
  for (src in c("2012 names", "ASW")) {
    hit <- code_candidates(code, dictionary[dictionary$source == src, ])
    if (nrow(hit) > 0) {
      out <- unique(hit$name)
      return(if (length(unique(hit$current)) == 1) out[1] else NA_character_)
    }
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

# harmonize with ASW (AmphiNom); asw_sync keeps the query order
sync <- function(x) {
  s <- suppressMessages(AmphiNom::asw_sync(unique(x)))
  stats::setNames(s$ASW_names, s$query)
}
query <- unique(stats::na.omit(labels$name_2012))
asw <- sync(query)
labels$asw_name <- unname(asw[labels$name_2012])

portik <- ape::read.tree(file.path("data", "raw", "portik2023_timetree.tre"))$tip.label
portik_names <- gsub("_", " ", portik)
portik_asw <- sync(portik_names)
in_portik <- unique(stats::na.omit(portik_asw))
labels <- labels |>
  dplyr::mutate(in_portik = !is.na(asw_name) & asw_name %in% in_portik,
                portik_tip = portik[match(asw_name, portik_asw[portik_names])])
readr::write_csv(labels, file.path(out_dir, "names.csv"))

review <- labels |>
  dplyr::filter(!in_portik) |>
  dplyr::mutate(candidates = ifelse(resolution == "full name", "",
                                    vapply(label, candidates_text, character(1))),
                decision = "")
readr::write_csv(review, file.path(out_dir, "names_to_review.csv"))
labels |>
  dplyr::count(dataset, resolution, in_portik) |>
  print(n = Inf)
message(nrow(review), " of ", nrow(labels), " labels need review; ",
        dplyr::n_distinct(labels$asw_name[labels$in_portik]), " species matched to Portik")
