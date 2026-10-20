# Species names of the empirical metacommunities: labels in the 2012 files
# (full names or abbreviated codes), harmonized with the Amphibian Species of
# the World (ASW) through AmphiNom. Used by analysis/21_regional_pools.R and
# analysis/23_case_study.R.

binomial_re <- "^[A-Z][a-z]+ [a-z-]+$"

#' Taxon labels of a NEXUS file (TAXLABELS block)
nexus_labels <- function(path) {
  x <- paste(readLines(path, warn = FALSE, encoding = "latin1"), collapse = " ")
  block <- regmatches(x, regexpr("(?i)TAXLABELS(.*?);", x, perl = TRUE))
  lab <- strsplit(sub("(?i)TAXLABELS", "", sub(";$", "", block), perl = TRUE),
                  "[[:space:]]+")[[1]]
  unique(gsub("'", "", lab[lab != ""]))
}

#' Names a code may stand for
#'
#' The 2012 binomials of the datasets, and the current and synonym ASW names of
#' species recorded in the study countries (AmphibiaWeb ranges); synonyms map to
#' their current name. A name may appear in both sources.
#' @param binomials_2012 full 2012 names ("Genus species").
#' @param countries ISO codes of the study countries.
#' @return tibble with `name`, `current` and `source` ("2012 names" or "ASW").
name_dictionary <- function(binomials_2012,
                            countries = c("BR", "AR", "PY", "UY", "BO", "PE")) {
  aweb <- AmphiNom::amphweb
  local_species <- paste(aweb$genus, aweb$species)[
    vapply(strsplit(aweb$isocc, ","), \(x) any(x %in% countries), logical(1))]
  syn <- AmphiNom::asw_synonyms
  syn <- syn[syn$species %in% local_species, ]
  dplyr::bind_rows(
    tibble::tibble(name = binomials_2012, current = binomials_2012, source = "2012 names"),
    tibble::tibble(name = syn$synonyms, current = syn$species, source = "ASW")
  ) |>
    dplyr::filter(grepl(binomial_re, name)) |>
    dplyr::distinct(name, source, .keep_all = TRUE)
}

#' Possible splits of a code into genus and epithet prefixes
#'
#' "A.albosignatus"; "Bokcir" (3 + 3), "Rorn" (1 + 3), "Hpu" (1 + 2),
#' "Halbopunctatus" (1 + epithet); a 6-letter code may be 3 + 3 ("Bokcir") or
#' 1 + epithet ("Dnanus").
code_splits <- function(code) {
  if (grepl("[.]", code)) {
    p <- strsplit(code, "[.]")[[1]]
    return(list(c(p[1], p[length(p)])))
  }
  split_at <- function(k) c(substr(code, 1, k), substr(code, k + 1, nchar(code)))
  if (nchar(code) == 6) list(split_at(3), split_at(1)) else list(split_at(1))
}

#' Rows of `dict` a code may stand for; the epithet prefix is compared on its
#' stem (without the last two letters), so gender changes still match
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

#' Resolve a code to a name
#'
#' Resolved when its candidates in the first set that has any point to a
#' single species: the 2012 names, then the valid ASW names (a code like
#' "P.boiei" also matches old synonyms such as Phyllomedusa boiei), then all
#' ASW names. NA otherwise.
resolve_code <- function(code, dictionary) {
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

#' Current ASW names: a name listed as a synonym maps to its current name;
#' otherwise asw_sync (NA when the name is not found)
#' @return named vector (names = unique queries).
sync_asw <- function(x) {
  syn <- AmphiNom::asw_synonyms
  x <- unique(x)
  out <- syn$species[match(x, syn$synonyms)]
  todo <- is.na(out)
  if (any(todo)) {
    s <- suppressMessages(AmphiNom::asw_sync(x[todo]))
    out[todo] <- ifelse(s$status == "name_not_found", NA_character_,
                        s$ASW_names)[match(x[todo], s$query)]
  }
  stats::setNames(out, x)
}

#' Closest ASW name to a name the ASW does not know (spelling errors), if
#' within 2 edits and pointing to a single current name; NA otherwise
fix_spelling <- function(name, dictionary) {
  asw <- dictionary[dictionary$source == "ASW", ]
  d <- utils::adist(name, asw$name)[1, ]
  hit <- unique(asw$current[d == min(d)])
  if (min(d) <= 2 && length(hit) == 1) hit else NA_character_
}

#' Undetermined taxa (sp., aff., cf., gr., and codes of "Genus sp.")
is_undetermined <- function(label) {
  grepl("(^|[ _.])(sp|aff|cf|gr)([ _.0-9]|$)", label) | grepl("^[A-Z]sp[0-9]*$", label)
}
