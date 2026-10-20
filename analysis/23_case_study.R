# Case study (Box 1): the tests and the clade-recovery methods applied to the
# empirical tadpole metacommunities that have at least 10 ponds and
# pond depth measured in every pond: Serra da Bocaina (Provete et al. 2014,
# Hydrobiologia), Both et al. (southern Brazil) and "vitao" (Prado &
# Rossa-Feres 2014, Austral Ecology); local copies in data/raw/nri_pcps/, not in
# git. Incidence data; depth as the common gradient, and canopy cover for
# Bocaina too, the gradient found to structure its phylogenetic composition in
# Provete (2018, bioRxiv 10.1101/252478). The tree is the Portik et al. 2023
# tree pruned to the species of each metacommunity (names from
# analysis/21_regional_pools.R for Bocaina and Both, resolved here for vitao).
# Writes data/derived/case_study/{tests,clades,metacommunities}.csv.
# Usage: Rscript analysis/23_case_study.R [nperm]

args <- commandArgs(trailingOnly = TRUE)
nperm <- if (length(args) >= 1) as.integer(args[[1]]) else 999L

for (f in list.files("R", "[.]R$", full.names = TRUE)) source(f)
raw_dir <- file.path("data", "raw", "nri_pcps")
out_dir <- file.path("data", "derived", "case_study")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20261019)

# tab-separated tables saved with CR line endings (Mac classic) and, for some,
# decimal commas
read_table_cr <- function(file, dec = ".") {
  x <- readLines(file.path(raw_dir, file), warn = FALSE, encoding = "latin1")
  x <- unlist(strsplit(x, "\r", fixed = TRUE))
  x <- x[nzchar(trimws(x))]
  utils::read.delim(text = x, row.names = 1, check.names = FALSE, dec = dec)
}

names_tab <- readr::read_csv(file.path("data", "derived", "pools", "names.csv"),
                             show_col_types = FALSE)
portik_tree <- ape::read.tree(file.path("data", "raw", "portik2023_timetree.tre"))
portik_asw <- sync_asw(gsub("_", " ", portik_tree$tip.label))

# Portik tip of each label: from the pools table, or resolved here for vitao
tips_from_pools <- function(dataset, labels) {
  tab <- names_tab[names_tab$dataset == dataset & startsWith(names_tab$status, "kept"), ]
  tab$portik_tip[match(labels, tab$label)]
}
tips_resolved <- function(labels) {
  binomials <- unique(stats::na.omit(names_tab$name_2012))
  dictionary <- name_dictionary(binomials[grepl(binomial_re, binomials)])
  name <- ifelse(is_undetermined(labels), NA_character_,
                 vapply(labels, resolve_code, character(1), dictionary = dictionary))
  asw <- unname(sync_asw(stats::na.omit(name))[name])
  fix <- !is.na(name) & is.na(asw)
  asw[fix] <- vapply(name[fix], fix_spelling, character(1), dictionary = dictionary)
  portik_tree$tip.label[match(asw, portik_asw)]
}

sets <- list(
  bocaina_depth = list(comm = t(read_table_cr("Bocaina_compos.txt")),
                       env = read_table_cr("Bocaina_envir.txt")[, "dep"],
                       env_sites = rownames(read_table_cr("Bocaina_envir.txt")),
                       tips = \(lab) tips_from_pools("bocaina", lab)),
  bocaina_canopy = list(comm = t(read_table_cr("Bocaina_compos.txt")),
                        env = read_table_cr("Bocaina_envir.txt")[, "canopy"],
                        env_sites = rownames(read_table_cr("Bocaina_envir.txt")),
                        tips = \(lab) tips_from_pools("bocaina", lab)),
  both_depth = list(comm = as.matrix(read_table_cr("both_compos.txt")),
              env = read_table_cr("Both_envir.txt", dec = ",")[, "Mean_depth"],
              env_sites = rownames(read_table_cr("Both_envir.txt", dec = ",")),
              tips = \(lab) tips_from_pools("both", lab)),
  vitao_depth = list(comm = as.matrix(read_table_cr("vitao_compos.txt")),
               env = read_table_cr("vitao_env.txt")[, "depth(cm)"],
               env_sites = rownames(read_table_cr("vitao_env.txt")),
               tips = tips_resolved)
)

results <- purrr::imap(sets, \(s, name) {
  comm <- s$comm
  env <- stats::setNames(s$env, s$env_sites)[rownames(comm)]
  stopifnot(!anyNA(env))
  tips <- s$tips(colnames(comm))
  keep <- !is.na(tips) & !duplicated(tips)
  comm <- (comm[, keep, drop = FALSE] > 0) * 1
  colnames(comm) <- tips[keep]
  comm <- comm[rowSums(comm) > 0, colSums(comm) > 0, drop = FALSE]
  env <- env[rownames(comm)]
  tree <- ape::keep.tip(portik_tree, colnames(comm))
  comm <- comm[, tree$tip.label]
  shape <- tree_shape(tree)
  info <- tibble::tibble(dataset = name, sites = nrow(comm), species = ncol(comm),
                         species_dropped = sum(!keep), gamma = shape[["gamma"]],
                         colless = shape[["colless"]], gradient_min = min(env),
                         gradient_max = max(env))
  tests <- dplyr::bind_rows(test_dp_all(comm, tree, env, axes = 1:2, nperm = nperm,
                                        extra_axes = TRUE),
                            test_pcc(comm, tree, env, nperm = nperm)) |>
    dplyr::mutate(dataset = name, .before = 1)
  # best-scoring clade of each method (no planted clade here: which clade each
  # method associates with depth)
  scores <- list(
    pcps = node_scores(tree, pcps_species_scores(comm, tree, env)),
    pcc = node_scores(tree, pcc_species_scores(comm, tree, env)),
    evopca = node_scores(tree, evopca_species_scores(comm, tree, env)),
    phylofactor = phylofactor_node_scores(comm, tree, env)
  )
  clades <- purrr::imap(scores, \(sc, method) {
    best <- sc$node[which.max(sc$score)]
    members <- ape::extract.clade(tree, best)$tip.label
    tibble::tibble(dataset = name, method = method, clade_size = length(members),
                   genera = paste(sort(unique(sub("_.*", "", members))), collapse = ", "),
                   species = paste(members, collapse = "; "))
  }) |>
    dplyr::bind_rows()
  list(info = info, tests = tests, clades = clades)
})

info <- dplyr::bind_rows(lapply(results, `[[`, "info"))
tests <- dplyr::bind_rows(lapply(results, `[[`, "tests"))
clades <- dplyr::bind_rows(lapply(results, `[[`, "clades"))
readr::write_csv(info, file.path(out_dir, "metacommunities.csv"))
readr::write_csv(tests, file.path(out_dir, "tests.csv"))
readr::write_csv(clades, file.path(out_dir, "clades.csv"))
print(info, width = Inf)
print(dplyr::select(tests, dataset, test, stat, p_site, p_taxa, p_max), n = Inf)
print(dplyr::select(clades, dataset, method, clade_size, genera), n = Inf)
