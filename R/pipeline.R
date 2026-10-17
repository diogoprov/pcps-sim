# Functions used by the `targets` pipeline (_targets.R): design profiles, tree
# grid construction and the replicate chunks.

#' Design profiles
#'
#' `smoke`: tiny run to check the pipeline end to end.
#' `bench`: time per replicate on a new machine (e.g. a small EC2 instance).
#' `pilot_glm`: mini-pilot of the D_P tests (ADONIS, Mantel, GLM on PCPS axes)
#' in three grid cells, incidence and abundance.
#' `full`: objective 1, all cells, scenarios and data types (docs/PLANO_SIMULACAO.md).
#' `axes`: as `full`, without PCC, adding the GLM with Procrustes and the RDA on
#' PCPS axes 1-2 (diagnosis of S2sigma, analysis/09_diag_s2sigma.R).
#' `extra_axes` defaults to FALSE in all other profiles.
#' @param name profile name; defaults to the environment variable PCPS_PROFILE.
design_profile <- function(name = Sys.getenv("PCPS_PROFILE", "pilot_glm")) {
  all_cells <- as.vector(outer(paste0("g_", c("yule", "mid", "low")),
                               paste0("c_", c("yule", "mid", "high")), paste, sep = "__"))
  switch(name,
    smoke = list(
      name = name, cells = c("g_yule__c_yule", "g_low__c_high"), n_spp = 50,
      scenarios = c("S0", "S3b"), data_types = c("incidence", "abundance"),
      reps = 4L, chunk = 2L, methods = c("pfw", "glm", "pcc"), nperm = 19L,
      nperm_pcc = 19L, k_trees = 4L
    ),
    bench = list(
      name = name, cells = "g_mid__c_mid", n_spp = c(50, 100, 200),
      scenarios = "S3b", data_types = "incidence", reps = 32L, chunk = 2L,
      methods = c("pfw", "glm", "pcc"), nperm = 999L, nperm_pcc = 999L, k_trees = 32L
    ),
    pilot_glm = list(
      name = name, cells = c("g_yule__c_yule", "g_low__c_yule", "g_yule__c_high"),
      n_spp = c(50, 100, 200), scenarios = scenarios_layerA()$scenario,
      data_types = c("incidence", "abundance"), reps = 100L, chunk = 25L,
      methods = c("pfw", "glm"), nperm = 999L, nperm_pcc = 999L, k_trees = 200L
    ),
    full = list(
      name = name, cells = all_cells,
      n_spp = c(50, 100, 200), scenarios = scenarios_layerA()$scenario,
      data_types = c("incidence", "abundance"), reps = 1000L, chunk = 50L,
      methods = c("pfw", "glm", "pcc"), nperm = 999L, nperm_pcc = 999L, k_trees = 1000L
    ),
    axes = list(
      name = name, cells = all_cells,
      n_spp = c(50, 100, 200), scenarios = scenarios_layerA()$scenario,
      data_types = c("incidence", "abundance"), reps = 1000L, chunk = 50L,
      methods = c("pfw", "glm"), nperm = 999L, nperm_pcc = 999L, k_trees = 1000L,
      extra_axes = TRUE
    ),
    stop("unknown profile: ", name)
  )
}

#' Design table without the tree dimension: one row per (scenario, data type,
#' chunk of replicates). The pipeline crosses it with the grid cells, so each
#' branch loads only the trees of its own cell.
design_table <- function(profile) {
  tidyr::crossing(
    scenario = profile$scenarios, data_type = profile$data_types,
    chunk = seq_len(ceiling(profile$reps / profile$chunk))
  ) |>
    dplyr::left_join(scenarios_layerA(), by = "scenario")
}

#' Empirical shape reference (analysis/05_tree_shape_empirical.R, clades only)
empirical_clade_shapes <- function(path) {
  tr <- ape::read.tree(path)
  genus <- sub("_.*", "", tr$tip.label)
  anura <- ape::extract.clade(
    tr, ape::getMRCA(tr, tr$tip.label[genus %in% c("Ascaphus", "Rana")])
  )
  sizes <- vapply(ape::prop.part(anura), length, 0)
  nodes <- ape::Ntip(anura) + which(sizes >= 40 & sizes <= 250)
  purrr::map(nodes, \(v) {
    t <- ape::extract.clade(anura, v)
    t$edge.length <- t$edge.length / max(ape::node.depth.edgelength(t))
    t <- phytools::force.ultrametric(t, method = "extend", message = FALSE)
    tibble::tibble(source = "clade", node = v, n = ape::Ntip(t),
                   gamma = treestats::gamma_statistic(t),
                   colless = treestats::colless(t, normalization = "yule"))
  }) |>
    dplyr::bind_rows()
}

#' Calibrated (beta, c) for each grid cell of one pool size
#' (analysis/06_tree_grid.R)
calibrate_cells <- function(levels_n, reps = 200) {
  n <- levels_n$n[1]
  gam <- calibration_curve(n, seq(-3, 0.2, 0.05), "gamma", reps)
  col <- calibration_curve(n, seq(-1.6, 0.4, 0.05), "colless", reps)
  levels_n |>
    dplyr::mutate(
      c = vapply(gamma_target, \(g) if (g == 0) 0 else invert_curve(gam, g), 0),
      beta = vapply(colless_target, \(cl) if (cl == 0) 0 else invert_curve(col, cl), 0),
      cell = paste0("g_", gamma_level, "__c_", colless_level)
    )
}

#' Trees for one grid cell
grid_cell_trees <- function(cell_row, k) {
  tr <- sim_tree_cell(k, cell_row$n, cell_row$beta, cell_row$c,
                      cell_row$gamma_target, cell_row$colless_target)
  stopifnot(length(tr) == k)
  tibble::tibble(n = cell_row$n, cell = cell_row$cell, trees = list(tr),
                 acceptance = attr(tr, "acceptance"))
}

#' Run one chunk of replicates
#'
#' Replicate r uses tree r of its cell, so the same tree is shared across
#' scenarios and data types (paired design). The random stream of each chunk is
#' set by `targets` (one seed per branch).
#' @param design_row one row of [design_table()].
#' @param cell_trees one row of the tree grid (n, cell, trees).
run_chunk <- function(design_row, cell_trees, profile) {
  trees <- cell_trees$trees[[1]]
  reps <- seq((design_row$chunk - 1) * profile$chunk + 1,
              min(design_row$chunk * profile$chunk, profile$reps))
  t0 <- Sys.time()
  purrr::map(reps, \(r) {
    mc <- sim_metacommunity_A(
      trees[[r]], design_row$lambda, design_row$sigma_signal, design_row$h_signal,
      design_row$env_original, incidence = design_row$data_type == "incidence"
    )
    analyse_metacommunity(mc, profile$methods, profile$nperm, profile$nperm_pcc,
                          extra_axes = isTRUE(profile$extra_axes)) |>
      dplyr::mutate(rep = r)
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(scenario = design_row$scenario, cell = cell_trees$cell,
                  n_spp = cell_trees$n, data_type = design_row$data_type,
                  .before = 1) |>
    dplyr::mutate(chunk_seconds = as.numeric(Sys.time() - t0, units = "secs"))
}
