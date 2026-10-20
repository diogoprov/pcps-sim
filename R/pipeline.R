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
#' `optima_uniform`: as `full` in the power scenarios, with the optima given the
#' same evenly spaced marginal distribution in every tree (separates the effect
#' of tree shape on the optima from that on the analysis), with the extra axes.
#' `layerB`: layer B (process-based dynamics) on the four corner cells and the
#' centre of the grid, the four archetypes and the eight scenarios, incidence
#' only; the cell sets the (beta, c) of the pool tree (option B of 2026-10-18,
#' docs/MEMORIA_PROJETO.md), and the shape of the pruned tree is recorded.
#' `layerB_spatial`: layer B on a spatial landscape (10 x 10 torus, autocorrelated
#' gradient, exponential dispersal kernel; [spatial_landscape()]) at N = 100,
#' with incidence and abundance data from the same dynamics, the extra PCPS
#' axes, and torus translations instead of the site shuffle (reviewer check).
#' `smoke_B`: tiny layer-B run.
#' `extra_axes` defaults to FALSE in all other profiles.
#' @param name profile name; defaults to the environment variable PCPS_PROFILE.
design_profile <- function(name = Sys.getenv("PCPS_PROFILE", "pilot_glm")) {
  all_cells <- as.vector(outer(paste0("g_", c("yule", "mid", "low")),
                               paste0("c_", c("yule", "mid", "high")), paste, sep = "__"))
  full <- list(
    name = name, cells = all_cells,
    n_spp = c(50, 100, 200), scenarios = scenarios_layerA()$scenario,
    data_types = c("incidence", "abundance"), reps = 1000L, chunk = 50L,
    methods = c("pfw", "glm", "pcc"), nperm = 999L, nperm_pcc = 999L, k_trees = 1000L
  )
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
    full = full,
    axes = utils::modifyList(full, list(name = name, methods = c("pfw", "glm"),
                                        extra_axes = TRUE)),
    optima_uniform = utils::modifyList(full, list(name = name,
                                                  scenarios = c("S3a", "S3b", "S3sigma"),
                                                  optima = "uniform", extra_axes = TRUE)),
    layerB = utils::modifyList(full, list(
      name = name, layer = "B",
      cells = c("g_yule__c_yule", "g_low__c_yule", "g_yule__c_high", "g_low__c_high",
                "g_mid__c_mid"),
      archetypes = archetypes_layerB()$archetype, data_types = "incidence",
      extra_axes = TRUE, pool_file = layerB_pool_file()
    )),
    layerB_spatial = utils::modifyList(full, list(
      name = name, layer = "B",
      cells = c("g_yule__c_yule", "g_low__c_yule", "g_yule__c_high", "g_low__c_high",
                "g_mid__c_mid"),
      n_spp = 100, archetypes = archetypes_layerB()$archetype,
      data_types = c("incidence", "abundance"), extra_axes = TRUE, spatial = TRUE,
      pool_file = layerB_pool_file()
    )),
    smoke_B = list(
      name = name, layer = "B", cells = c("g_yule__c_yule", "g_low__c_high"),
      n_spp = 50, scenarios = c("S0", "S3b"),
      archetypes = c("sorting__phylo", "mass_effects__none"), data_types = "incidence",
      reps = 4L, chunk = 2L, methods = c("pfw", "glm", "pcc"), nperm = 19L,
      nperm_pcc = 19L, extra_axes = TRUE, pool_file = layerB_pool_file()
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

#' Replicate numbers of one chunk
chunk_reps <- function(chunk, profile) {
  seq((chunk - 1) * profile$chunk + 1, min(chunk * profile$chunk, profile$reps))
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
  reps <- chunk_reps(design_row$chunk, profile)
  t0 <- Sys.time()
  purrr::map(reps, \(r) {
    mc <- sim_metacommunity_A(
      trees[[r]], design_row$lambda, design_row$sigma_signal, design_row$h_signal,
      design_row$env_original, incidence = design_row$data_type == "incidence",
      optima = if (is.null(profile$optima)) "bm" else profile$optima
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

#' Layer-B feasible pools of the profile ([feasible_pools()])
pools_B <- function(profile, pool_file) {
  feasible_pools(pool_file) |>
    dplyr::filter(archetype %in% profile$archetypes, scenario %in% profile$scenarios,
                  target %in% profile$n_spp)
}

#' Shape targets of the profile cells for each layer-B pool size
#'
#' Tree shape statistics change with the number of tips for fixed (beta, c), so
#' each pool size gets its own calibration (option B of 2026-10-18,
#' docs/MEMORIA_PROJETO.md): a pool of `n_pool` tips for target N is drawn with
#' the gamma and Colless targets of the cell at N, as layer-A trees of N tips.
#' @param levels output of [shape_levels()]; `pools` output of [pools_B()].
#' @return one row per (target, n_pool, cell), with `n = n_pool`, ready for
#'   [calibrate_cells()] grouped by (target, n).
pool_levels_B <- function(levels, pools, profile) {
  levels |>
    dplyr::filter(paste0("g_", gamma_level, "__c_", colless_level) %in% profile$cells) |>
    dplyr::rename(target = n) |>
    dplyr::inner_join(dplyr::distinct(pools, target, n_pool), by = "target",
                      relationship = "many-to-many") |>
    dplyr::mutate(n = n_pool)
}

#' Layer-B design table: one row per (cell, archetype, scenario, target N,
#' chunk of replicates), with the archetype, scenario, calibrated pool size and
#' the (beta, c) and shape targets of the cell for that pool size.
#' @param pool_cells output of [calibrate_cells()] on [pool_levels_B()].
design_table_B <- function(profile, pools, pool_cells) {
  stopifnot(!anyNA(pool_cells$beta), !anyNA(pool_cells$c))
  pool_cells <- dplyr::select(pool_cells, target, n_pool, cell, beta, c, gamma_target,
                              colless_target)
  pools |>
    dplyr::inner_join(archetypes_layerB(), by = "archetype") |>
    dplyr::inner_join(scenarios_layerA(), by = "scenario") |>
    dplyr::inner_join(pool_cells, by = c("target", "n_pool"),
                      relationship = "many-to-many") |>
    tidyr::crossing(chunk = seq_len(ceiling(profile$reps / profile$chunk)))
}

#' Run one chunk of layer-B replicates
#'
#' Each replicate draws a new pool tree of the calibrated pool size within the
#' shape tolerance of its cell ([sim_tree_cell()], as in layer A), so the cell
#' describes the regional pool; the shapes of the pool tree and of the tree
#' pruned by the dynamics (the one analysed) are recorded as covariates. Up to
#' 5000 tries per pool tree, since acceptance is only ~8-12% in some cells.
#' @param design_row one row of [design_table_B()].
run_chunk_B <- function(design_row, profile) {
  reps <- chunk_reps(design_row$chunk, profile)
  space <- if (isTRUE(profile$spatial)) spatial_landscape() else NULL
  data_types <- if (is.null(profile$data_types)) "incidence" else profile$data_types
  t0 <- Sys.time()
  purrr::map(reps, \(r) {
    pool <- sim_tree_cell(1, design_row$n_pool, design_row$beta, design_row$c,
                          design_row$gamma_target, design_row$colless_target,
                          max_tries = 5000)
    stopifnot(length(pool) == 1)
    pool <- pool[[1]]
    # counts are simulated once; incidence is derived from them, so both data
    # types share the metacommunity (the spatial profile analyses both)
    mc <- sim_metacommunity_B(pool, design_row$lambda, design_row$sigma_signal,
                              design_row$h_signal, design_row$env_original,
                              dispersal = design_row$dispersal,
                              competition = design_row$competition,
                              incidence = identical(data_types, "incidence"), space = space)
    site_env <- if (is.null(space)) NULL else torus_site_env(mc, space)
    shape_pool <- tree_shape(pool)
    shape_obs <- tree_shape(mc$tree)
    purrr::map(data_types, \(dt) {
      mc_dt <- mc
      if (dt == "incidence") mc_dt$comm <- (mc$comm > 0) * 1
      analyse_metacommunity(mc_dt, profile$methods, profile$nperm, profile$nperm_pcc,
                            extra_axes = isTRUE(profile$extra_axes), site_env = site_env) |>
        dplyr::mutate(data_type = dt)
    }) |>
      dplyr::bind_rows() |>
      dplyr::mutate(rep = r, n_realized = ncol(mc$comm),
                    gamma_pool = shape_pool[["gamma"]], colless_pool = shape_pool[["colless"]],
                    gamma = shape_obs[["gamma"]], colless = shape_obs[["colless"]])
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(scenario = design_row$scenario, archetype = design_row$archetype,
                  cell = design_row$cell, n_spp = design_row$target,
                  n_pool = design_row$n_pool, .before = 1) |>
    dplyr::relocate(data_type, .after = n_pool) |>
    dplyr::mutate(chunk_seconds = as.numeric(Sys.time() - t0, units = "secs"))
}
