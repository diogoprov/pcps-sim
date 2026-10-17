# Simulation pipeline (layer A). Profile and workers come from the environment:
#   PCPS_PROFILE = pilot_glm | full   (R/pipeline.R, design_profile())
#   PCPS_WORKERS = number of crew workers (default: cores - 1)
# Run with analysis/run_pipeline.R, which also selects a store per profile.

library(targets)
library(tarchetypes)
tar_source("R")

profile <- design_profile()
workers <- as.integer(Sys.getenv("PCPS_WORKERS", parallel::detectCores() - 1))

tar_option_set(
  packages = c("dplyr", "tibble"),
  format = "qs",
  seed = 20261018,
  controller = crew::crew_controller_local(workers = workers, seconds_idle = 60),
  storage = "worker",
  retrieval = "worker",
  memory = "transient",
  garbage_collection = TRUE,
  error = "continue"
)

list(
  tar_target(portik_file, "data/raw/portik2023_timetree.tre", format = "file"),
  tar_target(empirical_shapes, empirical_clade_shapes(portik_file)),
  tar_target(
    levels_by_n,
    shape_levels(empirical_shapes) |>
      dplyr::filter(n %in% profile$n_spp) |>
      dplyr::group_by(n) |>
      targets::tar_group(),
    iteration = "group"
  ),
  tar_target(cells_calibrated, calibrate_cells(levels_by_n), pattern = map(levels_by_n)),
  tar_target(
    cells_used,
    cells_calibrated |>
      dplyr::filter(cell %in% profile$cells) |>
      dplyr::group_by(n, cell) |>
      targets::tar_group(),
    iteration = "group"
  ),
  tar_target(tree_grid, grid_cell_trees(cells_used, profile$k_trees),
             pattern = map(cells_used)),
  tar_target(
    design,
    design_table(profile) |>
      dplyr::group_by(scenario, data_type, chunk) |>
      targets::tar_group(),
    iteration = "group"
  ),
  tar_target(results, run_chunk(design, tree_grid, profile),
             pattern = cross(tree_grid, design)),
  tar_target(
    results_csv,
    {
      path <- file.path("data", "derived", profile$name, "results.csv")
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(results, path)
      path
    },
    format = "file"
  ),
  tar_target(
    rejection,
    rejection_rates(results, test, scenario, cell, n_spp, data_type)
  )
)
