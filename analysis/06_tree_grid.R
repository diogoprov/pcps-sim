# Tree-shape grid: calibrate c (gamma) and beta (Colless) for each N and draw
# trees for each of the 9 cells.
# Usage: Rscript analysis/06_tree_grid.R <trees_per_cell> <workers>

args <- commandArgs(trailingOnly = TRUE)
k <- if (length(args) >= 1) as.integer(args[[1]]) else 200L
workers <- if (length(args) >= 2) as.integer(args[[2]]) else 12L

for (f in list.files("R", full.names = TRUE)) source(f)
future::plan(future::multisession, workers = workers)
opts <- furrr::furrr_options(seed = 20261017L, globals = TRUE)
out_dir <- file.path("data", "derived", "trees")

emp <- readr::read_csv(file.path(out_dir, "empirical_shape_portik2023.csv"),
                       show_col_types = FALSE)
levels <- shape_levels(emp)

# Calibration ---------------------------------------------------------------
t0 <- Sys.time()
cal_grid <- dplyr::bind_rows(
  tidyr::crossing(n = c(50, 100, 200), what = "gamma", value = seq(-3, 0.2, 0.05)),
  tidyr::crossing(n = c(50, 100, 200), what = "colless", value = seq(-1.6, 0.4, 0.05))
)
curves <- furrr::future_pmap(cal_grid, \(n, what, value)
  calibration_curve(n, value, what, reps = 200), .options = opts) |>
  dplyr::bind_rows()
message("Calibration: ", format(Sys.time() - t0))
readr::write_csv(curves, file.path(out_dir, "calibration_curves.csv"))

cells <- levels |>
  dplyr::mutate(
    c = purrr::map2_dbl(n, gamma_target, \(nn, g)
      if (g == 0) 0 else invert_curve(dplyr::filter(curves, n == nn, what == "gamma"), g)),
    beta = purrr::map2_dbl(n, colless_target, \(nn, cl)
      if (cl == 0) 0 else invert_curve(dplyr::filter(curves, n == nn, what == "colless"), cl)),
    cell = paste0("g_", gamma_level, "__c_", colless_level)
  )
print(cells, width = Inf)
stopifnot(!anyNA(cells$c), !anyNA(cells$beta))

# Trees per cell --------------------------------------------------------------
# Split each cell into chunks so the 27 cells balance across workers.
chunk <- 25L
jobs <- tidyr::crossing(dplyr::select(cells, n, cell, beta, c, gamma_target, colless_target),
                        part = seq_len(ceiling(k / chunk)))
t0 <- Sys.time()
drawn <- furrr::future_pmap(jobs, \(n, cell, beta, c, gamma_target, colless_target, part)
  sim_tree_cell(chunk, n, beta, c, gamma_target, colless_target), .options = opts)
message("Trees: ", format(Sys.time() - t0))

jobs$trees <- drawn
grid <- jobs |>
  dplyr::summarise(
    trees = list(purrr::list_flatten(trees)),
    acceptance = mean(purrr::map_dbl(drawn[dplyr::cur_group_rows()], attr, "acceptance")),
    shape = list(dplyr::bind_rows(purrr::map(drawn[dplyr::cur_group_rows()], attr, "shape"))),
    .by = c(n, cell, beta, c, gamma_target, colless_target)
  )
saveRDS(grid, file.path(out_dir, sprintf("tree_grid_k%d.rds", k)))

grid |>
  dplyr::mutate(
    n_trees = lengths(trees),
    gamma_mean = purrr::map_dbl(shape, \(s) mean(s$gamma)),
    colless_mean = purrr::map_dbl(shape, \(s) mean(s$colless)),
    r_gamma_colless = purrr::map_dbl(shape, \(s) stats::cor(s$gamma, s$colless))
  ) |>
  dplyr::select(-trees, -shape) |>
  print(n = Inf, width = Inf)
