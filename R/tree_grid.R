# 3 x 3 grid of tree shape (gamma x Colless) for each pool size N.
# Levels: Yule reference (0, the published designs) plus empirical anuran quantiles
# (analysis/05_tree_shape_empirical.R). Gamma is controlled by c and Colless by
# beta (R/trees.R); each is calibrated separately, then trees are drawn with the
# calibrated pair and kept only if both statistics fall within a tolerance band.

#' Grid levels from the empirical clade quantiles
#' @param emp output of analysis/05 (columns source, n, gamma, colless).
#' @param n_levels pool sizes; clades are binned around each size.
shape_levels <- function(emp, n_levels = c(50, 100, 200)) {
  emp <- dplyr::filter(emp, source == "clade") |>
    dplyr::mutate(n_bin = n_levels[findInterval(n, c(0, 75, 150, Inf))])
  emp |>
    dplyr::summarise(
      gamma_mid = stats::median(gamma), gamma_low = stats::quantile(gamma, 0.1),
      colless_mid = stats::median(colless), colless_high = stats::quantile(colless, 0.9),
      k = dplyr::n(), .by = n_bin
    ) |>
    dplyr::transmute(
      n = n_bin, k,
      gamma = purrr::pmap(list(gamma_low, gamma_mid),
                          \(lo, mid) c(yule = 0, mid = mid, low = lo)),
      colless = purrr::pmap(list(colless_mid, colless_high),
                            \(mid, hi) c(yule = 0, mid = mid, high = hi))
    ) |>
    tidyr::unnest_longer(gamma, indices_to = "gamma_level") |>
    tidyr::unnest_longer(colless, indices_to = "colless_level") |>
    dplyr::rename(gamma_target = gamma, colless_target = colless) |>
    dplyr::arrange(n)
}

#' Mean statistic over a parameter grid (one shape component at a time)
#'
#' Gamma depends only on the intervals and Colless only on the topology, so each
#' is calibrated holding the other process at its Yule value.
calibration_curve <- function(n, values, what = c("gamma", "colless"), reps = 200) {
  what <- match.arg(what)
  stat <- vapply(values, \(v) mean(replicate(reps, {
    tr <- if (what == "gamma") sim_tree_shape(n, beta = 0, c = v) else
      sim_tree_shape(n, beta = v, c = 0)
    tree_shape(tr)[[what]]
  })), 0)
  tibble::tibble(n = n, what = what, value = values, stat = stat)
}

#' Parameter value whose mean statistic hits the target (linear interpolation;
#' the curves are monotone: gamma increases with c, Colless decreases with beta)
invert_curve <- function(curve, target) {
  o <- order(curve$stat)
  stats::approx(curve$stat[o], curve$value[o], xout = target, rule = 1)$y
}

#' Draw `k` trees for one grid cell by rejection within tolerance bands
#' @return list of `phylo`, with attributes `shape` (tibble) and `acceptance`.
sim_tree_cell <- function(k, n, beta, c, gamma_target, colless_target,
                          tol_gamma = 0.5, tol_colless = 0.5, max_tries = 50 * k) {
  out <- vector("list", k)
  shape <- matrix(NA_real_, k, 2, dimnames = list(NULL, c("gamma", "colless")))
  got <- 0
  tries <- 0
  while (got < k && tries < max_tries) {
    tries <- tries + 1
    tr <- sim_tree_shape(n, beta = beta, c = c)
    s <- tree_shape(tr)
    if (abs(s[["gamma"]] - gamma_target) <= tol_gamma &&
        abs(s[["colless"]] - colless_target) <= tol_colless) {
      got <- got + 1
      out[[got]] <- tr
      shape[got, ] <- s
    }
  }
  if (got < k) warning(sprintf("cell n=%d: %d of %d trees after %d tries", n, got, k, tries))
  out <- out[seq_len(got)]
  attr(out, "shape") <- tibble::as_tibble(shape[seq_len(got), , drop = FALSE])
  attr(out, "acceptance") <- got / tries
  out
}
