# Layer A: static metacommunity simulation (Dray & Legendre 2008;
# Peres-Neto, Leibold & Dray 2012), as used by Duarte et al. 2016 and
# Peres-Neto & Kembel 2015.

#' Rescale a vector to [lo, hi]
rescale_to <- function(x, lo = -1, hi = 101) {
  lo + (x - min(x)) / (max(x) - min(x)) * (hi - lo)
}

#' Pure-birth tree as in Duarte et al. 2016 (Appendix S2)
#'
#' `geiger::sim.bdtree(b = 0.1, d = 0)`; zero-length edges set to 0.01.
sim_tree_duarte <- function(n_spp) {
  tree <- geiger::sim.bdtree(b = 0.1, d = 0, stop = "taxa", n = n_spp)
  tree$edge.length[tree$edge.length == 0] <- 0.01
  tree
}

#' Niche optima with Grafen's rho transformation (Duarte et al. 2016)
#'
#' Branch lengths are replaced by Grafen's (topology-based) lengths raised to
#' `rho`; rho = 1 gives Brownian-like signal, rho = 1e-4 almost none.
sim_optima_grafen <- function(tree, rho) {
  tr <- ape::compute.brlen(tree, method = "Grafen", power = rho)
  rescale_to(ape::rTraitCont(tr, model = "BM"))
}

#' Niche optima as in Peres-Neto & Kembel 2015 (Appendix A)
#'
#' Brownian trait (`phytools::fastBM`) if `signal`, otherwise normal deviates.
sim_optima_pnk <- function(tree, signal) {
  x <- if (signal) phytools::fastBM(tree) else
    stats::setNames(stats::rnorm(ape::Ntip(tree)), tree$tip.label)
  rescale_to(x)
}

#' Unimodal species responses along one gradient, with Poisson noise
#'
#' L_ij = h_j exp(-(E_i - u_j)^2 / (2 sigma_j^2)); h ~ U(0, 30),
#' sigma ~ N(mu_tol, 10). Neither source paper says how negative sigma is
#' handled; we use |sigma| (docs/ESPEC_PCC.md, section 4).
#' @param optima named vector of species optima (names = tip labels).
#' @param env environmental values of the sites.
#' @param incidence return presence-absence instead of abundances.
sim_comm_static <- function(optima, env, mu_tol = 10, incidence = TRUE) {
  k <- length(optima)
  h <- stats::runif(k, 0, 30)
  sigma <- abs(stats::rnorm(k, mu_tol, 10))
  lambda <- sweep(outer(env, optima, "-")^2, 2, 2 * sigma^2, "/")
  lambda <- sweep(exp(-lambda), 2, h, "*")
  comm <- matrix(stats::rpois(length(lambda), lambda), nrow = length(env),
                 dimnames = list(NULL, names(optima)))
  if (incidence) comm <- (comm > 0) * 1
  comm
}

#' Drop empty sites (P and the PCC are undefined for them)
drop_empty_sites <- function(comm, env) {
  keep <- rowSums(comm) > 0
  list(comm = comm[keep, , drop = FALSE], env = env[keep], n_dropped = sum(!keep))
}
