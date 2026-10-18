# Timing of the main computational steps under the BLAS in use (run once per
# container or BLAS setting; the printed statistics show whether the setting
# changes the results). Used on 2026-10-18 to compare rocker/r-ver's OpenBLAS
# with and without OPENBLAS_NUM_THREADS=1.
# Usage: Rscript analysis/12_bench_blas.R <reps>

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 3L

for (f in list.files("R", full.names = TRUE)) source(f)
message("BLAS: ", extSoftVersion()[["BLAS"]], " | LAPACK: ", La_library())

time_it <- function(label, expr_fn) {
  set.seed(20261018)
  t0 <- Sys.time()
  out <- lapply(seq_len(reps), \(i) expr_fn())
  secs <- as.numeric(Sys.time() - t0, units = "secs") / reps
  message(sprintf("%-28s %7.2f s per replicate", label, secs))
  out
}

set.seed(1)
tree <- sim_tree_shape(100)
mc_a <- sim_metacommunity_A(tree, lambda = 1, sigma_signal = FALSE, h_signal = FALSE,
                            env_original = TRUE)

dyn <- time_it("layer B dynamics (pool 180)", \() {
  mc <- sim_metacommunity_B(sim_tree_shape(180), lambda = 1, sigma_signal = FALSE,
                            h_signal = FALSE, env_original = TRUE, dispersal = 0.01,
                            competition = "phylo")
  ncol(mc$comm)
})
dp <- time_it("test_dp_all extra (999)", \() {
  test_dp_all(mc_a$comm, mc_a$tree, mc_a$env, nperm = 999, extra_axes = TRUE)
})
pcc <- time_it("test_pcc (999)", \() {
  test_pcc(mc_a$comm, mc_a$tree, mc_a$env, nperm = 999)
})

message("check: species left ", paste(unlist(dyn), collapse = " "))
print(as.data.frame(dp[[1]][, c("test", "stat", "p_site", "p_taxa")]), digits = 10)
print(as.data.frame(pcc[[1]][, c("test", "stat", "p_site", "p_taxa")]), digits = 10)
