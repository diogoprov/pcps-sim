# Diagnosis of the inflated taxa-shuffle type I error of the GLM on PCPS2 in
# scenario S2sigma (phylogenetic signal in niche breadth only; pilot_glm).
# For S2 (control) and S2sigma, two grid cells and N = 100/200 (incidence),
# records per replicate:
#   - taxa-shuffle p of the GLM on PCPS1 and PCPS2 as in the pipeline (one axis
#     at a time, no Procrustes), with symmetric Procrustes of null axes 1-2 onto
#     the observed axes (pcps.sig with choices = 1:2), for an RDA on axes 1-2,
#     and for a quadratic GLM on PCPS2 (env + env^2);
#   - correlations of PCPS1/PCPS2 with env, env^2 and the community-weighted
#     mean of sigma (CWM sigma), and of CWM sigma with env and env^2;
#   - share of axis 2 in the eigenvalues, observed and mean of the taxa null.
# Usage: Rscript analysis/09_diag_s2sigma.R <n_reps> <workers> <nperm>

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[[1]]) else 100L
workers <- if (length(args) >= 2) as.integer(args[[2]]) else 6L
nperm <- if (length(args) >= 3) as.integer(args[[3]]) else 499L

for (f in list.files("R", full.names = TRUE)) source(f)
future::plan(future::multisession, workers = workers)
opts <- furrr::furrr_options(seed = 20261017L, globals = TRUE,
                             packages = c("dplyr", "tibble"))
out_dir <- file.path("data", "derived", "diag_s2sigma")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

grid <- readRDS(file.path("data", "derived", "trees", "tree_grid_k200.rds"))
cells <- c("g_yule__c_yule", "g_low__c_yule")
trees <- dplyr::filter(grid, cell %in% cells)
tree_of <- \(cell, n_spp, rep)
  trees$trees[[which(trees$cell == cell & trees$n == n_spp)]][[rep]]

# PCoA axes 1-2 of D_P and the share of each in the positive eigenvalues
axes_dp <- function(comm, q) {
  g <- gower_centre(d_p(matrix_p(comm, q), sqrt = TRUE))
  ev <- eigen(g, symmetric = TRUE)
  pos <- pmax(ev$values, 0)
  list(y = sweep(ev$vectors[, 1:2], 2, sqrt(pos[1:2]), "*"),
       share = pos[1:2] / sum(pos))
}
# F of a multivariate regression (RDA) of y on the columns of x
f_rda <- function(y, x) {
  x <- cbind(1, x)
  fit <- x %*% qr.solve(x, y)
  yc <- sweep(y, 2, colMeans(y))
  ss_reg <- sum(sweep(fit, 2, colMeans(y))^2)
  ss_res <- sum((y - fit)^2)
  df1 <- ncol(x) - 1
  (ss_reg / df1) / (ss_res / (nrow(y) - df1 - 1))
}
procrustes_fit <- function(obs, null)
  stats::fitted(vegan::procrustes(obs, null, symmetric = TRUE))

diag_rep <- function(scenario, cell, n_spp, rep) {
  sc <- dplyr::filter(scenarios_layerA(), scenario == !!scenario)
  tree <- tree_of(cell, n_spp, rep)
  mc <- sim_metacommunity_A(tree, sc$lambda, sc$sigma_signal, sc$h_signal,
                            sc$env_original, incidence = TRUE)
  comm <- mc$comm
  env <- mc$env
  q <- phylo_q(mc$tree)
  sig <- mc$traits[colnames(comm), "sigma"]
  cwm_sigma <- drop(comm %*% sig) / rowSums(comm)
  env2 <- (env - mean(env))^2

  obs <- axes_dp(comm, q)
  y <- obs$y
  stat_set <- function(y_ax, y_adj) c(
    f_axes(y_ax, env),                               # per axis, no Procrustes
    f_axes(y_adj, env),                              # per axis, Procrustes on 1-2
    f_rda(y_ax, env),                                # RDA on axes 1-2
    f_rda(y_ax[, 2, drop = FALSE], cbind(env, env2)) # quadratic GLM on PCPS2
  )
  s_obs <- stat_set(y, y)
  k <- ncol(comm)
  null <- replicate(nperm, {
    j <- sample.int(k)
    nl <- axes_dp(comm, q[j, j])
    c(stat_set(nl$y, procrustes_fit(y, nl$y)), nl$share[2])
  })
  p <- vapply(seq_along(s_obs), \(a) perm_p(s_obs[a], null[a, ]), 0)
  tibble::tibble(
    scenario = scenario, cell = cell, n_spp = n_spp, rep = rep,
    p_pcps1 = p[1], p_pcps2 = p[2], p_pcps1_proc = p[3], p_pcps2_proc = p[4],
    p_rda12 = p[5], p_pcps2_quad = p[6],
    r_pcps1_env = stats::cor(y[, 1], env), r_pcps2_env = stats::cor(y[, 2], env),
    r_pcps2_env2 = stats::cor(y[, 2], env2),
    r_pcps1_cwm = stats::cor(y[, 1], cwm_sigma),
    r_pcps2_cwm = stats::cor(y[, 2], cwm_sigma),
    r_cwm_env = stats::cor(cwm_sigma, env), r_cwm_env2 = stats::cor(cwm_sigma, env2),
    share2_obs = obs$share[2], share2_null = mean(null[7, ]),
    r_sigma_tree = {
      # phylogenetic signal of sigma among the species kept (Blomberg's K)
      phytools::phylosig(mc$tree, setNames(sig, colnames(comm)), method = "K")[1]
    }
  )
}

jobs <- tidyr::crossing(scenario = c("S2", "S2sigma"), cell = cells,
                        n_spp = c(100, 200), rep = seq_len(n_reps))
t0 <- Sys.time()
res <- furrr::future_pmap(jobs, diag_rep, .options = opts) |> dplyr::bind_rows()
message("Diagnosis: ", format(Sys.time() - t0))
readr::write_csv(res, file.path(out_dir, "diag_s2sigma.csv"))

res |>
  dplyr::group_by(scenario, cell, n_spp) |>
  dplyr::summarise(dplyr::across(dplyr::starts_with("p_"), \(p) mean(p <= 0.05)),
                   dplyr::across(c(dplyr::starts_with("r_"), dplyr::starts_with("share")),
                                 \(x) mean(abs(x))),
                   .groups = "drop") |>
  as.data.frame() |>
  print(digits = 2)
