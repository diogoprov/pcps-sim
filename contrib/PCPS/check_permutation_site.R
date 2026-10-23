# Checks of the PCPS contribution (permut.torus() and the permutation.site
# argument), run against a build of PCPS with the patch applied. PCPS has no
# test suite, so this script is not part of the pull request.
# Usage: Rscript check_permutation_site.R   (with the patched PCPS installed)
library(PCPS)
stopifnot(exists("permut.torus", mode = "function"))

# permut.torus(): counts, permutations, toroidal distances preserved
side <- 10
perm <- permut.torus(side)
stopifnot(nrow(perm) == 8 * side^2 - 1, ncol(perm) == side^2,
          all(apply(perm, 1, function(p) all(sort(p) == seq_len(side^2)))),
          !any(duplicated(perm)))
stopifnot(nrow(permut.torus(side, symmetries = FALSE)) == side^2 - 1,
          nrow(permut.torus(4, 6)) == 4 * 4 * 6 - 1,
          nrow(permut.torus(side, runs = 99)) == 99)
rr <- rep(seq_len(side) - 1, each = side)
cc <- rep(seq_len(side) - 1, times = side)
torus_dist <- function(i, j) {
  dr <- abs(rr[i] - rr[j]); dc <- abs(cc[i] - cc[j])
  sqrt(pmin(dr, side - dr)^2 + pmin(dc, side - dc)^2)
}
d0 <- outer(seq_len(side^2), seq_len(side^2), torus_dist)
for (b in sample(nrow(perm), 20)) {
  p <- perm[b, ]
  stopifnot(isTRUE(all.equal(d0[p, p], d0)))
}

# a small metacommunity on a 10 x 10 grid
set.seed(1)
n_sp <- 15
tree <- ape::rcoal(n_sp)
phylodist <- ape::cophenetic.phylo(tree)
envir <- data.frame(E = as.vector(outer(1:side, 1:side, "+")))  # spatially structured
optima <- stats::runif(n_sp, 2, 20)
comm <- sapply(optima, function(u) stats::rpois(side^2, 5 * exp(-(envir$E - u)^2 / 20)))
colnames(comm) <- rownames(phylodist) <- colnames(phylodist) <- tree$tip.label
rownames(comm) <- rownames(envir) <- paste0("s", seq_len(side^2))
keep <- rowSums(comm) > 0
stopifnot(all(keep))  # permut.torus() needs the complete grid

# without permutation.site the results are unchanged
set.seed(2)
a <- pcps.sig(comm, phylodist, envir, FUN = FUN.GLM, formula = pcps.1 ~ E, choices = 1, runs = 49)
set.seed(2)
b <- pcps.sig(comm, phylodist, envir, FUN = FUN.GLM, formula = pcps.1 ~ E, choices = 1, runs = 49,
              permutation.site = NULL)
stopifnot(identical(a$p.site.shuffle, b$p.site.shuffle), identical(a$p.taxa.shuffle, b$p.taxa.shuffle))

# with torus translations: one null statistic per permutation, P value over them
tor <- permut.torus(side)
d <- pcps.sig(comm, phylodist, envir, FUN = FUN.GLM, formula = pcps.1 ~ E, choices = 1, runs = 49,
              permutation.site = tor)
stopifnot(nrow(d$statistic.null.site) == nrow(tor), nrow(d$statistic.null.taxa) == 49)
stopifnot(isTRUE(all.equal(d$p.site.shuffle,
                           (sum(d$statistic.null.site[, 1] >= d$obs.statistic[1]) + 1) / (nrow(tor) + 1))))
# (FUN.MANTEL: FUN.ADONIS calls vegan::adonis, removed from current vegan)
m <- matrix.p.sig(comm, phylodist, envir, FUN = FUN.MANTEL, method.p = "bray",
                  method.envir = "euclidean", runs = 49, permutation.site = tor)
stopifnot(nrow(m$statistic.null.site) == nrow(tor))

# invalid permutations are refused
bad <- tor[1:5, ]
bad[1, 1:2] <- 1
stopifnot(inherits(try(pcps.sig(comm, phylodist, envir, FUN = FUN.GLM, formula = pcps.1 ~ E,
                                choices = 1, runs = 9, permutation.site = bad), silent = TRUE),
                   "try-error"))
cat("All checks passed\n")
cat("Site-shuffle P: unrestricted", a$p.site.shuffle, "| torus", d$p.site.shuffle, "\n")
