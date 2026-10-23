# TODO before the PR: add a verified reference for the torus-translation test
# ([CITAR: torus-translation test]); none was checked yet.
#' @title Restricted site permutations for sampling units on a regular grid
#'
#' @description Generates the permutations of sampling units arranged on a
#' regular grid that preserve their spatial autocorrelation: all translations of
#' the grid wrapped as a torus, optionally combined with its reflections and
#' rotations. The result can be used as the argument \emph{permutation.site} of
#' \code{\link{pcps.sig}} and \code{\link{matrix.p.sig}}, which then replaces the
#' unrestricted site shuffle.
#'
#' @details When the sampling units are spatially autocorrelated, as along a
#' spatially structured environmental gradient, the unrestricted site shuffle
#' destroys this autocorrelation and the null distribution of the statistic
#' becomes too narrow, so the test of the environmental effect becomes liberal.
#' Torus translations move the whole pattern of the
#' environmental variables over the sampling units while keeping the distance
#' between any two units, and thus their spatial autocorrelation. The method
#' assumes that the sampling units form a complete regular grid with
#' \emph{nrow} rows and \emph{ncol} columns, and it treats opposite edges as
#' adjacent.
#'
#' The sampling units must be numbered row by row: unit \eqn{(r - 1) \times
#' ncol + c}{(r - 1) * ncol + c} is the one in row \eqn{r} and column
#' \eqn{c}, and the rows of the community data must follow this order. Each row
#' of the result is a vector \emph{p} such that unit \emph{i} receives the data of
#' unit \emph{p[i]} in that permutation. The identity is excluded. With
#' \emph{symmetries = TRUE}, the translations are combined with the reflections
#' of the grid and, for square grids, with its rotations and transpositions, which
#' gives up to \eqn{8 \times nrow \times ncol - 1}{8 * nrow * ncol - 1}
#' distinct permutations (\eqn{4 \times nrow \times ncol - 1}{4 * nrow * ncol - 1}
#' for rectangular grids). Since this number is limited by the size of the grid,
#' the smallest attainable P value is \eqn{1/(k + 1)}, where \emph{k} is the
#' number of permutations.
#'
#' @encoding UTF-8
#' @param nrow Number of rows of the grid.
#' @param ncol Number of columns of the grid (Default ncol = nrow).
#' @param symmetries Logical argument (TRUE or FALSE) to specify if the
#' reflections and rotations of the grid are combined with the translations
#' (Default symmetries = TRUE).
#' @param runs Number of permutations to return, sampled without replacement
#' among all distinct permutations. If NULL, all distinct permutations are
#' returned (Default runs = NULL).
#' @return A matrix with one permutation of the sampling units per row.
#' @seealso \code{\link{pcps.sig}}, \code{\link{matrix.p.sig}}
#' @keywords PCPS
#' @examples
#' # all 799 permutations of a 10 x 10 grid
#' perm <- permut.torus(10)
#' dim(perm)
#' # only the 99 translations
#' dim(permut.torus(10, symmetries = FALSE))
#' @export
permut.torus <- function(nrow, ncol = nrow, symmetries = TRUE, runs = NULL){
  if(length(nrow) != 1 || length(ncol) != 1 || nrow < 2 || ncol < 2 || nrow != round(nrow) || ncol != round(ncol)){
    stop("\n nrow and ncol must be single integers of at least 2 \n")
  }
  rr <- rep(seq_len(nrow) - 1, each = ncol)
  cc <- rep(seq_len(ncol) - 1, times = nrow)
  maps <- list(function(r, c) cbind(r, c))
  if(symmetries){
    maps <- c(maps, list(function(r, c) cbind(-r, c), function(r, c) cbind(r, -c),
                         function(r, c) cbind(-r, -c)))
    if(nrow == ncol){
      maps <- c(maps, list(function(r, c) cbind(c, r), function(r, c) cbind(-c, r),
                           function(r, c) cbind(c, -r), function(r, c) cbind(-c, -r)))
    }
  }
  perms <- list()
  for(f in maps){
    for(sr in seq_len(nrow) - 1){
      for(sc in seq_len(ncol) - 1){
        p <- f(rr, cc)
        perms[[length(perms) + 1]] <- ((p[, 1] + sr) %% nrow) * ncol + (p[, 2] + sc) %% ncol + 1
      }
    }
  }
  perms <- unique(do.call(rbind, perms))
  perms <- perms[!apply(perms, 1, function(p) all(p == seq_len(nrow * ncol))), , drop = FALSE]
  if(!is.null(runs)){
    if(runs > nrow(perms)){
      stop(paste("\n runs must not be larger than the number of distinct permutations (", nrow(perms), ") \n", sep = ""))
    }
    perms <- perms[sample.int(nrow(perms), runs), , drop = FALSE]
  }
  rownames(perms) <- NULL
  return(perms)
}
