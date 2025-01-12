#' Calculate the Distance Matrix of Multiple Functional or Time-Series Curves
#'
#' Computes pairwise distances between multiple time-series curves using a specified distance function, with optional parallel computing support.
#'
#' @param tsList A list of multiple time-series curves.
#' @param distFunc A distance function, such as \code{tsdiss.euclidean}, \code{tsdiss.manhattan}, or \code{tsdiss.correlation}.
#' @param cores The number of CPU cores to use for parallel computing. Defaults to the maximum available cores minus one.
#' @keywords distance matrix
#' @import TSclust
#' @import TSdist
#' @import compiler
#' @import parallel
#' @import igraph
#' @export
#'
#' @examples
#' tsList <- list(rnorm(100), rnorm(100))
#' distMat <- dist.parallel(tsList, distFunc = tsdiss.euclidean, cores = 2)
#' print(distMat)
#'
dist.parallel <- function(tsList, distFunc = tsdiss.euclidean, cores = NUM_CORES) {
  distFuncCompiled <- cmpfun(distFunc)
  tsListLength <- length(tsList)
  combs <- combn(tsListLength, 2)
  d <- mcmapply(
    dist.parallel2.compute, 
    x = combs[1,], 
    y = combs[2,], 
    MoreArgs = list(tsList = tsList, distFunc = distFuncCompiled), 
    mc.cores = cores
  )
  dist <- matrix(0, tsListLength, tsListLength)
  dist[lower.tri(dist)] <- d
  dist <- as.matrix(as.dist(dist))
}

#' Compute Pairwise Distance Between Two Time-Series Curves
#'
#' Computes the distance between two time-series curves using a specified distance function.
#'
#' @param x Index of the first time-series curve in the list.
#' @param y Index of the second time-series curve in the list.
#' @param tsList A list of time-series curves.
#' @param distFunc A compiled distance function.
#' @return A numeric distance value between two time-series curves.
#' @export
#'
dist.parallel2.compute <- function(x, y, tsList, distFunc) { 
  distFunc(tsList[[x]], tsList[[y]])
}

#' Correlation-Based Distance Between Time-Series Curves
#'
#' Computes the distance between two time-series curves based on their correlation.
#'
#' @param ts1 First time-series curve.
#' @param ts2 Second time-series curve.
#' @return A numeric distance value based on correlation.
#' @export
#'
tsdiss.correlation <- function(ts1, ts2) {
  CorDistance(ts1, ts2)
} 

#' Euclidean Distance Between Time-Series Curves
#'
#' Computes the Euclidean distance between two time-series curves.
#'
#' @param ts1 First time-series curve.
#' @param ts2 Second time-series curve.
#' @return A numeric Euclidean distance value.
#' @export
#'
tsdiss.euclidean <- function(ts1, ts2) {
  diss.EUCL(ts1, ts2)
} 

#' Manhattan Distance Between Time-Series Curves
#'
#' Computes the Manhattan (L1 norm) distance between two time-series curves.
#'
#' @param ts1 First time-series curve.
#' @param ts2 Second time-series curve.
#' @return A numeric Manhattan distance value.
#' @export
#'
tsdiss.manhattan <- function(ts1, ts2) {
  sum(abs(ts1 - ts2))
} 

#' Infinite Norm Distance Between Time-Series Curves
#'
#' Computes the infinite norm (maximum absolute difference) distance between two time-series curves.
#'
#' @param ts1 First time-series curve.
#' @param ts2 Second time-series curve.
#' @return A numeric infinite norm distance value.
#' @export
#'
tsdiss.infiniteNorm <- function(ts1, ts2) {
  max(abs(ts1 - ts2))
}






