#' Normalized Mutual Information (NMI) Between Two Partitions
#'
#' @description
#' Computes the Normalized Mutual Information (NMI) to measure the similarity between two partitions of a set.
#'
#' @param C A list representing one partition of a set S. Each element in the list is a subset of indices from the set.
#' @param D A list representing another partition of the same set S.
#' @param n An integer representing the cardinality (size) of the set S.
#'
#' @details
#' - **Normalized Mutual Information (NMI):** Measures the shared information between two partitions normalized by their entropies.
#' - NMI ranges from 0 to 1, where 1 indicates identical partitions and 0 indicates completely independent partitions.
#' - The calculation is based on mutual information (`Icd`) and entropy (`Hc`) of each partition.
#'
#' @return A numeric value representing the NMI score between the two partitions.
#'
#' @keywords normalized mutual information, clustering evaluation, partition similarity
#' @export
#'
#' @examples
#' n <- 20
#' C <- list(1:(n/2+1), (n/2+2):n, (n+1):(2*n))
#' D <- list(1:(n/2), (n/2+1):n, (n+1):(2*n))
#' res <- NMI(C, D, 2 * n)
#' print(res)

NMI <- function(C, D, n) {
  if (!is.list(C) || !is.list(D)) stop("'C' and 'D' must be lists representing partitions.")
  if (!is.numeric(n) || n <= 0) stop("'n' must be a positive numeric value representing set cardinality.")
  
  # Compute NMI based on mutual information and entropy
  mutual_info <- Icd(C, D, n)
  entropy_C <- Hc(C, n)
  entropy_D <- Hc(D, n)
  
  if ((entropy_C + entropy_D) == 0) {
    warning("The sum of entropies is zero, NMI is undefined. Returning NA.")
    return(NA)
  }
  
  res <- 2 * mutual_info / (entropy_C + entropy_D)
  return(res)
}

#' Mutual Information Between Two Partitions
#'
#' @description
#' Computes the mutual information shared between two partitions.
#'
#' @param C A list representing one partition.
#' @param D A list representing another partition.
#' @param n An integer representing the set size.
#'
#' @return A numeric value representing the mutual information between partitions.
Icd <- function(C, D, n) {
  K1 <- length(C)
  K2 <- length(D)
  mutual_info <- 0
  
  for (k in 1:K1) {
    Ck <- C[[k]]
    c1 <- length(Ck)
    if (c1 == 0) next
    
    for (j in 1:K2) {
      Dj <- D[[j]]
      d1 <- length(Dj)
      if (d1 == 0) next
      
      tmp <- length(intersect(Ck, Dj))
      if (tmp > 0) {
        mutual_info <- mutual_info + (tmp / n) * log((n * tmp) / (c1 * d1))
      }
    }
  }
  
  return(mutual_info)
}

#' Entropy of a Partition
#'
#' @description
#' Computes the entropy of a partition, representing the uncertainty within the partition.
#'
#' @param C A list representing a partition.
#' @param n An integer representing the set size.
#'
#' @return A numeric value representing the entropy of the partition.
Hc <- function(C, n) {
  entropy <- 0
  
  for (k in seq_along(C)) {
    Ck <- C[[k]]
    c1 <- length(Ck)
    if (c1 > 0) {
      entropy <- entropy - (c1 / n) * log(c1 / n)
    }
  }
  
  return(entropy)
}
