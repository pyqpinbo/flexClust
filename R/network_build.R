#' Build a Network from Distance Matrix
#'
#' @description
#' Constructs an undirected network graph based on a distance matrix, using a specified distance percentile threshold.
#'
#' @param distance_matrix A square distance matrix (n x n), where n represents the number of nodes.
#' @param percentile_value A numeric value between 0 and 1 representing the percentile threshold for edge formation.
#' @param normalize A logical value indicating whether to normalize the distance matrix before processing. Default is `TRUE`.
#'
#' @details
#' - The function normalizes the distance matrix (if `normalize = TRUE`) and identifies edges based on a percentile threshold.
#' - Nodes with distances below the given threshold are connected by edges in the resulting graph.
#'
#' @return
#' An undirected graph object of class `igraph`.
#'
#' @import igraph
#' @export
#'
#' @examples
#' set.seed(123)
#' distance_matrix <- matrix(runif(100), nrow = 10)
#' network <- network_build(distance_matrix, percentile_value = 0.5)
#' plot(network)
#'
network_build <- function(distance_matrix, percentile_value = 0.3, normalize = TRUE) {
  if (!is.matrix(distance_matrix) || nrow(distance_matrix) != ncol(distance_matrix)) {
    stop("'distance_matrix' must be a square matrix.")
  }
  
  if (!is.numeric(percentile_value) || percentile_value < 0 || percentile_value > 1) {
    stop("'percentile_value' must be a numeric value between 0 and 1.")
  }
  
  if (normalize) {
    distance_matrix <- dist.normalize(distance_matrix)
  }
  
  distance_threshold <- dist_percentile(distance_matrix, percentile_value)
  
  adjacency_matrix <- matrix(0, nrow = nrow(distance_matrix), ncol = ncol(distance_matrix))
  adjacency_matrix[distance_matrix < distance_threshold] <- 1
  
  graph <- graph.adjacency(adjacency_matrix, mode = "undirected", diag = FALSE)
  
  return(graph)
}

#' Compute Distance Percentile Threshold
#'
#' @description
#' Computes a distance threshold based on a given percentile from a square distance matrix.
#'
#' @param distance_matrix A square distance matrix.
#' @param percentile_value A numeric value between 0 and 1 representing the desired percentile threshold.
#'
#' @details
#' - The function extracts the upper triangle of the distance matrix (excluding the diagonal).
#' - Computes the percentile threshold based on the specified percentile value.
#'
#' @return
#' A numeric value representing the distance threshold at the given percentile.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' distance_matrix <- matrix(runif(100), nrow = 10)
#' threshold <- dist_percentile(distance_matrix, percentile_value = 0.5)
#' print(threshold)
#'
dist_percentile <- function(distance_matrix, percentile_value) {
  if (!is.matrix(distance_matrix) || nrow(distance_matrix) != ncol(distance_matrix)) {
    stop("'distance_matrix' must be a square matrix.")
  }
  
  if (!is.numeric(percentile_value) || percentile_value < 0 || percentile_value > 1) {
    stop("'percentile_value' must be a numeric value between 0 and 1.")
  }
  
  distance_matrix[is.na(distance_matrix)] <- +Inf
  upper_triangle <- distance_matrix[upper.tri(distance_matrix)]
  
  threshold <- quantile(upper_triangle, probs = percentile_value)
  return(threshold)
}

#' Normalize Distance Matrix
#'
#' @description
#' Normalizes a square distance matrix to a range of 0 to 1.
#'
#' @param dist A square distance matrix.
#'
#' @details
#' - Extracts the upper triangle of the distance matrix.
#' - Rescales the upper triangle values to fall between 0 and 1 using the `scales::rescale` function.
#' - Reflects the normalized values back into a symmetric distance matrix.
#'
#' @return
#' A normalized square distance matrix.
#'
#' @import scales
#' @export
#'
#' @examples
#' set.seed(123)
#' distance_matrix <- matrix(runif(100), nrow = 10)
#' normalized_matrix <- dist.normalize(distance_matrix)
#' print(normalized_matrix)
#'
dist.normalize <- function(dist) {
  dist_normalized <- matrix(0, nrow = nrow(dist), ncol = ncol(dist))
  
  upper_values <- dist[upper.tri(dist)]
  normalized_values <- scales::rescale(upper_values)
  
  dist_normalized[upper.tri(dist_normalized)] <- normalized_values
  dist_normalized <- dist_normalized + t(dist_normalized)
  
  return(dist_normalized)
}


#' Perform Network Clustering
#'
#' @description Clusters a network into communities.
#'
#' @param net An igraph object representing the network to be clustered.
#'
#' @return An object of class `communities`, representing the detected community structure.
#'
#' @export
#' 
network_clust <- function(net){
  net_layout = layout_components(net)
  communities = cluster_louvain(net)
  return(list(community = communities, layout = net_layout))
}


