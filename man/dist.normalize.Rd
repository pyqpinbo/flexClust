% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_build.R
\name{dist.normalize}
\alias{dist.normalize}
\title{Normalize Distance Matrix}
\usage{
dist.normalize(dist)
}
\arguments{
\item{dist}{A square distance matrix.}
}
\value{
A normalized square distance matrix.
}
\description{
Normalizes a square distance matrix to a range of 0 to 1.
}
\details{
\itemize{
\item Extracts the upper triangle of the distance matrix.
\item Rescales the upper triangle values to fall between 0 and 1 using the \code{scales::rescale} function.
\item Reflects the normalized values back into a symmetric distance matrix.
}
}
\examples{
set.seed(123)
distance_matrix <- matrix(runif(100), nrow = 10)
normalized_matrix <- dist.normalize(distance_matrix)
print(normalized_matrix)

}
