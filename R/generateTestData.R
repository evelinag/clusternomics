#' Generate simulated 1D dataset for testing
#'
#' Generate simple 1D dataset with two contexts, where the data are generated
#' from Gaussian distributions.
#' The generated output contains two datasets, where
#' each dataset contains 4 global clusters, originating from
#' two local clusters in each context.
#' @param groupCounts Number of data samples in each global cluster.
#'     It is assumed to be a vector of four elements: \code{ c(c11, c21, c12, c22)}
#'     where \code{cij} is the number of samples coming from cluster i in context 1
#'     and cluster j in context 2.
#' @param means Means of the simulated clusters.
#'     It is assumed to be a vector of two elements: \code{ c(m1, m2)}
#'     where \code{m1} is the mean of the first cluster in both contexts,
#'     and \code{m2} is the mean of the second cluster in both contexts.
#'
#' @return Returns the simulated datasets together with true assignmets.
#'   \item{data}{List of datasets for each context. This can be used as an input
#'      for the \code{contextCluster} function.}
#'   \item{groups}{True cluster assignments that were used to generate the data.}
#'
#' @examples
#' groupCounts <- c(50, 10, 40, 60)
#' means <- c(-1.5,1.5)
#' testData <- generateTestData_1D(groupCounts, means)
#' # Use the dataset as an input for the contextCluster function for testing
#' datasets <- testData$data
#'
#' @importFrom stats rnorm
#' @export
generateTestData_1D <- function(groupCounts, means) {
  vars <- 1
  groupDistributions <-
    data.frame(
      mean1 = c(means[1], means[2], means[1], means[2]),
      mean2 = c(means[1], means[1], means[2], means[2]),
      count = groupCounts)

  # generate data
  x1raw <- matrix(nrow=sum(groupCounts), ncol=1)
  x2raw <- matrix(nrow=sum(groupCounts), ncol=1)
  assignments <- vector(length=sum(groupCounts), mode='integer')
  for (g in 1:nrow(groupDistributions)) {
    (beg <- max(groupDistributions$count[1:g-1] %>% sum,0) + 1)
    (fin <- groupDistributions$count[1:g] %>% sum)
    if (beg <= fin) {
      x1raw[beg:fin] <- rnorm(groupDistributions$count[g], groupDistributions$mean1[g], vars)
      x2raw[beg:fin] <- rnorm(groupDistributions$count[g], groupDistributions$mean2[g], vars)
      assignments[beg:fin] <- g
    }
  }

  # permute data
  permutation <- sample(sum(groupCounts))
  x1 <- x1raw[permutation,] %>% matrix
  x2 <- x2raw[permutation,] %>% matrix
  groups <- assignments[permutation]
  return(list(data=list(x1, x2), groups=groups))
}

#' Generate simulated 2D dataset for testing
#'
#' Generate simple 2D dataset with two contexts, where the data are generated
#' from Gaussian distributions.
#' The generated output contains two datasets, where
#' each dataset contains 4 global clusters, originating from
#' two local clusters in each context.
#' @param groupCounts Number of data samples in each global cluster.
#'     It is assumed to be a vector of four elements: \code{ c(c11, c21, c12, c22)}
#'     where \code{cij} is the number of samples coming from cluster i in context 1
#'     and cluster j in context 2.
#' @param means Means of the simulated clusters.
#'     It is assumed to be a vector of two elements: \code{ c(m1, m2)}
#'     where \code{m1} is the mean of the first cluster in both contexts,
#'     and \code{m2} is the mean of the second cluster in both contexts.
#'     Because the data are two-dimensional, the mean is assumed to be the
#'     same in both dimensions.
#' @param variances Optionally, it is possible to specify different variance
#'     for each of the clusters. The variance is assumed to be a vector
#'     of two elements: \code{ c(v1, v2)}
#'     where \code{v1} is the variance of the first cluster in both contexts,
#'     and \code{v2} is the variance of the second cluster in both contexts.
#'     Because the data are two-dimensional, the variance is diagonal and the
#'     same in both dimensions.
#'
#' @return Returns the simulated datasets together with true assignmets.
#'   \item{data}{List of datasets for each context. This can be used as an input
#'      for the \code{contextCluster} function.}
#'   \item{groups}{True cluster assignments that were used to generate the data.}
#'
#' @examples
#' groupCounts <- c(50, 10, 40, 60)
#' means <- c(-1.5,1.5)
#' testData <- generateTestData_1D(groupCounts, means)
#' # Use the dataset as an input for the contextCluster function for testing
#' datasets <- testData$data
#'
#' @importFrom stats rnorm
#' @export
generateTestData_2D <- function(groupCounts, means, variances=NULL) {
  if (is.null(variances)) {
    vars <- rep(1,length(groupCounts))
  } else {
    vars <- variances
  }
  groupDistributions <-
    data.frame(
      mean1 = c(means[1], means[2], means[1], means[2]),
      mean2 = c(means[1], means[1], means[2], means[2]),
      count = groupCounts)

  # generate data
  x1raw <- matrix(nrow=sum(groupCounts), ncol=2)
  x2raw <- matrix(nrow=sum(groupCounts), ncol=2)
  assignments <- vector(length=sum(groupCounts), mode='integer')
  for (g in 1:nrow(groupDistributions)) {
    (beg <- max(groupDistributions$count[1:g-1] %>% sum,0) + 1)
    (fin <- groupDistributions$count[1:g] %>% sum)
    if (beg <= fin ) {
      x1raw[beg:fin,1] <- rnorm(groupDistributions$count[g], groupDistributions$mean1[g], vars[g])
      x2raw[beg:fin,1] <- rnorm(groupDistributions$count[g], groupDistributions$mean2[g], vars[g])
      x1raw[beg:fin,2] <- rnorm(groupDistributions$count[g], groupDistributions$mean1[g], vars[g])
      x2raw[beg:fin,2] <- rnorm(groupDistributions$count[g], groupDistributions$mean2[g], vars[g])
      assignments[beg:fin] <- g
    }
  }

  # permute data
  permutation <- sample(sum(groupCounts))
  x1 <- x1raw[permutation,] %>% as.matrix
  x2 <- x2raw[permutation,] %>% as.matrix
  groups <- assignments[permutation]
  return(list(data=list(x1, x2), groups=groups))
}


