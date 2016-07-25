#' Estimate number of clusters from global cluster assignments.
#' @param assignments Matrix of cluster assignments, where each row
#'     corresponds to cluster assignments sampled in one MCMC iteration
#'
#' @return Number of unique clusters in each MCMC iteration.
#'
#' @examples
#' # Generate simple test dataset
#' groupCounts <- c(50, 10, 40, 60)
#' means <- c(-1.5,1.5)
#' testData <- generateTestData_2D(groupCounts, means)
#' datasets <- testData$data
#'
#' # Fit the model
#' # 1. specify number of clusters
#' clusterCounts <- list(global=10, context=c(3,3))
#' # 2. Run inference
#' # Number of iterations is just for demonstration purposes, use
#' # a larger number of iterations in practice!
#' results <- contextCluster(datasets, clusterCounts,
#'      maxIter = 10, burnin = 5, lag = 1,
#'      dataDistributions = 'diagNormal',
#'      verbose = TRUE)
#'
#' # Extract only the sampled global assignments
#' samples <- results$samples
#' clusters <- plyr::laply(1:length(samples), function(i) samples[[i]]$Global)
#' numberOfClusters(clusters)
#'
#' @importFrom plyr laply
#' @importFrom magrittr %>%
#' @export
numberOfClusters <- function(assignments) {
  laply(1:nrow(assignments), function(i) {
    unique(assignments[i,]) %>% length
  })
}

#' Estimate sizes of clusters from global cluster assignments.
#' @param assignments Matrix of cluster assignments, where each row
#'     corresponds to cluster assignments sampled in one MCMC iteration
#'
#' @return Sizes of individual clusters in each MCMC iteration.
#'
#' @examples
#' # Generate simple test dataset
#' groupCounts <- c(50, 10, 40, 60)
#' means <- c(-1.5,1.5)
#' testData <- generateTestData_2D(groupCounts, means)
#' datasets <- testData$data
#'
#' # Fit the model
#' # 1. specify number of clusters
#' clusterCounts <- list(global=10, context=c(3,3))
#' # 2. Run inference
#' # Number of iterations is just for demonstration purposes, use
#' # a larger number of iterations in practice!
#' results <- contextCluster(datasets, clusterCounts,
#'      maxIter = 10, burnin = 5, lag = 1,
#'      dataDistributions = 'diagNormal',
#'      verbose = TRUE)
#'
#' # Extract only the sampled global assignments
#' samples <- results$samples
#' clusters <- plyr::laply(1:length(samples), function(i) samples[[i]]$Global)
#' clusterSizes(clusters)
#'
#' @importFrom magrittr %>%
#' @export
clusterSizes <- function(assignments) {
  clusterLabels <- unique(assignments %>% unlist)
  sizes <- matrix(nrow=nrow(assignments), ncol=length(clusterLabels))
  for (ci in 1:length(clusterLabels)) {
    sizes[,ci] <- rowSums(assignments == clusterLabels[ci])
  }

  sizes <- sizes %>% as.data.frame
  colnames(sizes) <- clusterLabels
  sizes
}

#' Compute the posterior co-clustering matrix from global cluster assignments.
#' @param assignments Matrix of cluster assignments, where each row
#'     corresponds to cluster assignments sampled in one MCMC iteration
#'
#' @return Posterior co-clustering matrix, where element \code{[i,j]} represents
#'     the posterior probability that data points \code{i} and \code{j} belong
#'     to the same cluster.
#'
#' @examples
#' # Generate simple test dataset
#' groupCounts <- c(50, 10, 40, 60)
#' means <- c(-1.5,1.5)
#' testData <- generateTestData_2D(groupCounts, means)
#' datasets <- testData$data
#'
#' # Fit the model
#' # 1. specify number of clusters
#' clusterCounts <- list(global=10, context=c(3,3))
#' # 2. Run inference
#' # Number of iterations is just for demonstration purposes, use
#' # a larger number of iterations in practice!
#' results <- contextCluster(datasets, clusterCounts,
#'      maxIter = 10, burnin = 5, lag = 1,
#'      dataDistributions = 'diagNormal',
#'      verbose = TRUE)
#'
#' # Extract only the sampled global assignments
#' samples <- results$samples
#' clusters <- plyr::laply(1:length(samples), function(i) samples[[i]]$Global)
#' coclusteringMatrix(clusters)
#'
#' @export
coclusteringMatrix <- function(assignments) {
  N <- ncol(assignments)
  coclust <- matrix(nrow=N, ncol=N)
  coclust[,] <- 0
  for (i in 1:N)  {
    for (j in (i+1):N) {
      if ((i < j) && (j <= N)) {  # indexing in R is awful
        coclust[i,j] <- sum(assignments[,i] == assignments[,j])
        coclust[j,i] <- coclust[i,j]
      }
    }
  }
  coclust <- coclust / nrow(assignments)
  coclust
}

# Get assignments to global and local clusters from raw MCMC samples
getClustersAssignments <- function(state) {
  nDataPoints <- length(state$dataAssignments)
  nContexts <- ncol(state$contextAssignments)
  clusterAssgn <- matrix(nrow=nDataPoints,ncol=nContexts)
  for (context in 1:nContexts) {
    clusterAssgn[,context] <- state$contextAssignments[state$dataAssignments, context]
  }

  # Map to global clusters
  globalAssgn <-
    aaply(clusterAssgn,1,function(row)
      as.character(row) %>% paste(collapse='-'))

  clusterAssgn <- as.data.frame(clusterAssgn)
  colnames(clusterAssgn) <- llply(1:nContexts, function (c) paste("Context", c))
  clusterAssgn$Global <- globalAssgn
  return(clusterAssgn)
}

# Get MCMC samples after burnin and with a specified lag
getMCMCSamples <- function(allSamples, burnin, lag = 5) {
  n <- length(allSamples)
  if (n > burnin) {
    idxs <- seq(burnin, n, by=lag)
    allSamples[idxs]
  } else {
    stop("Cannot discard burnin MCMC iterations because of insufficient number of samples.")
  }
}

