#' Generate a basic prior distribution for the datasets.
#'
#' Creates a basic prior distribution for the clustering model, assuming a unit
#' prior covariance matrix for clusters in each dataset.
#' @param datasets  List of data matrices where each matrix represents a
#'   context-specific dataset. Each data matrix has the size \emph{N} times \emph{M}, where
#'   \emph{N} is the number of data points and \emph{M} is the dimensionality of the data.
#'   The full list of matrices has length \emph{C}. The number of data points \emph{N} must
#'   be the same for all data matrices.
#' @param distributions Distribution of data in each dataset. Can be either a list of
#'   length \emph{C} where \code{dataDistributions[c]} is the distribution of dataset \emph{c},
#'   or a single string when all datasets have the same distribution. Currently implemented
#'   distribution is the \code{'diagNormal'} option for multivariate Normal distribution
#'   with diagonal covariance  matrix.
#' @param globalConcentration Prior concentration parameter for the global clusters. Small
#'   values of this parameter give larger prior probability to smaller number of clusters.
#' @param localConcentration Prior concentration parameter for the local context-specific
#'   clusters. Small values of this parameter give larger prior probability to smaller number of clusters.
#'
#' @return Returns the prior object that can be used as an input for the \code{contextCluster}
#'   function.
#'
#' @examples
#' # Example with simulated data (see vignette for details)
#' nContexts <- 2
#' # Number of elements in each cluster
#' groupCounts <- c(50, 10, 40, 60)
#' # Centers of clusters
#' means <- c(-1.5,1.5)
#' testData <- generateTestData_2D(groupCounts, means)
#' datasets <- testData$data
#'
#' # Generate the prior
#' fullDataDistributions <- rep('diagNormal', nContexts)
#' prior <- generatePrior(datasets, fullDataDistributions, 0.01, 0.1)
#'
#' # Fit the model
#' # 1. specify number of clusters
#' clusterCounts <- list(global=10, context=c(3,3))
#' # 2. Run inference
#' # Number of iterations is just for demonstration purposes, use
#' # a larger number of iterations in practice!
#' results <- contextCluster(datasets, clusterCounts,
#'      maxIter = 10, burnin = 5, lag = 1,
#'      dataDistributions = 'diagNormal', prior = prior,
#'      verbose = TRUE)
#'
#'
#' @export
generatePrior <- function(datasets, distributions='diagNormal',
                          globalConcentration=0.1, localConcentration=0.1) {
  nContexts <- length(datasets)

  dataDistributions <- vector(length=nContexts, mode="character")
  if (length(distributions) == 1) {
    dataDistributions[] <- distributions
  } else {
    if (length(distributions) == nContexts) dataDistributions <- distributions
    else stop(paste(
      "Number of distributions does not match number of datasets.\n",
      "Expected: ", nContexts, " or 1 (common for all datasets)\n",
      "Given: ",length(distributions)))
  }
  # concentration parameter value for each dataset
  alphas <- vector(length=nContexts, mode='numeric')
  if (length(localConcentration) == 1) {
    alphas[] <- localConcentration
  } else {
    if (length(localConcentration) == nContexts) alphas <- localConcentration
    else stop(paste(
      "Number of local concentration parameters does not match number of datasets.\n",
      "Expected: ", nContexts, " or 1 (common for all datasets)\n",
      "Given: ",length(localConcentration)))
  }

  getFullNormalPrior <- function(nDim) {
    m = rep(0,nDim)   # zero mean
    beta = 1.0
    W = diag(1/nDim, nDim, nDim)
    v = max(nDim, 1.001)  # degrees of freedom > number of dimensions - 1
    list(beta=beta, m=m, Sigma=W, df=v)
  }

  getDiagNormalPrior <- function(nDim) {
    beta <- 1.0
    m = rep(0.0, nDim)
    a = 1.0
    b = rep(1.0, nDim)
    list(beta=beta, m=m, a = a, b = b)
  }

  priors <-
    lapply(1:nContexts, function(context){
      nDim <- length(datasets[[context]][1,])

      theta <-
        switch(dataDistributions[[context]],
             normal=getFullNormalPrior(nDim),
             fullNormal=getFullNormalPrior(nDim),
             diagNormal=getDiagNormalPrior(nDim),
             normalDiagCovariance=getDiagNormalPrior(nDim))
      alpha <- alphas[[context]]
      distribution <-
        switch(dataDistributions[[context]],
               normal="fullNormal",
               fullNormal="fullNormal",
               diagNormal="diagNormal",
               normalDiagCovariance="diagNormal")
      list(theta=theta, alpha = alpha, distribution=distribution)
    })
  priors$gamma <- globalConcentration
  priors
}
