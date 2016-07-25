#' Fit an empirical Bayes prior to the data
#'
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
#' @param type Type of prior that is fitted to the data. The algorithm can fit either
#'   rate of the prior covariance matrix, or fit the full covariance matrix to the data.
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
#' prior <- empiricalBayesPrior(datasets, fullDataDistributions, 0.01, 0.1, 'fitRate')
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
#' @importFrom magrittr %>%
#' @importFrom MASS fitdistr
#' @importFrom plyr aaply
#' @importFrom stats var
#' @export

empiricalBayesPrior <- function(datasets, distributions='diagNormal',
                          globalConcentration=0.1, localConcentration=0.1, type='fitRate') {
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

  justFitGammaRate <- function(datasetIdx) {
    rateValues <- aaply(datasets[[datasetIdx]], 2, var) %>% as.numeric
    list(gammaShape=1, gammaRate=rateValues, beta = 1)
  }

  fitEverything <- function(datasetIdx) {
    nDim <- length(datasets[[datasetIdx]][1,])
    precisionValues <- aaply(datasets[[datasetIdx]], 2, function(observations) 1/(var(observations)))
    # precision parameters
    gammaParams <-
      tryCatch({
        fitdistr(precisionValues, 'gamma')$estimate %>% as.numeric
      }, error=function(e) { print(paste(datasetIdx, ' failed:', e))
                             c(1, 1/(mean(precisionValues)))}
      )

    # precision of the mean - we expect this to be relatively large - equal to the precision of the data?
    beta <- 1

    list(gammaShape=gammaParams[1], gammaRate=rep(gammaParams[2], nDim), beta=beta)
  }

  # Estimate values for diagonal normal distribution
  estimatePrior_diagNormal <- function(datasetIdx) {
    switch(type,
           fitRate=justFitGammaRate(datasetIdx),
           fitFullVariance=fitEverything(datasetIdx))
  }

  getDiagNormalPrior <- function(context) {
    nDim <- length(datasets[[context]][1,])
    m = rep(0.0, nDim)

    fitParams <- estimatePrior_diagNormal(context)

    beta <- fitParams$beta
    a = fitParams$gammaShape
    b = fitParams$gammaRate
    list(beta=beta, m=m, a = a, b = b)
  }

  getFullNormalPrior <- function(context) {
    nDim <- length(datasets[[context]][1,])
    m <- rep(0.0, nDim)
    beta <- 1.0
    W = diag(1/nDim, nDim, nDim)
    v = max(nDim, 1.001)  # degrees of freedom > number of dimensions - 1
    list(beta=beta, m=m, Sigma=W, v=v)
  }

  priors <-
    lapply(1:nContexts, function(context){
      nDim <- length(datasets[[context]][1,])

      theta <-
        switch(dataDistributions[[context]],
               diagNormal=getDiagNormalPrior(context),
               fullNormal=getFullNormalPrior(context))
      alpha <- alphas[[context]]
      distribution <-
        switch(dataDistributions[[context]],
               fullNormal="fullNormal",
               diagNormal="diagNormal")
      list(theta=theta, alpha = alpha, distribution=distribution)
    })
  priors$gamma <- globalConcentration
  priors

}
