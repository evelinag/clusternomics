#' Clusternomics: Context-dependent clustering
#'
#' This function fits the context-dependent clustering model
#' to the data using Gibbs sampling. It allows the user to specify a different
#' number of clusters on the global level, as well as on the local level.
#' @param datasets List of data matrices where each matrix represents a
#'   context-specific dataset. Each data matrix has the size \emph{N} times \emph{M}, where
#'   \emph{N} is the number of data points and \emph{M} is the dimensionality of the data.
#'   The full list of matrices has length \emph{C}. The number of data points \emph{N} must
#'   be the same for all data matrices.
#' @param clusterCounts Number of cluster on the global level and in each context.
#'   List with the following structure: \code{clusterCounts = list(global=global, context=context)} where
#'   \code{global} is the number of global clusters, and \code{context} is the list of
#'   numbers of clusters in the individual contexts (datasets) of length \emph{C} where
#'   \code{context[c]} is the number of clusters in dataset \emph{c}.
#' @param dataDistributions Distribution of data in each dataset. Can be either a list of
#'   length \emph{C} where \code{dataDistributions[c]} is the distribution of dataset \emph{c},
#'   or a single string when all datasets have the same distribution. Currently implemented
#'   distribution is the \code{'diagNormal'} option for multivariate Normal distribution
#'   with diagonal covariance  matrix.
#' @param prior Prior distribution. If \code{NULL} then the prior is estimated using
#'   the datasets. The \code{'diagNormal'} distribution uses the Normal-Gamma distribution
#'   as a prior for each dimension.
#' @param maxIter Number of iterations of the Gibbs sampling algorithm.
#' @param burnin Number of burn-in iterations that will be discarded. If not specified,
#'   the algorithm discards the first half of the \code{maxIter} samples.
#' @param lag Used for thinning the samples.
#' @param verbose Print progress, by default \code{FALSE}.
#'
#' @return Returns list containing the sequence of MCMC states and the log likelihoods of
#'   the individual states.
#'   \item{samples}{List of assignments sampled from the posterior,
#'      each state \code{samples[[i]]} is a data frame with local and global assignments
#'      for each data point.}
#'   \item{logliks}{Log likelihoods during MCMC iterations.}
#'
#' @examples
#' # Example with simulated data (see vignette for details)
#' # Number of elements in each cluster
#' groupCounts <- c(50, 10, 40, 60)
#' # Centers of clusters
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
#' # Extract results from the samples
#' # Final state:
#' state <- results$samples[[length(results$samples)]]
#' # 1) assignment to global clusters
#' globalAssgn <- state$Global
#' # 2) context-specific assignmnets- assignment in specific dataset (context)
#' contextAssgn <- state[,"Context 1"]
#'
#' @importFrom magrittr %>%
#' @importFrom plyr aaply
#' @importFrom plyr laply
#' @importFrom plyr llply
#' @importFrom stats runif
#' @export

contextCluster <- function(datasets, clusterCounts,
     dataDistributions="diagNormal", prior=NULL,
     maxIter=1000, burnin=NULL, lag=3,
     verbose = FALSE){

  nDataPoints <- dim(datasets[[1]])[1]
  nContexts <- length(datasets)
  # check if the number of contexts = length(clusterCounts)
  if (nContexts != length(clusterCounts)) stop("Number of datasets is different from number of contexts.")

  # check if the number of data points is the same in all contexts
  if (sum(laply(datasets, function(dt) nrow(dt) == nDataPoints)) != nContexts) {
    stop("Number of data points is different in individual contexts.")
  }

  if(verbose) message('Running initialisation')
  # Distributions for every dataset (context)
  if (length(dataDistributions) == 1) {
    fullDataDistributions <- rep(dataDistributions, nContexts)
  } else {
    fullDataDistributions <- dataDistributions
  }

  # Default burnin interval
  if (is.null(burnin)) {
    burnin <- maxIter/2
  }

  # Create an initial state
  if (is.null(prior)) {
    prior <- empiricalBayesPrior(datasets, fullDataDistributions)
  }

  state <- createGibbsState(datasets, clusterCounts, fullDataDistributions)
  dataStats <- createDataStats(datasets, prior, fullDataDistributions)

  #******************
  # Run Gibbs sampling
  #samples <- vector("list", length = maxIter)
  logliks <- rep(0, maxIter)
  assignSamples <- vector("list", length=maxIter)
  for (iter in 1:maxIter) {
    if (((iter%%10) == 0 ) && verbose) { message(paste("iter ",iter)) }
    state <- state %>%
      gibbsSampleZ(dataStats, prior, clusterCounts) %>%
      gibbsSampleContextK(prior, clusterCounts)

    logliks[iter] <- logJoint(state, prior, clusterCounts)
    #samples[[iter]] <- state
    assignSamples[[iter]]$dataAssignments <- state$Z
    assignSamples[[iter]]$contextAssignments <- state$contextK
  }
  #****************
  
  thinned_logliks <- getMCMCSamples(logliks, burnin, lag)
  thinned_samples <- getMCMCSamples(assignSamples, burnin, lag)

  DIC <- computeDIC(thinned_logliks, thinned_samples, nDataPoints, 
                    state$distributions, clusterCounts, prior)
  
  assignments <-
    thinned_samples %>%
    llply(getClustersAssignments)

  return(list(samples=assignments, logliks=logliks, DIC=DIC))
}

#========= Helper functions ========================

logsumexp <- function(xs){
  xMax <- max(xs)
  values <- (xs - xMax) %>% exp %>% sum %>% log
  values + xMax
}

normalizeProbs <- function(logliks) {
  normalizer <- logsumexp(logliks)
  probs <- exp(logliks - normalizer)
  probs
}

sampleCategorical <- function(logliks) {
  normalizer <- logsumexp(logliks)
  probs <- exp(logliks - normalizer)
  cdf <- cumsum(probs)
  min(which(runif(1) < cdf))
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

# Function to combine cluster assignments into global clusters
generateMapping <- function(clusterCounts) {
  nGlobalClusters <- prod(clusterCounts)
  nContexts <- length(clusterCounts)
  mapping <- matrix(nrow=nGlobalClusters,ncol=nContexts+1)

  # indices of global clusters
  mapping[,nContexts+1] <- 1:nGlobalClusters

  # indices of local clusters
  localClusters <- list()
  for(c in 1:nContexts) localClusters[[c]] <- 1:clusterCounts[c]
  grid <- expand.grid(localClusters)
  for(c in 1:nContexts) mapping[,c] <- grid[[c]]

  mapping
}

getGlobalClusters <- function(state, mapping) {
  nDataPoints <- length(state$Z)
  nContexts <- length(state$distributions)
  clusterAssgn <- matrix(nrow=nDataPoints,ncol=nContexts)
  for (context in 1:nContexts) {
    clusterAssgn[,context] <- state$contextK[state$Z, context]
  }

  # Map to global clusters
  aaply(clusterAssgn,1,function(row)
    as.character(row) %>% paste(collapse=''))
}

saveSample <- function(state, mapping, iter, filename) {
  nDataPoints <- length(state$Z)
  nContexts <- length(state$distributions)

  if (!file.exists(filename)) {
    # Write header
    contextHeaders <-
      laply(1:nContexts, function(context) {
        laply(1:nDataPoints, function(n) paste("Context", context, " x", n, sep='')) %>% paste(collapse=',')
      }) %>% paste(collapse=',')
    # Global assignments
    globalHeaders1 <-
      laply(1:nDataPoints, function(n) paste("Global_assgn x", n, sep='')) %>% paste(collapse=',')
    # Underlying global clusters
    globalHeaders2 <-
      laply(1:nDataPoints, function(n) paste("Global x", n, sep='')) %>% paste(collapse=',')

    # Put all headers together
    header <- paste('Sample', contextHeaders, globalHeaders1, globalHeaders2, sep=',')

    fileConn<-file(filename,open = "a")
    writeLines(header, fileConn)
  } else {
    fileConn<-file(filename, open = "a")
  }

  # 1) context-specific clusters
  clusterAssgn = matrix(nrow=nDataPoints,ncol=nContexts)
  for (context in 1:nContexts) {
    clusterAssgn[,context] <- state$contextK[state$Z, context]
  }
  contextClusters <- clusterAssgn %>% as.vector %>% as.matrix

  # 2) global clusters
  assignments <- state$Z %>% as.matrix
  globalClusters <- getGlobalClusters(state, mapping)

  writeLines(c(iter, t(contextClusters), t(assignments), t(globalClusters)) %>% paste(collapse=','), fileConn)
  close(fileConn)
}

computeDIC <- function(thinned_logliks, thinned_samples, nDataPoints, distributions, clusterCounts, prior) {
  D_bar <- -2 * mean(thinned_logliks)   # -2 * average posterior log likelihoods
  theta_bar <- c() # mean of posterior parameter samples
  # take the posterior mode - most common assignment of each data point
  
  # most common global assignment for each data point
  modeDataAssignments <-
    1:nDataPoints %>%
    laply(function(id) {thinned_samples %>% laply(function(s) s$dataAssignments[id])}) %>%
    aaply(1, function(column) which.max(tabulate(column)))
  idx2d <- # helper variable to get indices into 2d array
    dim(thinned_samples[[1]]$contextAssignments) %>%
    lapply(seq) %>%
    expand.grid() %>%
    as.matrix
  modeContextAssignmentsTmp <-
    idx2d %>%
    aaply(1, function(x) {
      i <- x[1]; j <- x[2]
      thinned_samples %>%
        laply(function(s) s$contextAssignments[i,j]) %>%
        tabulate %>%
        which.max
    })
  # reshape results
  modeContextAssignments <-
    matrix(nrow=nrow(thinned_samples[[1]]$contextAssignments), ncol=ncol(thinned_samples[[1]]$contextAssignments))
  for (i in 1:nrow(modeContextAssignments)) {
    for (j in 1:ncol(modeContextAssignments)) {
      modeContextAssignments[i,j] <- modeContextAssignmentsTmp[idx2d[,1] == i & idx2d[,2] == j]
    }
  }
  
  modeState <- c()
  modeState$Z <- modeDataAssignments
  modeState$contextK <- modeContextAssignments
  modeState$distributions <- distributions
  modeState$clusterStats <-
    precomputeClusterStatistics(datasets, clusterCounts, modeState$Z, modeState$contextK, modeState$distributions)
  
  # Compute log likelihood of the resulting state
  D_hat <- -2 * logJoint(modeState, prior, clusterCounts)
  p_D <- D_bar - D_hat
  DIC <- p_D + D_bar
  DIC
}

#========= Initialization ================

precomputeClusterStatistics <- function(datasets, clusterCounts, Z, contextAssgn, fullDataDistributions) {
  N <- laply(1:clusterCounts$global, function(s) sum(Z == s))
  globalStats <- llply(1:length(datasets), function(context) {
    initGlobalClusterStatistics(Z, clusterCounts$global, datasets[[context]], fullDataDistributions[[context]])
  })
  globalStats$N <- N

  contextStats <- llply(1:length(datasets), function(context) {
    initContextClusterStatistics(Z, contextAssgn[,context], clusterCounts$context[context],
                                 datasets[[context]], fullDataDistributions[[context]])
  })
  list(globalStats = globalStats, contextStats = contextStats)
}

createDataStats <- function(datasets, prior, fullDataDistributions) {
  dataStats <- precomputeDataStatistics(datasets, fullDataDistributions)
  dataMarginals <- precomputeDataMarginals(dataStats, prior, fullDataDistributions)
  list(Stats = dataStats, Marginals = dataMarginals)
}

precomputeDataStatistics <- function(datasets, fullDataDistributions) {
  nContexts <- length(datasets)
  N <- nrow(datasets[[1]])
  llply(1:N, function(n) {
    llply(1:nContexts, function(context) {
      getDataStatistics(datasets[[context]][n,, drop=F], fullDataDistributions[[context]])
    })
  })
}

precomputeDataMarginals <- function(dataStatistics, prior, fullDataDistributions) {
  N <- length(dataStatistics)
  nContexts <- length(fullDataDistributions)
  llply(1:N, function(n) {
    llply(1:nContexts, function(context) {
      getMarginal(dataStatistics[[n]][[context]], 1, prior[[context]], fullDataDistributions[[context]])
    })
  })
}

# Generate a random initial state
createGibbsState <- function(datasets, clusterCounts, fullDataDistributions) {
  nContexts <- length(datasets)

  # Z ... assignment of data points to global clusters, with values 1..S
  N <- nrow(datasets[[1]])
  #Z <- sample(1:clusterCounts$global, N, replace=T)
  Z <- rep(1, N)

  # contextK ... assignment of global clusters to context clusters
  # => matrix of size [globalClusterCount, nContexts]
  # contextK[s,c] = index of cluster in context c for global cluster s
  contextK <- matrix(nrow=clusterCounts$global, ncol=nContexts)
  for (ic in 1:nContexts) {
    contextK[,ic] <- sample(1:clusterCounts$context[ic], clusterCounts$global, replace=T)
    #contextK[,ic] <- 1
  }

  # Sufficient statistics for clusters
  clusterStats <- precomputeClusterStatistics(datasets, clusterCounts, Z, contextK, fullDataDistributions)

  state <- list(Z = Z, contextK = contextK, distributions=fullDataDistributions,
                clusterStats = clusterStats)
  state
}


#========= Sampling Z =============================

getLogliksForGlobalCluster <- function(n, N, state, dataStats, prior) {
  S <- nrow(state$contextK)
  nContexts <- length(state$distributions)
  N_s <- state$clusterStats$globalStats$N
  N_s[state$Z[n]] <- N_s[state$Z[n]] - 1
  globalCurrentCluster <- state$Z[n]

  contextLogliks <- matrix(nrow=nContexts, ncol=S)
  for (context in 1:nContexts) {
    currentCluster <- state$contextK[globalCurrentCluster, context]
    clusterStats <- state$clusterStats$contextStats[[context]]

    cl <-
      predictiveLoglik(dataStats$Stats[[n]][[context]],
                       dataStats$Marginals[[n]][[context]], 1,
                       currentCluster, clusterStats,
                       prior[[context]], state$distributions[[context]])

    contextLogliks[context, ] <- cl[state$contextK[,context]]
  }

  logliks <- log(prior$gamma/S + N_s) - log(prior$gamma + N - 1) + colSums(contextLogliks)
  logliks
}


gibbsSampleZ <- function(state, dataStats, prior, clusterCounts) {
  N <- length(state$Z)
  nContexts <- length(state$distributions)

  #for (n in 1:N) {
  for (n in (sample(1:N,N,replace = F))) {   # randomized
    logliks <- getLogliksForGlobalCluster(n, N, state, dataStats, prior)
    probs <- normalizeProbs(logliks)
    newZ <- sample(length(probs), 1, prob=probs)
    oldZ <- state$Z[n]

    # update parameters
    if (newZ != oldZ) {
      # Update global cluster counts
      state$clusterStats$globalStats$N[newZ] <- state$clusterStats$globalStats$N[newZ] + 1
      state$clusterStats$globalStats$N[oldZ] <- state$clusterStats$globalStats$N[oldZ] - 1

      for (context in 1:nContexts) {
        state$clusterStats$globalStats[[context]] <-
          updateGlobalClusterStats(oldZ, newZ, state$clusterStats$globalStats[[context]], state$clusterStats$globalStats$N,
                                   dataStats$Stats[[n]][[context]], 1, state$distributions[[context]])


        oldK <- state$contextK[oldZ, context]
        newK <- state$contextK[newZ, context]

        if(oldK != newK) {
          state$clusterStats$contextStats[[context]] <-
            updateContextClusterStats(oldK, newK, state$clusterStats$contextStats[[context]],
                                      dataStats$Stats[[n]][[context]], 1, state$distributions[[context]])
        }
      }
      state$Z[n] <- newZ
    }
  }
  state
}

#========= Sampling K ===================================

getLogliksForContextCluster <- function(s, S, context, K, state, prior) {
  xStats <- getGlobalClusterStats(s, state$clusterStats$globalStats[[context]],
                                  state$distributions[[context]])
  clusterMarginal <- getMarginal(xStats, state$clusterStats$globalStats$N[s],
                                 prior[[context]], state$distribution[[context]])

  # Number of global clusters assigned to each local cluster
  m_l <- rep(0,K)
  for (l in 1:K) {
    kIdxs <- which(state$contextK[,context] == l)
    m_l[l] <- if (s %in% kIdxs) {length(kIdxs) - 1} else {length(kIdxs)}
  }

  clusterLoglik <-
    predictiveLoglik(xStats, clusterMarginal, state$clusterStats$globalStats$N[s],
                     state$contextK[s,context], state$clusterStats$contextStats[[context]],
                     prior[[context]], state$distributions[[context]])

  log(prior[[context]]$alpha/K + m_l) -
    log(prior[[context]]$alpha + S - 1) +
    clusterLoglik
}

gibbsSampleContextK <- function(state, prior, clusterCounts) {
  nContexts <- length(clusterCounts$context)
  S <- clusterCounts$global

  for (context in 1:nContexts) {
    K <- clusterCounts$context[context]
    #for (s in 1:S) {
    for (s in (sample(1:S, S, replace=F))) {
      logliks <- getLogliksForContextCluster(s, S, context, K, state, prior)
      probs <- normalizeProbs(logliks)
      oldK <- state$contextK[s,context]
      newK <- sample(length(probs), 1, prob=probs)

      # Update parameters
      if (oldK != newK) {
        globalClusterStats <- getGlobalClusterStats(s, state$clusterStats$globalStats[[context]], state$distributions[[context]])
        xN <- state$clusterStats$globalStats$N[s]

        contextClusterStats <- state$clusterStats$contextStats[[context]]
        contextClusterStats$xN <- state$clusterStats$contextStats[[context]]$xN

        state$clusterStats$contextStats[[context]] <-
          updateContextClusterStats(oldK, newK, contextClusterStats,
                                    globalClusterStats, xN, state$distributions[[context]])
      }

      state$contextK[s,context] <- newK
    }
  }
  state
}

#========= Helper functions for different distributions =========================
# - all distributions must implement the following functions

predictiveLoglik <- function(x, xMarginals, xN, currentCluster, clusterStats,
                             priorParams, distribution) {
  switch(distribution,
         diagNormal=predictiveLogLik_diagNormal(x, xMarginals, xN, currentCluster, clusterStats, priorParams))
}

getMarginal <- function(dataStats, dataN, priorParams, distribution) {
  switch(distribution,
         diagNormal=marginal_diagNormal(dataStats, dataN, priorParams))
}

getDataStatistics <- function(xs, distribution) {
  switch(distribution,
         diagNormal=getDataStatistics_diagNormal(xs))
}

initContextClusterStatistics <- function(Z, contextK, K, contextData, distribution) {
  switch(distribution,
         diagNormal = initContextClusterStatistics_diagNormal(Z, contextK, K, contextData))
}

initGlobalClusterStatistics <- function(Z, S, contextData, distribution) {
  switch(distribution,
         diagNormal = initGlobalClusterStatistics_diagNormal(Z, S, contextData))
}


updateContextClusterStats <- function(oldK, newK, clusterStats, dataStats, xN, distribution) {
  switch(distribution,
         diagNormal=updateContextClusterStats_diagNormal(oldK, newK, clusterStats, dataStats, xN))
}

updateGlobalClusterStats <- function(oldZ, newZ, clusterStats, clusterN, dataStats, xN, distribution) {
  switch(distribution,
         diagNormal=updateGlobalClusterStats_diagNormal(oldZ, newZ, clusterStats, clusterN, dataStats, xN))
}

getGlobalClusterStats <- function(globalClusterIdxs, clusterStats, distribution) {
  switch(distribution,
         diagNormal=getGlobalClusterStats_diagNormal(globalClusterIdxs, clusterStats))
}


getContextClusterMarginals <- function(clusterStats, prior, distribution) {
  switch(distribution,
         diagNormal = getContextClusterMarginals_diagNormal(clusterStats, prior))
}

#========= Normal distribution with diagonal covariance ====================

getParams_diagNormal <- function(sumX, sumX2, N, priorParams) {
  D <- length(priorParams$theta$m)
  S <- nrow(sumX)
  sXX <- sumX*sumX

  beta <- priorParams$theta$beta + N
  a <- priorParams$theta$a + N/2

  #   m1 <- (priorParams$theta$beta * 1/N * sXX )/
  #     rep.col(2*(priorParams$theta$beta + N), D)
  #   m2 <- (priorParams$theta$beta * 1/N * sXX ) *
  #     1/(2*(priorParams$theta$beta + N))
  invN <- 1/N
  invN[is.infinite(invN)] <- 0

  if (all(priorParams$theta$m == 0)) {
    dataContributions <- 0.5 * sumX2 - sXX * (1/2*invN) +
      (priorParams$theta$beta * invN * sXX ) * 1/(2*(priorParams$theta$beta + N))
    b <- sweep(dataContributions, 2, priorParams$theta$b, '+')
  } else {
    dataContributions <- 0.5 * sumX2 - rep.col(1/(2*N), D) * sXX +
      (priorParams$theta$beta * (
        rep.col(1/N, D) * sXX - 2*rep.row(priorParams$theta$m,S) * sumX +
          rep.col(N, D) * rep.row(priorParams$theta$m^2, S)))/rep.col(2*(priorParams$theta$beta + N), D)
    b <- sweep(dataContributions, 2, priorParams$theta$b, '+')
  }
  list(a = a, b = b, beta = beta)
}

predictiveLogLik_diagNormal <- function(xStats, xMarginals, xN, currentCluster,
                                        clusterStats, priorParams) {

  # Remove x from its current cluster
  clusterStats$xN[currentCluster] <- clusterStats$xN[currentCluster] - xN
  clusterStats$sumX[currentCluster,] <- clusterStats$sumX[currentCluster,] - xStats$sumX
  clusterStats$sumX2[currentCluster,] <- clusterStats$sumX2[currentCluster,] - xStats$sumX2

  paramsOther <- getParams_diagNormal(clusterStats$sumX, clusterStats$sumX2,
                                      clusterStats$xN, priorParams)
  paramsAll <- getParams_diagNormal(sweep(clusterStats$sumX, 2, xStats$sumX, '+'),
                                    sweep(clusterStats$sumX2, 2, xStats$sumX2, '+'),
                                    clusterStats$xN + xN, priorParams)
  D <- length(priorParams$theta$m)

  logliks <-
    D * (lgamma(paramsAll$a) - lgamma(paramsOther$a)) +
    paramsOther$a * rowSums(log(paramsOther$b)) -
    paramsAll$a * rowSums(log(paramsAll$b)) +
    D/2 * (log(paramsOther$beta) - log(paramsAll$beta)) -
    (xN * D)/2 * log(2*pi)
  logliks[clusterStats$xN == 0] <- xMarginals
  logliks
}

marginal_diagNormal <- function(xStats, xN, priorParams) {
  params <- getParams_diagNormal(xStats$sumX %>% as.matrix, xStats$sumX2 %>% as.matrix, xN, priorParams)
  D <- length(priorParams$theta$m)

  loglik <-
    D* (lgamma(params$a) - lgamma(priorParams$theta$a)) +
    D/2 * (log(priorParams$theta$beta) - log(params$beta)) -
    (xN * D)/2 * log(2*pi) +
    sum( priorParams$theta$a * log(priorParams$theta$b) - params$a * log(params$b))
  loglik[xN == 0] <- 0
  loglik
}

getDataStatistics_diagNormal <- function(xs) {
  list(sumX = t(as.matrix(colSums(xs))), sumX2 = t(as.matrix(colSums(xs*xs))))
}

initGlobalClusterStatistics_diagNormal <- function(Z, S, contextData) {
  sumX <- laply(1:S, function(s) contextData[Z == s,,drop=F] %>% colSums, .drop=F)
  sumX2 <- laply(1:S, function(s) {
    xs <- contextData[Z == s,,drop=F]
    (xs * xs) %>% colSums }, .drop=F)
  list(sumX = sumX, sumX2 = sumX2)
}

# Compute statistics for local clusters in a specific context
# inefficient but it's called only during initialisation
initContextClusterStatistics_diagNormal <- function(Z, contextK, K, contextData) {
  sumX <- laply(1:K, function(k) {
    globalIdxs <- which(contextK == k)
    contextData[Z %in% globalIdxs,,drop=F] %>% colSums
  }, .drop = F)
  sumX2 <- laply(1:K, function(k) {
    globalIdxs <- which(contextK == k)
    xs <- contextData[Z %in% globalIdxs,,drop=F]
    (xs * xs) %>% colSums
  }, .drop=F)
  xN <- laply(1:K, function(k) {
    globalIdxs <- which(contextK == k)
    sum(Z %in% globalIdxs)
  })

  list(sumX = sumX, sumX2 = sumX2, xN = xN)
}

updateContextClusterStats_diagNormal <- function(oldK, newK, clusterStats, dataStats, xN) {
  # 1) Add x to new cluster
  clusterStats$xN[newK] <- clusterStats$xN[newK] + xN
  clusterStats$sumX[newK,] <- clusterStats$sumX[newK,] + dataStats$sumX
  clusterStats$sumX2[newK,] <- clusterStats$sumX2[newK,] + dataStats$sumX2

  # 2) remove x from the previous cluster
  clusterStats$xN[oldK] <- clusterStats$xN[oldK] - xN
  if (clusterStats$xN[oldK] == 0) {
    clusterStats$sumX[oldK,] <- 0
    clusterStats$sumX2[oldK,] <- 0
  } else {
    clusterStats$sumX[oldK,] <- clusterStats$sumX[oldK,] - dataStats$sumX
    clusterStats$sumX2[oldK,] <- clusterStats$sumX2[oldK,] - dataStats$sumX2
  }
  clusterStats
}

updateGlobalClusterStats_diagNormal <- function(oldZ, newZ, clusterStats, clusterN, dataStats, xN) {
  clusterStats$sumX[newZ,] <- clusterStats$sumX[newZ,] + dataStats$sumX
  clusterStats$sumX2[newZ,] <- clusterStats$sumX2[newZ,] + dataStats$sumX2

  if (clusterN[oldZ] == 0) {
    clusterStats$sumX[oldZ,] <- 0
    clusterStats$sumX2[oldZ,] <- 0
  } else {
    clusterStats$sumX[oldZ,] <- clusterStats$sumX[oldZ,] - dataStats$sumX
    clusterStats$sumX2[oldZ,] <- clusterStats$sumX2[oldZ,] - dataStats$sumX2
  }
  clusterStats
}

getGlobalClusterStats_diagNormal <- function(s, clusterStats) {
  sumX <- clusterStats$sumX[s,,drop=F]
  sumX2 <- clusterStats$sumX2[s,,drop=F]
  list(sumX = sumX, sumX2 = sumX2)
}

#========= Joint ==============

logZ <- function(state, S, N, prior) {
  n_s <- state$clusterStats$globalStats$N
  lgamma(prior$gamma) - (S * lgamma(prior$gamma/S)) +
    sum(lgamma(prior$gamma/S + n_s)) - lgamma(prior$gamma + N)
}

logK <- function(state, context, K, prior) {
  m_l <- rep(0,K)
  for (l in 1:K) {
    m_l[l] <- sum(state$contextK[,context]==l)
  }

  M <- state$contextK %>% nrow
  lgamma(prior[[context]]$alpha) - K * lgamma(prior[[context]]$alpha/K) +
    sum(lgamma(prior[[context]]$alpha/K + m_l)) - lgamma(prior[[context]]$alpha + M)
}

logJoint <- function(state, prior, clusterCounts) {
  S <- clusterCounts$global
  N <- state$Z %>% length
  nContexts <- length(clusterCounts$context)

  # Data log likelihood
  dataLoglik <- 0
  loglikK <- 0
  for (context in 1:nContexts) {
    dataLoglik <-
      dataLoglik + getContextClusterMarginals(state$clusterStats$contextStats[[context]], prior[[context]],
                                              state$distributions[[context]]) %>%
      sum
    loglikK <- loglikK + logK(state, context, clusterCounts$context[context], prior)
  }


  logZ(state, S, N, prior) + loglikK +
    dataLoglik
}

getContextClusterMarginals_diagNormal <- function(clusterStats, priorParams) {
  params <- getParams_diagNormal(clusterStats$sumX, clusterStats$sumX2,
                                 clusterStats$xN, priorParams)
  D <- length(priorParams$theta$m)

  loglik <-
    D * (lgamma(params$a) - lgamma(priorParams$theta$a)) +
    D/2 * (log(priorParams$theta$beta) - log(params$beta)) -
    (clusterStats$xN * D)/2 * log(2*pi) +
    rowSums( priorParams$theta$a * log(priorParams$theta$b) - params$a * log(params$b))
  loglik[clusterStats$xN == 0] <- 0
  loglik
}
