# Various test for the following:
#
# a) Testing that more efficient implementation yields the same results as 
#    slow step-by-step computation
#
# b) Gibbs sampling is consistent (conditional distributions are consistent
#    with the joint distribution)


predictiveLogLik <- function(x, xN, clusterStats, clusterN, containsX, priorParams, distribution) {
  switch(distribution,
         diagNormal=predictiveLogLik_diagNormal(x, xN, clusterStats, clusterN, containsX, priorParams))  
}

predictiveLogLik_diagNormal <- function(x, xN, clusterStats, clusterN, containsX, priorParams) {
  D <- x %>% ncol
  rownames(x) <- NULL
  
  # I shouldn't need to compute this
  sumX <- colSums(x) 
  sumX2 <- colSums(x*x)
  
  if (clusterN == 0 || (containsX && xN == clusterN)) {
    marginal_diagNormal(list(sumX, sumX2), xN, priorParams)
  } else {
    if (!containsX) {
      paramsOther <- getParams_diagNormal(clusterStats$sumX, clusterStats$sumX2, clusterN, priorParams)
      paramsAll <- getParams_diagNormal(clusterStats$sumX + sumX, 
                                        clusterStats$sumX2 + sumX2, clusterN + xN, priorParams)
    } else {
      paramsOther <- getParams_diagNormal(clusterStats$sumX - sumX, clusterStats$sumX2 - sumX2,
                                          clusterN - xN, priorParams)
      paramsAll <- getParams_diagNormal(clusterStats$sumX, clusterStats$sumX2, clusterN, priorParams)
    }
    
    D * (lgamma(paramsAll$a) - lgamma(paramsOther$a)) + 
      paramsOther$a * sum(log(paramsOther$b)) - 
      paramsAll$a * sum(log(paramsAll$b)) + 
      D * 0.5 * (log(paramsOther$beta) - log(paramsAll$beta)) - (xN*D)/2 * log(2*pi)
  }  
}

getParams_diagNormal <- function(sumX, sumX2, N, priorParams) {
  beta <- priorParams$theta$beta + N
  a <- priorParams$theta$a + N/2
  sXX <- sumX * sumX
  
  if (N > 0) {
    b <- priorParams$theta$b + 0.5 * sumX2 - 1/(2*N) * sXX + 
      (priorParams$theta$beta *(1/N * sXX - 2*priorParams$theta$m * sumX + N * priorParams$theta$m^2))/
      (2*(priorParams$theta$beta + N))
  } else {
    b <- priorParams$theta$b
  }
  list(a = a, b = b, beta = beta)  
}

test_that("Log likelihoods for sampling z are correct", {
  groupCounts <- c(100, 70, 0, 100)
  means <- c(-1.5,1.5)
  testData <- generateTestData_2D(groupCounts, means)
  datasets <- testData$data
  fullDataDistributions = rep('diagNormal', 2)
  clusterCounts = list(global=10, context=c(2,2))
  nContexts <- length(clusterCounts$context)
  
  prior <- generatePrior(datasets, fullDataDistributions)
  prior$gamma <- 3.0
  for (context in 1:nContexts) {
    prior[[context]]$alpha <- 2.0
  }
  state <- createGibbsState(datasets, clusterCounts, fullDataDistributions)
  dataStats <- createDataStats(datasets, prior, fullDataDistributions)
  # ^^^^ end initialisation
  
  getParams <- function(x, priorParams) {
    N <- nrow(x)
    beta <- priorParams$theta$beta + N
    a <- priorParams$theta$a + N/2
    avgX <- colMeans(x)
    diff2 <- aaply(x, 1, function(row) (row - avgX)^2) %>% rbind %>% colSums
    b <- priorParams$theta$b + 
      0.5 * diff2 +
      (priorParams$theta$beta * N * (avgX - priorParams$theta$m)^2)/(2*(priorParams$theta$beta + N))
    list(a = a, b = b, beta = beta)  
  }
  
  marginal <- function(x, priorParams) {
    D <- x %>% ncol
    params <- getParams(x, priorParams)
    D * (lgamma(params$a) - lgamma(priorParams$theta$a)) +
      D * 0.5 * (log(priorParams$theta$beta) - log(params$beta)) -
      (nrow(x) * D)/2 * log(2*pi) + 
      sum( priorParams$theta$a * log(priorParams$theta$b) - params$a * log(params$b))
  }
  
  ploglik <- function(x, otherX, priorParams, distribution) {
    D <- x %>% ncol
    rownames(x) <- NULL
    m <- nrow(x)
    
    if (nrow(otherX) == 0) {
      marginal(x, priorParams)
    } else {
      paramsOther <- getParams(otherX, priorParams)
      paramsAll <- getParams(rbind2(otherX,x), priorParams)
      
      D * (lgamma(paramsAll$a) - lgamma(paramsOther$a)) + 
        paramsOther$a * sum(log(paramsOther$b)) - 
        paramsAll$a * sum(log(paramsAll$b)) + 
        D * 0.5 * (log(paramsOther$beta) - log(paramsAll$beta)) - (m*D)/2 * log(2*pi)
    }
  }
  
  # vvvv compute log likelihoods
  zSamples <- state$Z
  N <- length(zSamples)
  for (n in 1:N) { 
    S <- clusterCounts$global
    nContexts <- length(clusterCounts$context)
    
    logliks <- vector("numeric", S)
    for (s in 1:S) {   
      n_s <- sum(zSamples == s)
      if (zSamples[n] == s) n_s <- n_s - 1
      
      loglikfx_context <- rep(0, nContexts)
      for (context in 1:nContexts) {
        globalIdxs <- which(state$contextK[,context] == state$contextK[s, context])
        zIdxs <- zSamples %in% globalIdxs
        zIdxs[n] <- FALSE    

        x <- datasets[[context]][n,] %>% rbind  # create row vector
        
        otherX <- datasets[[context]][zIdxs, ,drop=F] %>% rbind
        
        loglikfx_context[context] <- ploglik(x, otherX, prior[[context]], state$distributions[[context]])
      }      
      
      logliks[s] <- 
        log(prior$gamma/clusterCounts$global + n_s) -
        log(prior$gamma + N - 1) + 
        sum(loglikfx_context)
    }
    expect_equal(logliks, getLogliksForGlobalCluster(n, N, state, dataStats, prior) %>% as.numeric,
                 tolerance=1e-10)
  }
})


getSamplingDistribution <- function(n, N, S, nContexts, state, prior, datasets) {
  logliks <- rep(0, S)
  for (s in 1:S) {
    gamma0 <- prior$gamma
    n_s <- 0
    for (nn in 1:N) {
      if ((nn != n) && (state$Z[nn] == s))  {
        n_s <- n_s + 1
      } 
    }
    
    contextLogliks <- rep(0, nContexts)
    for (context in 1:nContexts) {
      x <- datasets[[context]][n,, drop=F]
      D <- ncol(x)
      
      kIdxs <- which(state$contextK[,context] == state$contextK[s,context])
      zIdxs <- state$Z %in% kIdxs
      zIdxs[n] <- F
      otherXs <- datasets[[context]][zIdxs,, drop=F]
      xN <- nrow(otherXs)
      
      a0 <- prior[[context]]$theta$a
      b0 <- prior[[context]]$theta$b
      beta0 <- prior[[context]]$theta$beta
      m0 <- prior[[context]]$theta$m
      
      aOthers <- a0 + xN/2
      betaOthers <- beta0 + xN
      avgOthers <- colMeans(otherXs)
      if (xN > 0) {
        diffOthers <- laply(1:xN, function(o) (otherXs[o,,drop=F] - avgOthers)^2 ) %>% rbind %>% colSums
        bOthers <- b0 + 0.5 * diffOthers + (beta0 * xN * (avgOthers - m0)^2)/(2*(beta0 + xN))
      } else {
        bOthers <- b0
      }
      
      xAll <- rbind(otherXs, x)
      aAll <- a0 + (xN + 1)/2
      betaAll <- beta0 + xN + 1
      avgAll <- colMeans(xAll)
      diffAll <- laply(1:(xN+1), function(o) (xAll[o,,drop=F] - avgAll)^2 ) %>% rbind %>% colSums
      bAll <- b0 + 0.5 * diffAll + (beta0 * (xN + 1) * (avgAll - m0)^2)/(2*(beta0 + xN + 1))
      
      logpOtherXs <- D * (lgamma(aOthers) - lgamma(a0)) + sum(a0 * log(b0) - aOthers * log(bOthers)) + 
        D * 0.5 * (log(beta0) - log(betaOthers)) - (xN * D)/2 * log(2*pi)
      
      logpAllXs <- D * (lgamma(aAll) - lgamma(a0)) + sum(a0 * log(b0) - aAll * log(bAll)) + 
        D * 0.5 * (log(beta0) - log(betaAll)) - ((xN + 1) * D)/2 * log(2*pi)
      
      contextLogliks[context] <- logpAllXs - logpOtherXs
    }
    logliks[s] <- log(gamma0/S + n_s) - log(gamma0 + N - 1) + sum(contextLogliks)
  }
  logliks
}


test_that("Log likelihoods for sampling z are correct (v2)", {
  set.seed(4)
  groupCounts <- c(100, 70, 0, 100)
  means <- c(-1.5,1.5)
  testData <- generateTestData_2D(groupCounts, means)
  datasets <- testData$data
  fullDataDistributions = rep('diagNormal', 2)
  clusterCounts = list(global=4, context=c(2,2))
  nContexts <- length(clusterCounts$context)
  
  prior <- generatePrior(datasets, fullDataDistributions)
  prior$gamma <- 3.0
  for (context in 1:nContexts) {
    prior[[context]]$alpha <- 2.0
  }
  state <- createGibbsState(datasets, clusterCounts, fullDataDistributions)
  dataStats <- createDataStats(datasets, prior, fullDataDistributions)
  
  N <- sum(groupCounts)
  S <- clusterCounts$global
  for (n in 1:N) {
    logliks <- getSamplingDistribution(n, N, S, nContexts, state, prior, datasets)
    expect_equal(logliks, getLogliksForGlobalCluster(n, N, state, dataStats, prior) %>% as.numeric, 
                 tolerance=1e-10)
  }
  
})

testLogJoint <- function(state, clusterCounts, datasets, prior) {
  nContexts <- length(clusterCounts$context)
  # Compute joint log likelihood - explicitly
  # 1) log p(Z) 
  S <- clusterCounts$global
  gamma <- prior$gamma
  n_s <- rep(0, S)
  N <- length(state$Z)
  for (n in 1:N) n_s[state$Z[n]] <- n_s[state$Z[n]] + 1
  
  loglikZ <- 
    lgamma(gamma) - (S * lgamma(gamma/S)) + 
    sum(lgamma(gamma/S + n_s)) - lgamma(gamma + N)
  
  # 2) log p(K)
  loglikK <- 
    laply(1:nContexts, function(context) {
      alpha <- prior[[context]]$alpha
      K <- clusterCounts$context[context]
      m_l <- laply(1:K, function(l) {
        sum(state$contextK[,context] == l)
      })
      lgamma(alpha) - K * lgamma(alpha/K) + 
        sum(lgamma(alpha/K + m_l)) - lgamma(alpha + S)
    }) %>% sum
  
  # 3) log p(X|K,Z)
  loglikX <-
    laply(1:nContexts, function(context) {
      laply(1:clusterCounts$context[context], function(k) {
        globalIdxs <- which(state$contextK[,context] == k)
        xs <- datasets[[context]][state$Z %in% globalIdxs,, drop=F]
        xN <- nrow(xs)
        D <- ncol(xs)
        
        # Compute marginal likelihood (assume Gaussian)
        a0 <- prior[[context]]$theta$a
        b0 <- prior[[context]]$theta$b
        m0 <- prior[[context]]$theta$m
        beta0 <- prior[[context]]$theta$beta
        
        a <- a0 + xN/2
        beta <- beta0 + xN
        
        if (xN == 0) {
          xMean <- rep(0, length(m0)) 
          diff <- rep(0, length(m0))
        } else {
          xMean <- colMeans(xs)
          diff <- colSums((xs - rep.row(xMean, xN))^2)
        }
        
        b <- b0 + 0.5 * diff +
          (beta0 * xN * (xMean - m0)^2)/(2 * (beta0 + xN))
        
        D*(lgamma(a) - lgamma(a0)) + sum(a0 * log(b0) - a * log(b)) +
          D * 0.5 * (log(beta0) - log(beta)) - (xN * D)/2 * log(2*pi)
      }) %>% sum
    }) %>% sum
  
  jointLoglik <- loglikZ + loglikK + loglikX 
  jointLoglik
}



test_that("Gibbs sampling is consistent", {
  groupCounts <- c(50, 70, 0, 30)
  means <- c(-1.5,1.5)
  testData <- generateTestData_2D(groupCounts, means)
  datasets <- testData$data
  fullDataDistributions = rep('diagNormal', 2)
  clusterCounts = list(global=10, context=c(7,9))
  nContexts <- length(clusterCounts$context)
  
  prior <- generatePrior(datasets, fullDataDistributions)
  prior$gamma <- 3.0
  for (context in 1:nContexts) {
    prior[[context]]$alpha <- 2.0
  }
  state <- createGibbsState(datasets, clusterCounts, fullDataDistributions)
  dataStats <- createDataStats(datasets, prior, fullDataDistributions)
  
  N <- sum(groupCounts)
  S <- clusterCounts$global
  
  
  getJointDistribution <- function(st) { testLogJoint(st, clusterCounts, datasets, prior) }
  getCondDistribution <- function(n, st) { getSamplingDistribution(n, N, S, nContexts, st, prior, datasets)}
    
  for (n in 1:N) {
    stateA <- state
    stateB <- state
    newZ <- sample((1:clusterCounts$global)[-stateB$Z[n]], 1)  
    stateB$Z[n] <- newZ
    
    p1 <- getCondDistribution(n, stateA)
    p2 <- getCondDistribution(n, stateB)
    q1 <- getJointDistribution(stateA)
    q2 <- getJointDistribution(stateB)
    
    probsA <- normalizeProbs(p1)
    probsB <- normalizeProbs(p2)
    condDiff <- probsA[stateA$Z[n]]/probsB[stateB$Z[n]]

    expect_equal(condDiff, exp(q1-q2))
  }
})

test_that("Log likelihoods for sampling k are correct", {
  set.seed(2)
  groupCounts <- c(100, 70, 0, 100)
  means <- c(-1.5,1.5)
  testData <- generateTestData_2D(groupCounts, means)
  datasets <- testData$data
  fullDataDistributions = rep('diagNormal', 2)
  clusterCounts = list(global=4, context=c(2,2))
  nContexts <- length(clusterCounts$context)
  
  prior <- generatePrior(datasets, fullDataDistributions)
  prior$gamma <- 3.0
  for (context in 1:nContexts) {
    prior[[context]]$alpha <- 2.0
  }
  state <- createGibbsState(datasets, clusterCounts, fullDataDistributions)
  # ^^^^ end initialisation
  
  getParams <- function(x, priorParams) {
    N <- nrow(x)
    beta <- priorParams$theta$beta + N
    a <- priorParams$theta$a + N/2
    avgX <- colMeans(x)
    diff2 <- aaply(x, 1, function(row) (row - avgX)^2) %>% rbind %>% colSums
    b <- priorParams$theta$b + 
      0.5 * diff2 +
      (priorParams$theta$beta * N * (avgX - priorParams$theta$m)^2)/(2*(priorParams$theta$beta + N))
    list(a = a, b = b, beta = beta)  
  }
  
  marginal <- function(x, priorParams) {
    params <- getParams(x, priorParams)
    D <- ncol(x)
    D * (lgamma(params$a) - lgamma(priorParams$theta$a)) +
      D * 0.5 * (log(priorParams$theta$beta) - log(params$beta)) -
      (D * nrow(x))/2 * log(2*pi) + 
      sum( priorParams$theta$a * log(priorParams$theta$b) - params$a * log(params$b))
  }
  
  ploglik <- function(x, otherX, priorParams, distribution) {
    D <- x %>% ncol
    rownames(x) <- NULL
    m <- nrow(x)
    
    if (nrow(otherX) == 0) {
      marginal(x, priorParams)
    } else {
      paramsOther <- getParams(otherX, priorParams)
      paramsAll <- getParams(rbind2(otherX,x), priorParams)
      
      D * (lgamma(paramsAll$a) - lgamma(paramsOther$a)) + 
        paramsOther$a * sum(log(paramsOther$b)) - 
        paramsAll$a * sum(log(paramsAll$b)) + 
        D * 0.5 * (log(paramsOther$beta) - log(paramsAll$beta)) - (D*m)/2 * log(2*pi)
    }
  }
  
  nContexts <- length(clusterCounts$context)
  S <- clusterCounts$global
  
  for (context in 1:nContexts) {
    K <- clusterCounts$context[context]
    for (s in 1:S) {
      zIdxs <- state$Z == s
      x <- datasets[[context]][zIdxs, ,drop=F]
      
      logliks <- vector("numeric", K)
      for (l in 1:K) {
        kIdxs <- which(state$contextK[,context] == l)
        zkIdxs <- (laply(kIdxs, function(k) state$Z == k) %>% rbind %>% colSums ) > 0
        oIdxs <- zkIdxs & !zIdxs
        otherX <- datasets[[context]][oIdxs,]
        
        if (s %in% kIdxs) {
          m_l <- length(kIdxs) - 1
        } else {
          m_l <- length(kIdxs)
        }
            
        clusterLoglik <- ploglik(x, otherX, prior[[context]], state$distribution[[context]])
        logliks[l] <- 
          log(prior[[context]]$alpha/K + m_l) - 
          log(prior[[context]]$alpha + S - 1) + 
          clusterLoglik
      }
      
      expect_equal(logliks, getLogliksForContextCluster(s, S, context, K, state, prior))
    }
  }
  
})


test_that("Cluster statistics are consistent after sampling Z", {
  groupCounts <- c(100, 70, 0, 100)
  means <- c(-1.5,1.5)
  testData <- generateTestData_2D(groupCounts, means)
  datasets <- testData$data
  fullDataDistributions = rep('diagNormal', 2)
  clusterCounts = list(global=10, context=c(3,6))
  nContexts <- length(clusterCounts$context)
  
  prior <- generatePrior(datasets, fullDataDistributions)
  state <- createGibbsState(datasets, clusterCounts, fullDataDistributions)
  dataStats <- createDataStats(datasets, prior, fullDataDistributions)

  newState <- gibbsSampleZ(state, dataStats, prior, clusterCounts)
  expect_equal(newState$clusterStats, 
               precomputeClusterStatistics(datasets, clusterCounts, newState$Z, newState$contextK, fullDataDistributions) )
})

test_that("Cluster statistics are consistent after sampling K", {
  groupCounts <- c(100, 70, 0, 100)
  means <- c(-1.5,1.5)
  testData <- generateTestData_2D(groupCounts, means)
  datasets <- testData$data
  fullDataDistributions = rep('diagNormal', 2)
  clusterCounts = list(global=10, context=c(3,6))
  nContexts <- length(clusterCounts$context)
  
  prior <- generatePrior(datasets, fullDataDistributions)
  prior$gamma <- 3.0
  for (context in 1:nContexts) {
    prior[[context]]$alpha <- 2.0
  }
  state <- createGibbsState(datasets, clusterCounts, fullDataDistributions)
  dataStats <- createDataStats(datasets, prior, fullDataDistributions)
  
  newState <- gibbsSampleContextK(state, prior, clusterCounts)
  expect_equal(newState$clusterStats, 
               precomputeClusterStatistics(datasets, clusterCounts, newState$Z, newState$contextK, fullDataDistributions) )
})

test_that("Joint log likelihood is correct", {
  groupCounts <- c(100, 70, 0, 100)
  means <- c(-1.5,1.5)
  testData <- generateTestData_2D(groupCounts, means)
  datasets <- testData$data
  fullDataDistributions = rep('diagNormal', 2)
  clusterCounts = list(global=100, context=c(30,100))
  nContexts <- length(clusterCounts$context)
  verbose = TRUE
  prior = NULL
  previousState = NULL
  max.iter = 20
  verbose = FALSE
  
  prior <- generatePrior(datasets, fullDataDistributions)
  prior$gamma <- 3.0
  for (context in 1:nContexts) {
    prior[[context]]$alpha <- 2.0
  }
  state <- createGibbsState(datasets, clusterCounts, fullDataDistributions)
  
  jointLoglik <- testLogJoint(state, clusterCounts, datasets, prior)
  
  expect_equal(jointLoglik, logJoint(state, prior, clusterCounts))
})

test_that("Sampling distribution for Z is consistent with the joint", {
  set.seed(1)
  groupCounts <- c(100, 70, 0, 100)
  means <- c(-1.5,1.5)
  testData <- generateTestData_2D(groupCounts, means)
  datasets <- testData$data
  fullDataDistributions = rep('diagNormal', 2)
  clusterCounts = list(global=4, context=c(2,2))
  nContexts <- length(clusterCounts$context)
  
  prior <- generatePrior(datasets, fullDataDistributions)
  prior$gamma <- 3.0
  for (context in 1:nContexts) {
    prior[[context]]$alpha <- 2.0
  }
  state <- createGibbsState(datasets, clusterCounts, fullDataDistributions)
  dataStats <- createDataStats(datasets, prior, fullDataDistributions)
  
  # sample a new assignment for a random data point
  n <- sample(1:sum(groupCounts), 1)
  stateA <- state
  stateB <- state
  stateB$Z[n] <- sample((1:clusterCounts$global)[-stateB$Z[n]], 1)
  stateB$clusterStats <- precomputeClusterStatistics(datasets, clusterCounts, stateB$Z, stateB$contextK, fullDataDistributions)
    
  jointDiff <- testLogJoint(stateA, clusterCounts, datasets, prior) - 
    testLogJoint(stateB, clusterCounts, datasets, prior)  
  
  N <- sum(groupCounts)
  logCondA <- getLogliksForGlobalCluster(n, N, stateA, dataStats, prior)
  logCondB <- getLogliksForGlobalCluster(n, N, stateB, dataStats, prior)

  probsA <- normalizeProbs(logCondA)
  probsB <- normalizeProbs(logCondB)
  condDiff <- probsA[stateA$Z[n]]/probsB[stateB$Z[n]]
  
  expect_equal(condDiff %>% as.numeric, exp(jointDiff), tolerance=1e-10)
})

test_that("Sampling distribution for K is consistent with the joint", {
  groupCounts <- c(100, 70, 0, 100)
  means <- c(-1.5,1.5)
  testData <- generateTestData_2D(groupCounts, means)
  datasets <- testData$data
  fullDataDistributions = rep('diagNormal', 2)
  clusterCounts = list(global=4, context=c(2,2))
  nContexts <- length(clusterCounts$context)
  
  prior <- generatePrior(datasets, fullDataDistributions)
  prior$gamma <- 3.0
  for (context in 1:nContexts) {
    prior[[context]]$alpha <- 2.0
  }
  state <- createGibbsState(datasets, clusterCounts, fullDataDistributions)
  dataStats <- createDataStats(datasets, prior, fullDataDistributions)
  
  # sample a new assignment for a random global cluster
  s <- sample(1:clusterCounts$global, 1)
  stateA <- state
  stateB <- state
  context <- 1
  
  stateB$contextK[s,context] <- sample((1:clusterCounts$context[context])[-stateB$contextK[s, context]], 1)
  stateB$clusterStats <- precomputeClusterStatistics(datasets, clusterCounts, stateB$Z, stateB$contextK, fullDataDistributions)
  
  jointDiff <- testLogJoint(stateA, clusterCounts, datasets, prior) - 
    testLogJoint(stateB, clusterCounts, datasets, prior)
  
  logCondA <- getLogliksForContextCluster(s, clusterCounts$global, context, clusterCounts$context[context], stateA, prior)
  logCondB <- getLogliksForContextCluster(s, clusterCounts$global, context, clusterCounts$context[context], stateB, prior)
  probsA <- normalizeProbs(logCondA)
  probsB <- normalizeProbs(logCondB)
  condDiff <- probsA[stateA$contextK[s,context]]/probsB[stateB$contextK[s,context]]
  
  expect_equal(condDiff, exp(jointDiff), tolerance=1e-10)
})
