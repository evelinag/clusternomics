
#============================================================
# Global behaviour - tests
#============================================================

test_that("Clusternomics identifies well-separated clusters",{
  set.seed(1)
  groupCounts <- c(100, 70, 5, 100)
  means <- c(-2.5,2.5)
  testData <- generateTestData_2D(groupCounts, means)
  datasets <- testData$data
  fullDataDistributions <- rep('diagNormal', 2)
  clusterCounts <- list(global=10, context=c(2,2))
  nContexts <- length(clusterCounts$context)

  result <- contextCluster(datasets,
              clusterCounts, "diagNormal", maxIter=500,
              prior=NULL, verbose = T)

  samples <- result$samples
  clusters <-
    laply(1:length(samples),
          function(i) samples[[i]]$Global)

  # Correct number of clusters
  k <-
    numberOfClusters(clusters) %>%
      table %>%
      which.max %>%
      names %>%
      as.numeric
  expect_equal(k, sum(groupCounts > 0))

  # Correct sizes of clusters
  assgnEstim <-
    clusters[nrow(clusters),] %>%
    table %>% as.numeric %>% sort
  assgnTrue <-
    testData$groups %>%
    table %>% as.numeric %>% sort
  expect_equal(assgnEstim, assgnTrue)
})



test_that("BIC is lower for better model",{
  set.seed(1)
  groupCounts <- c(100, 70, 5, 100)
  means <- c(-2.5,2.5)
  testData <- generateTestData_2D(groupCounts, means)
  datasets <- testData$data
  fullDataDistributions <- rep('diagNormal', 2)
  clusterCounts <- list(global=10, context=c(2,2))
  nContexts <- length(clusterCounts$context)

  result <- contextCluster(datasets,
                           clusterCounts, "diagNormal", maxIter=500,
                           prior=NULL, verbose = T)
  resultWorse <- contextCluster(datasets,
                                list(global=2, context=c(2,1)), "diagNormal", maxIter=500, prior=NULL, verbose = T)

  expect_lt(result$DIC, resultWorse$DIC)
})
