context("testing Clustering with Matching")

test_that("test that matchClustering() works and outputs the same results when used in parallel and sequentially", {

  ns <- 10
  runs <- 50
  features <- 20
  run1 <- 1
  run2 <- 2
  p_boot <- matrix(runif(features*runs*ns),nrow = features)
  
  message("building distMatrix...")
  t1 <- Sys.time()
  distMatrix <- 1 - computeCorrelation(p_boot)
  t2 <- Sys.time()
  secs_building <- as.double(difftime(t2,t1),units="secs")
  
  message("running matched clustering parallel=FALSE...")
  t1 <- Sys.time()
  m1 <- matchedClustering(distMatrix,ns,maxMatch = TRUE)
  t2 <- Sys.time()
  secs_m1 <- as.double(difftime(t2,t1),units="secs")
  
  message("running matched clustering parallel=TRUE...")
  t1 <- Sys.time()
  m2 <- matchedClustering(distMatrix,ns,maxMatch = TRUE,parallel = TRUE,nparallel = 4)
  t2 <- Sys.time()
  secs_m2 <- as.double(difftime(t2,t1),units="secs")
  
  message("building distMatrix: ",sprintf("%.1f",secs_building)," seconds")
  message("sequential: ",sprintf("%.1f",secs_m1)," seconds")
  message("parallel: ",sprintf("%.1f",secs_m2)," seconds")
  message("saved ",sprintf("%.1f",(secs_m1-secs_m2)/secs_m1*100),"% of time")
  message("parallel and sequential solutions are identical: ",all(m1==m2))
  expect_equal( m1, m2 )
})