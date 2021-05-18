context("testing functions Fit and FitMS, should take a few seconds")

test_that("test that Fit() works on a random matrix, using useBootstrap=FALSE", {
  #test substitutions
  mut.order <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
  n_row <- 96
  n_col <- 15
  rnd_matrix <- round(matrix(runif(n_row*n_col,min = 0,max = 50),nrow = n_row,ncol = n_col))
  colnames(rnd_matrix) <- paste0("C",1:n_col)
  row.names(rnd_matrix) <- mut.order
  res <- Fit(catalogues = rnd_matrix,
             signatures = cosmic30,
             useBootstrap = FALSE,
             threshold_percent = 0.1)
  
  #here goes the test
  expect_true(TRUE)
})

test_that("test that Fit() works on a random matrix, using useBootstrap=TRUE", {
  #test substitutions
  mut.order <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
  n_row <- 96
  n_col <- 15
  rnd_matrix <- round(matrix(runif(n_row*n_col,min = 0,max = 50),nrow = n_row,ncol = n_col))
  colnames(rnd_matrix) <- paste0("C",1:n_col)
  row.names(rnd_matrix) <- mut.order
  res <- Fit(catalogues = rnd_matrix,
             signatures = cosmic30,
             useBootstrap = TRUE,
             nboot = 5,
             threshold_percent = 0.1,
             threshold_p.value = 0.1)

  #here goes the test
  expect_true(TRUE)
})

test_that("test that Fit() works on a random matrix with useBootstrap=TRUE, and that random seed works, no parallel", {
  #test substitutions
  mut.order <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
  n_row <- 96
  n_col <- 5
  rnd_matrix <- round(matrix(runif(n_row*n_col,min = 0,max = 50),nrow = n_row,ncol = n_col))
  colnames(rnd_matrix) <- paste0("C",1:n_col)
  row.names(rnd_matrix) <- mut.order
  res1 <- Fit(catalogues = rnd_matrix,
              signatures = cosmic30,
              useBootstrap = TRUE,
              nboot = 5,
              threshold_percent = 0,
              threshold_p.value = 0.1,
              nparallel = 1,
              randomSeed = 1)
  res2 <- Fit(catalogues = rnd_matrix,
              signatures = cosmic30,
              useBootstrap = TRUE,
              nboot = 5,
              threshold_percent = 0,
              threshold_p.value = 0.1,
              nparallel = 1,
              randomSeed = 1)
  
  #here goes the test
  expect_equal(res1$exposures,res2$exposures)
})

test_that("test that Fit() works on a random matrix with useBootstrap=TRUE, and that random seed works, using parallel", {
  #test substitutions
  mut.order <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
  n_row <- 96
  n_col <- 5
  rnd_matrix <- round(matrix(runif(n_row*n_col,min = 0,max = 50),nrow = n_row,ncol = n_col))
  colnames(rnd_matrix) <- paste0("C",1:n_col)
  row.names(rnd_matrix) <- mut.order
  res1 <- Fit(catalogues = rnd_matrix,
              signatures = cosmic30,
              useBootstrap = TRUE,
              nboot = 5,
              threshold_percent = 0,
              threshold_p.value = 0.1,
              nparallel = 2,
              randomSeed = 1)
  res2 <- Fit(catalogues = rnd_matrix,
              signatures = cosmic30,
              useBootstrap = TRUE,
              nboot = 5,
              threshold_percent = 0,
              threshold_p.value = 0.1,
              nparallel = 2,
              randomSeed = 1)
  
  #here goes the test
  expect_equal(res1$exposures,res2$exposures)
})

test_that("test that plotFit() works on a random matrix", {
  #test substitutions
  mut.order <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
  n_row <- 96
  n_col <- 15
  rnd_matrix <- round(matrix(runif(n_row*n_col,min = 0,max = 50),nrow = n_row,ncol = n_col))
  colnames(rnd_matrix) <- paste0("C",1:n_col)
  row.names(rnd_matrix) <- mut.order
  res <- Fit(catalogues = rnd_matrix,
             signatures = cosmic30,
             useBootstrap = TRUE,
             nboot = 5,
             threshold_percent = 0.1,
             threshold_p.value = 0.1)
  plotFit(res,outdir = "test_fit/")
  
  unlink("test_fit/", recursive = TRUE)
  
  #here goes the test
  expect_true(TRUE)
})

test_that("test that FitMS() with multiStepMode=constrainedFit works with useBootstrap=FALSE then plots", {
  #test substitutions
  catalogues <- read.table("FitAndFitMStests/FitandFitMStestCatalogues.tsv",sep = "\t",check.names = F,header = T,stringsAsFactors = F)
  res <- FitMS(catalogues = catalogues,
               organ = "Breast",
               multiStepMode = "constrainedFit",
               threshold_percent = 0.1)
  plotFitMS(res,outdir = "test_fit/")
  
  unlink("test_fit/", recursive = TRUE)
  
  #here goes the test
  expect_true(TRUE)
})

test_that("test that FitMS() with multiStepMode=constrainedFit with useBootstrap=TRUE and randomSeed, no parallal", {
  #test substitutions
  catalogues <- read.table("FitAndFitMStests/FitandFitMStestCatalogues.tsv",sep = "\t",check.names = F,header = T,stringsAsFactors = F)
  res1 <- FitMS(catalogues = catalogues,
               organ = "Breast",
               multiStepMode = "constrainedFit",
               useBootstrap = TRUE,
               nboot = 5,
               threshold_percent = 0.1,
               nparallel = 1,
               randomSeed = 1)
  res2 <- FitMS(catalogues = catalogues,
                organ = "Breast",
                multiStepMode = "constrainedFit",
                useBootstrap = TRUE,
                nboot = 5,
                threshold_percent = 0.1,
                nparallel = 1,
                randomSeed = 1)
  
  #here goes the test
  expect_equal(res1$exposures,res2$exposures)
})

test_that("test that FitMS() with multiStepMode=constrainedFit with useBootstrap=TRUE and randomSeed, use parallal", {
  #test substitutions
  catalogues <- read.table("FitAndFitMStests/FitandFitMStestCatalogues.tsv",sep = "\t",check.names = F,header = T,stringsAsFactors = F)
  res1 <- FitMS(catalogues = catalogues,
                organ = "Breast",
                multiStepMode = "constrainedFit",
                useBootstrap = TRUE,
                nboot = 5,
                threshold_percent = 0.1,
                nparallel = 2,
                randomSeed = 1)
  res2 <- FitMS(catalogues = catalogues,
                organ = "Breast",
                multiStepMode = "constrainedFit",
                useBootstrap = TRUE,
                nboot = 5,
                threshold_percent = 0.1,
                nparallel = 2,
                randomSeed = 1)
  
  #here goes the test
  expect_equal(res1$exposures,res2$exposures)
})
