context("testing functions Fit and FitMS (DBS), should take a few seconds")

test_that("test that Fit() works on a random matrix, using useBootstrap=FALSE", {
  #test substitutions
  mut.order <- c("AC>CA","AC>CG","AC>CT","AC>GA","AC>GG","AC>GT","AC>TA","AC>TG","AC>TT",
                 "AT>CA","AT>CC","AT>CG","AT>GA","AT>GC","AT>TA",
                 "CC>AA","CC>AG","CC>AT","CC>GA","CC>GG","CC>GT","CC>TA","CC>TG","CC>TT",
                 "CG>AT","CG>GC","CG>GT","CG>TA","CG>TC","CG>TT",
                 "CT>AA","CT>AC","CT>AG","CT>GA","CT>GC","CT>GG","CT>TA","CT>TC","CT>TG",
                 "GC>AA","GC>AG","GC>AT","GC>CA","GC>CG","GC>TA",
                 "TA>AT","TA>CG","TA>CT","TA>GC","TA>GG","TA>GT",
                 "TC>AA","TC>AG","TC>AT","TC>CA","TC>CG","TC>CT","TC>GA","TC>GG","TC>GT",
                 "TG>AA","TG>AC","TG>AT","TG>CA","TG>CC","TG>CT","TG>GA","TG>GC","TG>GT",
                 "TT>AA","TT>AC","TT>AG","TT>CA","TT>CC","TT>CG","TT>GA","TT>GC","TT>GG")
  n_row <- 78
  n_col <- 15
  rnd_matrix <- round(matrix(runif(n_row*n_col,min = 0,max = 50),nrow = n_row,ncol = n_col))
  colnames(rnd_matrix) <- paste0("C",1:n_col)
  row.names(rnd_matrix) <- mut.order
  res <- Fit(catalogues = rnd_matrix,
             signatures = getReferenceSignatures(typemut = "DNV"),
             useBootstrap = FALSE,
             threshold_percent = 0.1)

  #here goes the test
  expect_true(TRUE)
})

test_that("test that Fit() works on a random matrix with useBootstrap=TRUE, and that random seed works, no parallel", {
  #test substitutions
  mut.order <- c("AC>CA","AC>CG","AC>CT","AC>GA","AC>GG","AC>GT","AC>TA","AC>TG","AC>TT",
                 "AT>CA","AT>CC","AT>CG","AT>GA","AT>GC","AT>TA",
                 "CC>AA","CC>AG","CC>AT","CC>GA","CC>GG","CC>GT","CC>TA","CC>TG","CC>TT",
                 "CG>AT","CG>GC","CG>GT","CG>TA","CG>TC","CG>TT",
                 "CT>AA","CT>AC","CT>AG","CT>GA","CT>GC","CT>GG","CT>TA","CT>TC","CT>TG",
                 "GC>AA","GC>AG","GC>AT","GC>CA","GC>CG","GC>TA",
                 "TA>AT","TA>CG","TA>CT","TA>GC","TA>GG","TA>GT",
                 "TC>AA","TC>AG","TC>AT","TC>CA","TC>CG","TC>CT","TC>GA","TC>GG","TC>GT",
                 "TG>AA","TG>AC","TG>AT","TG>CA","TG>CC","TG>CT","TG>GA","TG>GC","TG>GT",
                 "TT>AA","TT>AC","TT>AG","TT>CA","TT>CC","TT>CG","TT>GA","TT>GC","TT>GG")
  n_row <- 78
  n_col <- 5
  rnd_matrix <- round(matrix(runif(n_row*n_col,min = 0,max = 50),nrow = n_row,ncol = n_col))
  colnames(rnd_matrix) <- paste0("C",1:n_col)
  row.names(rnd_matrix) <- mut.order
  res1 <- Fit(catalogues = rnd_matrix,
              signatures = getReferenceSignatures(typemut = "DNV"),
              useBootstrap = TRUE,
              nboot = 5,
              threshold_percent = 0,
              threshold_p.value = 0.1,
              nparallel = 1,
              randomSeed = 1)
  res2 <- Fit(catalogues = rnd_matrix,
              signatures = getReferenceSignatures(typemut = "DNV"),
              useBootstrap = TRUE,
              nboot = 5,
              threshold_percent = 0,
              threshold_p.value = 0.1,
              nparallel = 1,
              randomSeed = 1)

  #here goes the test
  expect_equal(res1$exposures,res2$exposures)
})

test_that("test that plotFit() works on a random matrix", {
  #test substitutions
  mut.order <- c("AC>CA","AC>CG","AC>CT","AC>GA","AC>GG","AC>GT","AC>TA","AC>TG","AC>TT",
                 "AT>CA","AT>CC","AT>CG","AT>GA","AT>GC","AT>TA",
                 "CC>AA","CC>AG","CC>AT","CC>GA","CC>GG","CC>GT","CC>TA","CC>TG","CC>TT",
                 "CG>AT","CG>GC","CG>GT","CG>TA","CG>TC","CG>TT",
                 "CT>AA","CT>AC","CT>AG","CT>GA","CT>GC","CT>GG","CT>TA","CT>TC","CT>TG",
                 "GC>AA","GC>AG","GC>AT","GC>CA","GC>CG","GC>TA",
                 "TA>AT","TA>CG","TA>CT","TA>GC","TA>GG","TA>GT",
                 "TC>AA","TC>AG","TC>AT","TC>CA","TC>CG","TC>CT","TC>GA","TC>GG","TC>GT",
                 "TG>AA","TG>AC","TG>AT","TG>CA","TG>CC","TG>CT","TG>GA","TG>GC","TG>GT",
                 "TT>AA","TT>AC","TT>AG","TT>CA","TT>CC","TT>CG","TT>GA","TT>GC","TT>GG")
  n_row <- 78
  n_col <- 15
  rnd_matrix <- round(matrix(runif(n_row*n_col,min = 0,max = 50),nrow = n_row,ncol = n_col))
  colnames(rnd_matrix) <- paste0("C",1:n_col)
  row.names(rnd_matrix) <- mut.order
  res <- Fit(catalogues = rnd_matrix,
             signatures = getReferenceSignatures(typemut = "DNV"),
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
  catalogues <- read.table("FitAndFitMStests/FitandFitMStestCatalogues_DBS.tsv",sep = "\t",check.names = F,header = T,stringsAsFactors = F)
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
  catalogues <- read.table("FitAndFitMStests/FitandFitMStestCatalogues_DBS.tsv",sep = "\t",check.names = F,header = T,stringsAsFactors = F)
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


test_that("test that FitMS() with different common and rare signatures tiers", {
  #test substitutions
  catalogues <- read.table("FitAndFitMStests/FitandFitMStestCatalogues_DBS.tsv",sep = "\t",check.names = F,header = T,stringsAsFactors = F)
  res <- FitMS(catalogues = catalogues,
               organ = "Breast",
               commonSignatureTier = "T2",
               rareSignatureTier = "T0")

  res <- FitMS(catalogues = catalogues,
               organ = "Breast",
               commonSignatureTier = "T3",
               rareSignatureTier = "T1")

  res <- FitMS(catalogues = catalogues,
               organ = "Breast",
               commonSignatureTier = "T1",
               rareSignatureTier = "T3")

  res <- FitMS(catalogues = catalogues,
               organ = "Breast",
               commonSignatureTier = "T2",
               rareSignatureTier = "T4")

  #here goes the test
  expect_true(TRUE)
})

test_that("test that FitMS() with different rare signatures selection criteria", {
  #test substitutions
  catalogues <- read.table("FitAndFitMStests/FitandFitMStestCatalogues_DBS.tsv",sep = "\t",check.names = F,header = T,stringsAsFactors = F)
  res <- FitMS(catalogues = catalogues,
               organ = "Breast",
               rareCandidateSelectionCriteria = "MinError")

  #here goes the test
  expect_true(TRUE)
})
