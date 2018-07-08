toPyr <- function(nctds) {
    nctds.pyr <- nctds
    nctds.pyr[nctds=='A'] <- 'T'
    nctds.pyr[nctds=='T'] <- 'A'
    nctds.pyr[nctds=='G'] <- 'C'
    nctds.pyr[nctds=='C'] <- 'G'
    as.character(nctds.pyr)
}
