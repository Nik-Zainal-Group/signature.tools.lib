

classifyRearrangementsFromBedpe <- function(sv_bedpe){
  svclass <- c()
  for (i in 1:nrow(sv_bedpe)){
    if(sv_bedpe[i,"chrom1"]!=sv_bedpe[i,"chrom2"]){
      svclass <- c(svclass,"translocation")
    }else if(sv_bedpe[i,"strand1"]!=sv_bedpe[i,"strand2"]){
      svclass <- c(svclass,"inversion")
    }else if(sv_bedpe[i,"strand1"]=="+"){
      svclass <- c(svclass,"deletion")
    }else if(sv_bedpe[i,"strand1"]=="-"){
      svclass <- c(svclass,"tandem-duplication")
    }
  }
  sv_bedpe[,"svclass"] <- svclass
  #return updated df
  sv_bedpe
}