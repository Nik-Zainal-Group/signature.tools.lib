
# micro-homology between the deleted sequence and 3' sequence (number of bases)
mhcaller <- function(d, prime3) {
                                        # d is the change of a deletion
                                        # prime 3 is the sequence 3' of the deleted segment
  countmh = 0
  seq = '-'
  # Dealing with microhomology or lack of microhmomology in the first position

  
  for (i in 1:(nchar(d))) { # for all substrings
     
      if (substr(d, 1, i)==substr(prime3,1,i)) {
          countmh <- countmh + 1
          seq <- substr(d, 1, i)
      }
  }
  result <- list()
  result$countmh <- countmh # number of basepairs matching between the deleted fragment and the 3' segment
  result$seq <- seq
  result
}


# counts how many times the pattern is repeated within the string
tandemcount <- function(pat,string) {
  sum = 0

  start.pos <- seq(from=1, to=nchar(string)-nchar(pat)+1, by=nchar(pat))

  for (s in start.pos) {
      if (substr(string, s, s+nchar(pat)-1) == pat) {
          sum <- sum + 1
      }
  }
  sum
  
}

# get all factors of an integer
get_all_factors <- function(n)
{
  prime_factor_tables <- lapply(
    setNames(n, n), 
    function(i)
    {
      if(i == 1) return(data.frame(x = 1L, freq = 1L))
      plyr::count(as.integer(gmp::factorize(i)))
    }
  )
  lapply(
    prime_factor_tables, 
    function(pft)
    {
      powers <- plyr::alply(pft, 1, function(row) row$x ^ seq.int(0L, row$freq))
      power_grid <- do.call(expand.grid, powers)
      sort(unique(apply(power_grid, 1, prod)))
    }
  )
}


# ----------------------------------------------------------
# This finds the smallest repeating subunit of the deletion
                                        # ----------------------------------------------------------
findsmallestrep <- function(d) {
    d.factors <- as.numeric(sort(unlist(get_all_factors(nchar(d))), decreasing=TRUE))
    rep.unit <- ''
    
    for (f in d.factors) {
        no.repeats <- nchar(d)/f
        rep.string <- paste0(rep(substring(d,1,f), no.repeats), collapse='')
        if (d==rep.string) {
            rep.unit <- substring(d,1,f)
        }
    }
    rep.unit
}


# ------------------------------------------------ 
# This is the Repeat caller
# ------------------------------------------------      

# Taken in the deletion sequence, the 3' context and the length of the deletion
# I have joined deletion and 3' context. To remove this, the rules for microhom/rep calling need to change
# 
repcaller <- function(d, prime3, prime5, l) {
                                        # d : deletion
                                        # prime3 : 3' context
                                        # prime5: 5' context
                                        # l : length of change

    result <- list() 
    result$countrep = 0
    result$rep<-''
                                        # This is for counting single base repeats
                                        # if the length of change is 1
    if (l==1) {
        countrep <- 0
        i <- 1
        while (substr(d,1,1)==substr(prime3,i,i)) {
            countrep <- countrep +1
            i <- i + 1
        }
        result$countrep <- countrep
        result$rep <- d
        return(result)
    } else if (d==substr(prime3, 1, nchar(d))) { # This is for counting whole deletion/DI repeats that are present in the flanking region  
       
        countrep = tandemcount(d,prime3)
        rep = findsmallestrep(d)
        countrep = max(countrep,tandemcount(rep,prime3))
        result$countrep <- countrep
        result$rep <- rep
        return(result)
    } else {   # This is for counting anything in between 1bp and whole del repetition                                     
        rep = '-'
        
        for (t in seq(from=(nchar(d)-1), to=2)) {  # Look for repeats of 2bp to n-1 bp of the length of the indel
            if (grepl(substring(d,1,t), prime3, perl=TRUE)) {
                  countrep = tandemcount(substr(d, 1, t),prime3)
                  rep = findsmallestrep(substr(d, 1, t))
                  unit = tandemcount(rep,d)*nchar(rep)
                                        # The false calls arise in examples such as these : del = AACCCCATCTCTACT; 3' = AAAATTACAAACAAAT; rep = 'A'; repcount in 3' = 4 which is greater than MH = 2; Therefore, it is called Repeat-mediated
                                        # In fact, it should be repeat count = 0; Therefore, call should be microhomology mediated. 
                                        # To do this, compare check how far the repeat stretched into the indel itself. Eg: 'A' is counted twice in the deletion. So compare del[:2] to del[2:4]. If they are the same,then keep it, else false
                  if (substr(d,1,unit) == substr(d, unit+1,unit*2)) {
                                        #print countrep, rep, tandemcount(rep,d), tandemcount(rep,prime3), "Repeat", unit, d[:unit], d[unit:unit*2]
                      countrep = max(countrep, tandemcount(rep,prime3))         
                      result$countrep <- countrep
                      result$rep <- rep
                      return(result)
                  } else {
                      result$countrep <- 0
                      result$rep <- '-'
                      return(result)
                  }
              } # endif
        } # end for
        result$countrep <- 0 # in case no repeat 2bp or longer is fount
        result$rep <- '-'
        return(result)
    } # end else
} # end repcaller



                                        # Given microhomology is in bases, repeat should be in bases too
                                        # If there is a single repeat of the indel 3' of it, then it should be labelled as Repeat-mediated not MH. Eg : TTTA	TTTATTATTAAGATTTTTAAATTTTAATT has 4bp MH and 1 repeat of TTTA. Counting itself, this is a repeat of 2, so it is repeat-mediated.15.05.14
                                        # Except if it is a single base from a longer indel that is repeating, then it is treated as MH
  
finalcaller <- function(mhcount, replength, rept) {
                                        # mhcount
                                        # replength
                                        # rept

    
    if (replength >= mhcount) {
        if ((replength/nchar(rept)) >= 1) {
            return("Repeat-mediated")
        } else if (mhcount > 0) {
            return("Microhomology-mediated")
        } else {
            return("None")
        }
    } else {
        if (mhcount > 0) {
            return("Microhomology-mediated")
        } else { 
            return("None")
        }
    }

}

mh <- function(indel.df) {

                                        # indel.df needs following columns:
                                        # indel.type
                                        # change
                                        # slice3
                                        # slice5
                                        # indel.length
    
    
    classification <- rep (NA, nrow(indel.df))

    if (nrow(indel.df)>0) {
        for (i in 1:nrow(indel.df)) {

            if (indel.df$indel.type[i]=='D') { # the classification is only for deletions
                
                as = as.character(indel.df$change[i]) # The actual deletion
                bs = as.character(indel.df$slice3[i])            
                cs = as.character(indel.df$slice5[i]) # Sequence 5' to deletion            
                
                mhcount = mhcaller(as,bs)$countmh  #  number of microhomology bases            
                
                                        #Look for Microhomology first and then for repeats - tandem/normal
                r = repcaller(as,bs,cs,indel.df$indel.length[i])
                repcount <- r$countrep # number of times there is a repeat
                rept <- r$rep # the repeat sequence
                
                finalcall = finalcaller(mhcount, repcount*nchar(rept), rept)
                classification[i] <- finalcall
                                        #print a, b, mhcount, mh, repcount*len(repeat), repeat, finalcall
            } # end if
        } # end for
        
        classification[is.na(classification)] <- '-'
        
        indel.df$classification <- classification
    } else {
    }
    
    indel.df
}  # end the main mh function  




    
