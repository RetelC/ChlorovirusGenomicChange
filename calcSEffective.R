calcSEffective <- function(afs, tps, includezeroes=TRUE){
  ## from input numeric vectors giving observed allele frequencies and
  ## time points (in generations), calculate effective selection coefficient
  ## following (Kosheleva et al., MBE 2017)
  ## output is a numeric vector of length {length(afs) - 1)} 
  ## includezeroes: when FALSE, output NA if allele frequency on either
  ##   side is 0 or 1
  ##   when TRUE, change 0 to 1/1000 and 1 to 1-(1/1000), proceed as above.
  ##   (these will be underestimations)
  
  ## check input
  if(!is.numeric(afs) | !is.numeric(tps)){
    stop("calcSeff: both afs and tps should be numeric")
  }
  if(length(afs) != length(tps)){
    stop("calcSeff: afs and tps don't have the same length")
  }
  if(length(afs) < 2 ){
    stop("calcSeff: length of afs needs to be larger than 1")
  }
  
  ## deal with including observed 0's or not
  if(includezeroes){
    afs[afs==0] <- 1/1000; afs[afs==1] <- 1-(1/1000)
  }
  
  ## number of time points
  nt <- length(afs)
  ## calculate length of time intervals
  delta_t <- tps[2:nt] - tps[1:(nt-1)]
  ## 
  af_1 <- afs[1:(nt-1)]; af_2 <- afs[2:nt]
  out <- (1 / delta_t) * log(
    (af_2 * (1 - af_1)) / (af_1 * (1 - af_2))
  )
  
  ## If locus has no variation at both time points, 
  ## set output to NA. The else-code can be improved. 
  if(includezeroes){
    out[(af_1<=1/1000 & af_2<=1/1000) | 
          (af_1>=1-(1/1000) & af_2>=1-(1/1000))] <- NA
  }else{
    out[(af_1==0 | af_2==0) | (af_1==1 | af_2==1)] <- NA
  }
  
  ## make the time intervals recognizable
  names(out) <- paste0("Gen", tps[1:(nt-1)], "-", tps[2:nt])
  return(out)
}