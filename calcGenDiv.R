weirsD <- function(afs){
  ## calculate average expected heterozygosity, given a vector of 
  ## derived (or ancestral) frequencies at a set of bi-allelic loci. 
  ## ** actually equivalent to calculating mean( 2 * p * (1-p) ) **
  
  ## number of informative loci
  m <- sum(!is.na(afs))
  ## extract available af values
  yt <- afs[!is.na(afs)]
  return(1 - (sum((yt^2) + (1-yt)^2) / m))
}
