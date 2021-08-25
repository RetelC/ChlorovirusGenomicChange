extractFieldFromFbsVCF <- function(vcf_in, field = "DP", field.numeric = T){
  ## from a .vcf file produced by freebayes without header (!), 
  ## extract the {field} value from the INFO column. 
  ## if {field.numeric}, the values are transformed to numeric before
  ## returning
  
  ## check if the information field is present in every row
  if(sum(grepl(field, x = vcf_in$INFO)) != nrow(vcf_in)){
    stop(paste0("Input field {", field, "} is not present in every row of the input file"))
  }
  ## extract it
  value_out <- strsplit(vcf_in$INFO, split = ";") %>% 
    (function(x1) sapply(x1, (function(info) grep(paste0("^", field, "="), x = info, value = T)))) %>%
    (function(x2) strsplit(x2, split = "=")) %>%
    (function(x3) sapply(x3, "[[", 2)) 
  
  ## for loci with more than one alternative allele, for now simply
  ## take the value corresponding to the first allele (but produce 
  ## a warning)
  if(any(grepl(",", value_out))){
    warning(paste0("The following entries have more than one alternative allele and ", field, " of first allele is returned: ", paste0(which(grepl(",", value_out)), collapse = ", ")))
    value_out <- strsplit(value_out, split = ",") %>%
      (function(x1) sapply(x1, "[[", 1))
  }
  
  
  ## convert to numeric and return
  if(field.numeric) value_out <- as.numeric(value_out)
  return(value_out)
  
}

extractFieldFromVCFColumn <- function(
  vcf_in, sample_name = "1A", field = "DP", field.numeric = T
){
  ## from an input vcf file, extract the {field} value from the
  ## sample given by {sample_name}. This is the equivalent of 
  ## extractFieldFromFbsVCF(), which does the same for the across-sample
  ## information. 
  ## check if {field} is present in format
  if(!(sum(grepl(paste0(field), x = vcf_in$FORMAT)) == nrow(vcf_in))){
    stop(paste0("{ ", field, " } is not present in every row of vcf_in"))
  }
  ## check if {sample_name} is present as a column name
  if(sum(colnames(vcf_in) == sample_name) != 1){
    stop(paste0("{ ", sample_name, " } is not equal to precisely one element of colnames(vcf_in)"))
  }
  ## corner case: if vcf_in is empty, return a length-0 vector
  if(nrow(vcf_in) < 1){
    if(field.numeric) return(numeric()) else return(character())
  }

  ## determine which field of the FORMAT column corresponds to {field}
  idx_fld <- strsplit(vcf_in$FORMAT, split = ":") %>%
    (function(obj) sapply(obj, function(x1) (which(grepl(paste0("^", field, "$"), x = x1)))))
  ## split the concerning column from vcf_in
  sample_split <- strsplit(vcf_in[, sample_name], split = ":")
  
  ## extract the {idx_fld}
  value_out <- character(nrow(vcf_in))
  for(i in 1:nrow(vcf_in)){
    value_out[i] <- sample_split[[i]][idx_fld[i]]
  }
  
  ## for loci with more than one alternative allele, for now simply
  ## take the value corresponding to the first allele (but produce 
  ## a warning)
  if(any(grepl(",", value_out))){
    warning(paste0("The following entries have more than one alternative allele and ", field, " of first allele is returned: ", paste0(which(grepl(",", value_out)), collapse = ", ")))
    value_out <- strsplit(value_out, split = ",") %>%
      (function(x1) sapply(x1, "[[", 1))
  }
  
  ## for loci with more than one alternative allele, for now simply
  ## take the value corresponding to the first allele (but produce 
  ## a warning)
  if(any(grepl(",", value_out))){
    idx_multiallelic <- which(grepl(",", value_out))
    warning(paste0("The following entries have more than one alternative allele and ", field, " of first allele is returned: ", paste0(idx_multiallelic, collapse = ", ")))
    value_out[idx_multiallelic] <- strsplit(
      value_out[idx_multiallelic], split = ","
    ) %>% (function(x1) sapply(x1, "[[", 1))
  }
  ## convert to numeric; in case of missing information, value_out consists
  ## of a '.'. In this case, as.numeric() returns an NA (which is desired) 
  ## as well as a warning (which I suppress). 
  if(field.numeric) suppressWarnings(expr = value_out <- as.numeric(value_out))
  
  ## return
  return(value_out)
}


extractAlternateAFFromVCFColumn <- function(
  vcf_in, sample_name = "1A"
){
  ## from an input vcf file, extract the {AD (Allelic Depth)} value from the
  ## sample given by {sample_name}. Divide the maximum of the elements
  ## that is not the first element by the overall sum of the elements. 
  ## This function solves the problem that AO / RO does not work in 
  ## case of multi-allelic loci. 
  ## check if {AD} is present in format
  if(!(sum(grepl("AD", x = vcf_in$FORMAT)) == nrow(vcf_in))){
    stop(paste0("{ AD } is not present in every row of vcf_in"))
  }
  ## check if {sample_name} is present as a column name
  if(sum(colnames(vcf_in) == sample_name) != 1){
    stop(paste0("{ ", sample_name, " } is not equal to precisely one element of colnames(vcf_in)"))
  }
  ## corner case: if vcf_in is empty, return a length-0 vector
  if(nrow(vcf_in) < 1){
    return(numeric())
  }
  ## corner case: when a position is not covered by any reads, 
  ## freebayes outputs a "."
  
  ## determine which field of the FORMAT column corresponds to {AD}
  idx_ad <- strsplit(vcf_in$FORMAT, split = ":") %>%
    (function(obj) sapply(obj, function(x1) (which(grepl("^AD$", x = x1)))))
  ## split the concerning column from vcf_in
  sample_split <- strsplit(vcf_in[, sample_name], split = ":")
  
  ## extract the {idx_ad}'t column of sample_split
  ad <- character(nrow(vcf_in))
  for(i in 1:nrow(vcf_in)){
    ad[i] <- sample_split[[i]][idx_ad[i]] 
  }
  
  ## now, split by ',', take the maximum of the elements that's not
  ## the first element (the first corresponds to reference depth), 
  ## and divide by the sum of all elements. 
  ## ! In case of missing information, freebayes outputs a ".", and 
  ## { max() / sum() } gives a warning before returning NA. 
  ## I suppressed the warning in this line code because the behavious is 
  ## as intended. !
  aaf_out <- strsplit(ad, split = ",") %>% 
    (function(obj1) lapply(obj1, as.numeric)) %>%
    (function(obj2) suppressWarnings(
      expr = sapply(obj2, (function(x2) max(x2[-1]) / sum(x2)))
    ))
  
  ## output a warning in case of multi-allelic sites: 
  if(any(grepl(",[0-9]+,", x = ad, perl = T))){
    warning(paste0("The following entries have more than one alternative allele and the frequency of the most common allele (calculated per sample) is returned: ", paste0(which(grepl(",[0-9]+,", x = ad, perl = T)), collapse = ", ")))
  }
  
  ## return
  return(aaf_out)
}

calcProbOneBinom <- function(AO = 0, RO = 0){
  ## For a set of RO:AO values, calculate the probability that it 
  ## is a draw from a binomial distribution with a probability equal to the
  ## overall allele frequency 
  
  ## check if input is numeric, and both vectors have the same length
  if((!is.numeric(AO)) | (!is.numeric(RO))){
    stop("AO and RO are not both numeric")
  }
  if(length(AO) != length(RO)){
    stop("AO and RO do not have the same length")
  }
  
  ## corner case: 
  if(length(AO) == 0) return(numeric(0))
  
  ## calculate overall allele frequency
  AFoverall <- sum(AO) / sum(AO + RO)
  
  ## for all observation pairs, calculate the probability that they
  ## come from a binomial distribution with { p == AFoverall }. 
  out <- pmin(
    pbinom(q = AO, prob = AFoverall, size = AO + RO), 
    pbinom(q = AO, prob = AFoverall, size = AO + RO, lower.tail = F)
  )
  return(out)
}
calcProbOneBinom <- function(AO = 0, RO = 0){
  ## For a set of RO:AO values, calculate the probability that they are
  ## draws from a binomial distribution with a probability equal to the
  ## overall allele frequency 
  
  ## check if input is numeric, and both vectors have the same length
  if((!is.numeric(AO)) | (!is.numeric(RO))){
    stop("AO and RO are not both numeric")
  }
  if(length(AO) != length(RO)){
    stop("AO and RO do not have the same length")
  }
  
  ## corner case: 
  if(length(AO) == 0) return(numeric(0))
  
  ## calculate overall allele frequency
  AFoverall <- sum(AO) / sum(AO + RO)
  
  ## for all observation pairs, calculate the probability that they
  ## come from a binomial distribution with { p == AFoverall }. 
  out <- pmin(
    pbinom(q = AO, prob = AFoverall, size = AO + RO), 
    pbinom(q = AO, prob = AFoverall, size = AO + RO, lower.tail = F)
  )
  return(out)
}
