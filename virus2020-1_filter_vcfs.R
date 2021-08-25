#!/bin/R
#######################################
## Process freebayes pooled genotyping output
## Cas Retel, cas.retel@eawag.ch
## 2020.01.15
#######################################
## After running freebayes on pooled sequencing samples
## and filtering only the RO (Reference Observations) and
## AO (Alternate Observations) column. 
#######################################
rm(list=ls())
require(tidyverse)
pkgs <- list("RColorBrewer", "circlize")
sapply(pkgs, require, character.only=T)

source('~/Documents/HVInt/scripts/processFbsVCF.R')
source('~/Documents/HVInt/scripts/plotAfs.R')

## some directories
dir_proj3 <- '~/Documents/HVInt/Project3_chlorovirusmuts/'
setwd(dir_proj3)
dir_v20 <- '~/Documents/HVInt/Project3_chlorovirusmuts/genom/Virus2020-1/'

## some colours
source('~/Documents/Functions/gg_multiplot.R')
oldmai <- par("mai")
col_h <- "#009E73"
col_v <- "#E69F00"

# col_trts <- rep(c("#00A2FF", "#EF5FA7"), each = 3)
# col_vd <- "#9E4A00"
# col_v2 <- alpha(col_v, .8)
col_ann <- c("#ffffff", "#bfbfbf", "#7b7d7b", "#000000")

## miscellaneous; sample identifiers etc. 
pop.name <- readLines(paste0(dir_v20, "Virus2020-1.repls"))
n.pops <- length(pop.name)


treatment <- str_sub(pop.name, 1, 3) %>% 
  as_factor() %>% fct_relevel(levels = c("SEL", "DMS"))
repl.name <- str_sub(pop.name, 1, 5) %>%
  as_factor %>% 
  fct_relevel(levels = c("SEL_1", "SEL_2", "SEL_3", "DMS_1", "DMS_2", "DMS_3"))
repl.num <- str_replace(
  pop.name, pattern = "^\\w+_(\\d)_\\d+$", replacement = "\\1"
) %>% as.numeric()
repl.day <- str_replace(
  pop.name, pattern = "^\\w+_\\d_(\\d+)$", replacement = "\\1"
) %>% as.numeric()

col_trts <- c("#00A2FF", "#EF5FA7")

## read in freebayes vcf
v20_fbsfilename <- str_c(dir_v20, "fbs_dir/virus2020-1_fbsvanilla.vcf")
v20_fbsfilename_qual20 <- str_c(dir_v20, "fbs_dir/virus2020-1_fbs_qual20.vcf")

v20_fbsheader <- readLines(v20_fbsfilename_qual20) %>% 
  str_subset(pattern = "^##")
v20_fbsvcf_qual20 <- read_delim(
  v20_fbsfilename_qual20, comment = "#", col_names = F, delim = "\t"
)
colnames(v20_fbsvcf_qual20) <- readLines(v20_fbsfilename_qual20) %>% 
  str_subset(pattern = "^#CHROM")  %>%
  str_split(pattern = "\t", simplify = T) %>%
  str_replace(pattern = "#", replacement = "")


## Let's check if everything is okay. 
all(pop.name %in% names(v20_fbsvcf_qual20)[-(1:9)])
all(names(v20_fbsvcf_qual20)[-(1:9)] %in% pop.name)

## start by extracting the DP and TYPE info from the vcf
## (only variants with a quality of 20)
v20_fbsdaf_qual20 <- v20_fbsvcf_qual20 %>%
  select(c(1, 2, 4, 5, 6)) %>% mutate(
    CHROMPOS = str_c(CHROM, POS, sep = "::"), 
    DEPTH = extractFieldFromFbsVCF(v20_fbsvcf_qual20, field = "DP"),
    TYPE = extractFieldFromFbsVCF(v20_fbsvcf_qual20, field = "TYPE", field.numeric = F) %>%
      as_factor %>% 
      fct_relevel(levels = c("snp", "mnp", "ins", "del", "complex"))
  )

## the warning messages occur because some loci have more than one
## alternative state. I'm dealing with this below, ignore warnings for now. 

## focus on SNPs
v20_fbsvcf_qual20snp <- v20_fbsvcf_qual20 %>% 
  filter(v20_fbsdaf_qual20$TYPE == "snp")


## pop.name is not in the same order as in the .vcf. 
## fix that and reorder the factor levels
pop.fbsord <- colnames(v20_fbsvcf_qual20snp)[-(1:9)]
repl.fbsord <- str_sub(colnames(v20_fbsvcf_qual20snp)[-(1:9)], 1, 5) %>%
  as_factor %>% fct_relevel(c("SEL_1", "SEL_2", "SEL_3", "DMS_1", "DMS_2", "DMS_3"))
treatment.fbsord <- str_sub(colnames(v20_fbsvcf_qual20snp), 1, 3) %>%
  as_factor %>% fct_relevel(c("SEL", "DMS"))
day.fbsord <- str_replace(
  colnames(v20_fbsvcf_qual20snp)[-(1:9)], pattern = "\\w+_\\d_(\\d+)", 
  replace = "\\1"
) %>% as.numeric


###################
## extract reference and alternative observations
str_split(v20_fbsvcf_qual20snp$INFO, pattern = ";") %>%
  (function(x1) sapply(x1, (function(info) grep("^TYPE=", x = info, value = T))))
## because most loci are multi-allelic, I need to take a detour 
## and figure out which is the dominant allele. 

## generate objects to store Reference and Alternative observations
## (for AA1 and)
RO <- AO1 <- AO2 <- matrix(
  numeric(n.pops * nrow(v20_fbsvcf_qual20snp)), 
  nrow = nrow(v20_fbsvcf_qual20snp)
)
colnames(RO) <- colnames(AO1) <- colnames(AO2) <- 
  colnames(v20_fbsvcf_qual20snp)[-(1:9)]

## figure out which element of the INFO field contains AO and RO
idx_ro <- strsplit(v20_fbsvcf_qual20snp$FORMAT, split = ":") %>%
  (function(obj) sapply(obj, function(x1) (which(grepl("^RO$", x = x1))))) %>%
  (function(x1) head(x1, 1))
idx_ao <- strsplit(v20_fbsvcf_qual20snp$FORMAT, split = ":") %>%
  (function(obj) sapply(obj, function(x1) (which(grepl("^AO$", x = x1))))) %>%
  (function(x1) head(x1, 1))
## ! this only functions correctly if all loci have identical
## structure, i.e. if the code below returns TRUE: 
# strsplit(v20_fbsvcf_qual20snp$FORMAT, split = ":") %>%
#   (function(obj) sapply(obj, function(x1) (which(grepl("^RO$", x = x1))))) %>%
#   (function(x1) (length(unique(x1)) == 1))

## extract RO and AO's from the vcf: 
for(ipop in 1:n.pops){
  sample_split <- pull(v20_fbsvcf_qual20snp, 9 + ipop) %>% 
    str_split(pattern = ":")
  for(ir in 1:nrow(v20_fbsvcf_qual20snp)){
    if(grepl(",", sample_split[[ir]][idx_ao])){
      ## if the locus has multiple alternative states, store them both: 
      AO1[ir, ipop] <- sample_split[[ir]][idx_ao] %>%
        (function(x1) strsplit(x1, split = ",")) %>% 
        (function(obj1) sapply(obj1, "[[", 1)) %>% as.numeric
      AO2[ir, ipop] <- sample_split[[ir]][idx_ao] %>%
        (function(x1) strsplit(x1, split = ",")) %>% 
        (function(obj1) sapply(obj1, "[[", 2)) %>% as.numeric
    }else{
      ## set AO2 to NA
      AO1[ir, ipop] <- sample_split[[ir]][idx_ao] %>% as.numeric
      AO2[ir, ipop] <- NA
    }
    
    RO[ir, ipop] <- sample_split[[ir]][idx_ro] %>% as.numeric
  }
}

## Now calculate allele frequency
AF1 <- AO1 / (AO1 + AO2 + RO)
AF2 <- AO2 / (AO1 + AO2 + RO)
colnames(AF1) <- colnames(AF2) <- colnames(v20_fbsvcf_qual20snp)[-(1:9)]

## how many loci are multi-allelic? 
plot(AF1 ~ AF2, pch = 16, col = alpha(1, .66))
## This seems pretty continuous (which is bad)

## How much loci actually exist where only AO1 or AO2 is > 0.05? 
idx_notma_pop <- apply(((AF1 > 0.05) & (AF2 > 0.05)), 1, (function(x1) !any(x1, na.rm = T)))

sum(idx_notma_pop)
## So out of the 82 SNPs with a quality of > 20, 61 are never 
## multi-allelic within a population. That's good news. 

## Instead of assessing this per population, I might be a little more 
## stringent and figure out if there is one dominant allele per replicate. 
idx_notma_repl <- matrix(
  logical(nrow(v20_fbsvcf_qual20snp) * nlevels(repl.fbsord)), ncol = nlevels(repl.fbsord)
)
colnames(idx_notma_repl) <- levels(repl.fbsord)

for(irepl in 1:nlevels(repl.fbsord)){
  idx_repl <- which(repl.fbsord == levels(repl.fbsord)[irepl])
  idx_notma_repl[, irepl] <- !(
    apply(AF1[, idx_repl], 1, (function(x1) any(x1 > 0.05, na.rm = TRUE))) & 
      apply(AF2[, idx_repl], 1, (function(x1) any(x1 > 0.05, na.rm = TRUE)))
  )
  
  rm(idx_repl)
}

## Now I need to decide to remove them from all replicates or from one. 
sum(idx_notma_repl, na.rm = TRUE); mean(idx_notma_repl, na.rm = TRUE); 
rowMeans(idx_notma_repl, na.rm = TRUE) %>% table
## Note that 1's are 'good' loci, i.e. never multi-allelic. 
## This tells us that when a locus is multi-allelic in one 
## replicate, is it usually (19 out of 22 times) multi-allelic 
## in at least half of the replicates. 
## Get rid of them! 
idx_notma <- apply(idx_notma_repl, 1, all)
sum(idx_notma)
## This actually leaves 60 loci instead of 61 in idx_notma_pop. 

## So now I only need to make sure that I extract the major alternative allele 
## per replicate. 
## create an object to store major allele frequency values
AOmaj <- AFmaj <- AOmin <- AFmin <- matrix(
  numeric(length(repl.name) * nrow(v20_fbsvcf_qual20snp)), 
  nrow = nrow(v20_fbsvcf_qual20snp)
)
colnames(AOmaj) <- colnames(AFmaj) <- colnames(AOmin) <- colnames(AFmin) <- 
  colnames(v20_fbsvcf_qual20snp)[-(1:9)]

## create a logical object (n.loci * n.replicates) to store if the first
## allele is major or not (meaning that it has the highest average frequency)
idx_firstallelemajor <- matrix(
  logical(nrow(v20_fbsvcf_qual20snp) * nlevels(repl.fbsord)), ncol = nlevels(repl.fbsord)
)
colnames(idx_firstallelemajor) <- levels(repl.fbsord)

for(irepl in 1:nlevels(repl.fbsord)){
  ## determine which columns correspond to this replicate
  idx_repl <- which(repl.fbsord == levels(repl.fbsord)[irepl])
  ## determine whether or not the first alternative allele is the dominant one; 
  ## this can be done for all loci simultaneously
  idx_firstallelemajor[, irepl] <- 
    (rowMeans(AF1[, idx_repl], na.rm = T) > rowMeans(AF2[, idx_repl], na.rm = T))
  ## set NA values to TRUE; it's fine because corresponding AF values
  ## are NA as well
  idx_firstallelemajor[is.na(idx_firstallelemajor)] <- TRUE
  
  ## store the major allele in AFmaj
  ## I don't see a way to prevent the for-loop without creating
  ## additional objects containing the combined information of AF1 and AF2
  for(ir in 1:nrow(v20_fbsvcf_qual20snp)){
    if(idx_firstallelemajor[ir, irepl]){
      AOmaj[ir, idx_repl] <- AO1[ir, idx_repl]
      AFmaj[ir, idx_repl] <- AF1[ir, idx_repl]
      AOmin[ir, idx_repl] <- AO2[ir, idx_repl]
      AFmin[ir, idx_repl] <- AF2[ir, idx_repl]
    }else{
      AOmaj[ir, idx_repl] <- AO2[ir, idx_repl]
      AFmaj[ir, idx_repl] <- AF2[ir, idx_repl]
      AOmin[ir, idx_repl] <- AO1[ir, idx_repl]
      AFmin[ir, idx_repl] <- AF1[ir, idx_repl]
    }
  }
  
  
  rm(idx_repl)
}

## finally, set observations based on less than ten reads to NA
AFmaj[(AOmaj + RO <= 10)] <- NA


# ## Now I need to check if this code does what it is intended to
# idx_firstallelemajor
# ## let's look at locus 32, 46 and 64 
# rbind(colnames(idx_firstallelemajor[32, ]), idx_firstallelemajor[32, ])
# cbind(AF1[32, ], AF2[32, ], AFmaj[32, ]); all(AF2[32, ] == AFmaj[32, ])
# rbind(colnames(idx_firstallelemajor[46, ]), idx_firstallelemajor[46, ])
# all(AF1[46, ] == AFmaj[46, ])
# 
# rbind(colnames(idx_firstallelemajor[64, ]), idx_firstallelemajor[64, ])
# cbind(AF1[64, ], AF2[64, ], AFmaj[64, ]); 
# ## all good. 

## how much loci have different major alleles in different replicates?
idx_qual20snp_multiallelic <- !(rowSums(idx_firstallelemajor) %in% c(0, 6))
## how does this compare to the number of multi-allelic snps? 
grepl(",", x = v20_fbsvcf_qual20snp$ALT)
pull(v20_fbsvcf_qual20snp, ALT) %>% str_detect(pattern = ",")

## create an object of these frequencies and the misc. information
## per locus
v20_fbsdaf_qual20snp <- cbind(
    v20_fbsvcf_qual20snp[, c(1, 2, 4, 5, 6)], 
    CHROMPOS = with(v20_fbsvcf_qual20snp, paste(CHROM, POS, sep = "::")), 
    DEPTH = extractFieldFromFbsVCF(v20_fbsvcf_qual20snp, field = "DP"), 
    AFmaj
)


###################
## NA values per locus
apply(
  v20_fbsdaf_qual20snp[, -(1:7)], 
  1, (function(x1) sum(is.na(x1)))
) %>% hist(
  breaks = 0:35, main = "Number of NA values per locus;\n; freebayes Q20 SNPs, nsb & nap"
)
## Given the high overall coverage and the distribution of NA's per 
## locus, we remove everything with more than three NA values 
## (which is equivalent to more than 10% missing information)

idx_qual20snp_na3 <- apply(
  v20_fbsdaf_qual20snp[, -(1:7)], 
  1, (function(x1) sum(is.na(x1)) <= 3)
)

## number of SNPs removed, counting the multi-allelic ones double: 
sum(!idx_qual20snp_na3) + sum((!idx_qual20snp_na3) & idx_qual20snp_multiallelic)
## NA values per locus
###################

###################
## Filter loci with allele frequencies that are all binomial samples
## from one distribution; PER REPLICATE
## These are likely artefacts
## first create objects to store the probability per population
## per frequency estimate (RO-AO combination)

## PSD stands for "probability of same distribution"
## create an object to store the results; dimensions {n.loci * n.replicate},
## and an object to store the 'notsinglebinomialperrepl' logical value

PSD <- PSDcorr <- matrix(
  rep(TRUE, length(repl.name) * nrow(v20_fbsvcf_qual20snp)),
  nrow = nrow(v20_fbsvcf_qual20snp)
)
colnames(PSD) <- colnames(PSDcorr) <- colnames(v20_fbsvcf_qual20snp)[-(1:9)]

idx_qual20snp_nsbperrepl <- matrix(
  rep(TRUE, nrow(v20_fbsvcf_qual20snp) * length(levels(repl.name))),
  nrow = nrow(v20_fbsvcf_qual20snp)
)
colnames(idx_qual20snp_nsbperrepl) <- levels(repl.name)

## For every RO:AO combination, calculate the probability that it
## is a draw from a binomial distribution with a probability equal to the
## overall allele frequency. I am doing this per replicate:

for(irepl in 1:nlevels(repl.name)){
  idx_repl <- (repl.fbsord == levels(repl.name)[irepl])
  for(iloc in 1:nrow(v20_fbsvcf_qual20snp)){
    ## for every RO:AO combo, calculate its probability assuming a
    ## binomial distribution with overall.AF
    PSD[iloc, idx_repl] <- calcProbOneBinom(
      AO = AOmaj[iloc, idx_repl],
      RO = AOmin[iloc, idx_repl] + RO[iloc, idx_repl]
    )
    
    ## If any allele frequencies > 0.05 exist and the RO:AO's all come
    ## from a single binomial distribution, set nsbperrepl to FALSE
    if(
      any(AFmaj[iloc, idx_repl] > 0.05, na.rm = TRUE) &
      all(p.adjust(PSD[iloc, idx_repl], method = "BH") > 0.05, na.rm = TRUE)
    ){ idx_qual20snp_nsbperrepl[iloc, irepl] <- FALSE }
  }
  
  rm(idx_repl)
}


## If a row of idx_qual20snp_nsbperrepl contains at least one FALSE,
## we remove the locus from further analysis.
idx_qual20snp_nsb_calcperrepl <- apply(idx_qual20snp_nsbperrepl, 1, all)
sum(idx_qual20snp_nsb_calcperrepl)
sum(idx_qual20snp_nsb_calcperrepl & idx_qual20snp_na3)
## _nsb stands for 'not single binomial'
## 52 left

## number of SNPs removed by _nsb, counting the multi-allelic ones double: 
sum(!idx_qual20snp_nsb_calcperrepl) + sum((!idx_qual20snp_nsb_calcperrepl) & idx_qual20snp_multiallelic)


# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_afsperlocus_fbssnpsQ20_onebinomdis.pdf"), width=6, height=4)
for(i in which(!idx_qual20snp_nsb_calcperrepl)){
  hist(
    as.numeric(v20_fbsdaf_qual20snp[i, -(1:7)]), breaks = .05*(0:20),
    xlab = "Alternative allele frequency", ylab = "Number of samples",
    main = paste0(
      "Frequencies of SNPs at ", v20_fbsdaf_qual20snp[i, "CHROMPOS"],
      "\nfreebayes Q20 SNPs, draws from one binomial dist"
    )
  )
}
# dev.off(); system(paste0("open ", dir_v20, "figures/variantscalled_freebayes/fig_afsperlocus_fbssnpsQ20_onebinomdis.pdf"))
## Filter loci with allele frequencies that are all binomial samples
## from one distribution; PER REPLICATE
###################


###################
## Check out consistently polymorphic loci 
## _nap for "not always polymorphic"

idx_qual20snp_nap <- apply(
  v20_fbsdaf_qual20snp, 1, (function(x1) 
    any(x1[grep("SEL_", colnames(v20_fbsdaf_qual20snp))] <= .01, na.rm = T) & 
      any(x1[grep("DMS_", colnames(v20_fbsdaf_qual20snp))] <= .01, na.rm = T)
  )
)

# idx_qual20snp_nap2 <- apply(
#   AOmaj, 1, (function(x1) 
#     any(x1[grep("SEL_", colnames(v20_fbsdaf_qual20snp))] <= 2, na.rm = T) & 
#       any(x1[grep("DMS_", colnames(v20_fbsdaf_qual20snp))] <= 2, na.rm = T)
#   )
# )
## I am as yet unsure if I should impose a filter based on frequency 
## or based on count. I think counts make more sense if I want to 
## remove sequencing errors. But they're more or less equivalent 
## of course. 

## to what extent do these overlap with _nsb? 
sum(idx_qual20snp_nsb_calcperrepl); sum(idx_qual20snp_nap); 
table((idx_qual20snp_nsb_calcperrepl) + (idx_qual20snp_nap))


## 30 of these overlap. 41 SNPs left. 

# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_afsperlocus_fbssnpsQ20_alwayspolymorphic.pdf"), width=6, height=4)
# for(i in which(!idx_qual20snp_nap)){
#   hist(
#     as.numeric(v20_fbsdaf_qual20snp[i, -(1:7)]), breaks = .05*(0:20), 
#     xlab = "Alternative allele frequency", ylab = "Number of samples", 
#     main = paste0(
#       "Frequencies of SNPs at ", v20_fbsdaf_qual20snp[i, "CHROMPOS"], 
#       "\nfreebayes Q20 SNPs, consistently polymorphic loci"
#     )
#   )
# }
# dev.off(); system(paste0("open ", dir_v20, "figures/variantscalled_freebayes/fig_afsperlocus_fbssnpsQ20_alwayspolymorphic.pdf"))


## plot allele frequencies of alleles removed by _nap but not by
## _nap2; kept code because it can be nicely recycled. 
# tmpdafs <- v20_fbsdaf_qual20snp[
#   which(((idx_qual20snp_nsb_calcperrepl) & (idx_qual20snp_nap)) & 
#           !((idx_qual20snp_nsb_calcperrepl) & (idx_qual20snp_nap2))), 
# ]
# for(repli in levels(repl.name)){
#   plotAfs(
#     af_mat = as.matrix(tmpdafs[, grep(repli, names(tmpdafs))]), 
#     tps = strsplit(
#       grep(repli, names(tmpdafs), value = TRUE), 
#       split = "_"
#     ) %>% (function(obj1) sapply(obj1, "[[", 3)) %>% as.numeric, 
#     # col = c("yellow", "orange", "red", "light blue", "navy"), 
#     title_v = paste0(
#       repli, ": consistently polymorphic based on allele count\nbut not based on allele freq"
#     )
#   )
#   
# }
# rm(tmpdafs)

## combine all three filters
idx_qual20snp_na3_nsb_nap <- apply(idx_qual20snp_nsbperrepl, 1, all) & 
  apply(
    v20_fbsdaf_qual20snp, 1, (function(x1) 
      any(x1[grep("SEL_", colnames(v20_fbsdaf_qual20snp))] <= .01, na.rm = T) & 
        any(x1[grep("DMS_", colnames(v20_fbsdaf_qual20snp))] <= .01, na.rm = T)
    )
  ) & apply(
    v20_fbsdaf_qual20snp[, -(1:7)], 
    1, (function(x1) sum(is.na(x1)) <= 3)
  )


v20_fbsvcf_qual20snp_na3_nsb_nap <- 
  v20_fbsvcf_qual20snp[idx_qual20snp_na3_nsb_nap, ]
v20_fbsdaf_qual20snp_na3_nsb_nap <- v20_fbsdaf_qual20snp[idx_qual20snp_na3_nsb_nap, ]

## number of SNPs removed by _nap, counting the multi-allelic ones double: 
sum(!idx_qual20snp_nap) + sum((!idx_qual20snp_nap) & idx_qual20snp_multiallelic)


## Check out consistently polymorphic loci 
###################

###################
## WRITE FILTERED OUTPUT TO FILE
## These files go into bcftools norm, then into snpeff. 
# outfile_v20filt_vcf <- paste0(dir_v20, "fbs_dir/virus2020-1_fbs_filtered.vcf")
# outfile_v20filt_daf <- paste0(dir_v20, "fbs_dir/virus2020-1_fbs_filtered.daf")
# 
# writeLines(
#   v20_fbsheader, con = outfile_v20filt_vcf
# )
# write_delim(
#   v20_fbsvcf_qual20snp[idx_qual20snp_na3_nsb_nap, ], 
#   file = outfile_v20filt_vcf, append = T
# )
# 
# write_delim(
#   v20_fbsdaf_qual20snp, 
#   file = outfile_v20filt_daf
# )
# 
## WRITE FILTERED OUTPUT TO FILE
###################


###################
## Some statistics on the numbers and distribution of polymorphic sites
## Fixed differences from the ancestor
apply(
  v20_fbsdaf_qual20snp[idx_qual20snp_na3_nsb_nap, -(1:7)], 1, 
  (function(x1) all(x1 > .99))
) %>% sum
## no snps like this exist. 

## the tidyverse variant: 
v20_fbsdaf_qual20snp[idx_qual20snp_na3_nsb_nap, ] %>% 
  rowwise() %>% summarise(
    n.fixed = (all(c_across(SEL_3_51:SEL_1_15) > .99, na.rm = T))
  ) %>% summarise(variants.fixed = sum(as.numeric(n.fixed)))



## Number of polymorphic sites per population
apply(
  v20_fbsdaf_qual20snp[idx_qual20snp_na3_nsb_nap, -(1:7)], 2, 
  (function(x1) sum((x1 > .01) & (x1 <= .99), na.rm = T))
)

## number of polymorphic sites per replicate
for(repi in levels(repl.name)){
  print(cbind(
    repi, 
    apply(
      v20_fbsdaf_qual20snp_na3_nsb_nap[
        , grepl(repi, names(v20_fbsdaf_qual20snp_na3_nsb_nap))
      ], 1, function(x1) (sum((x1 > .01) & (x1 <= .99), na.rm = TRUE) > 1)
    ) %>% sum
  ))
}
## Some statistics on the numbers and distribution of polymorphic sites
###################

###################
## plot site-frequency spectra; for the heck of it
# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_vafs_fbssnpsQ20_na3_nsb_nap.pdf"), width=6, height=4)
# for(popi in pop.name){
#   hist(
#     v20_fbsdaf_qual20snp[[popi]][idx_qual20snp_na3_nsb_nap], 
#     breaks = .05*(0:20), xlab = "Allele frequency", 
#     main = paste0("VAF of ", popi, "; freebayes Q20 SNPs, nsb & nap & NA3")
#   )
# }
# dev.off(); system(paste0("open ", dir_v20, "figures/variantscalled_freebayes/fig_vafs_fbssnpsQ20_na3_nsb_nap.pdf"))
# ## A lot of 0's, as expected. This looks much better than for Chlorella. 

## plot site-frequency spectra; for the heck of it
###################



###################
## Allele frequency trajectories over time. These plots are 
## created per replicate. 
## Let's simply create a colour palette for the 41 SNPs left
set.seed(1212)
col_allloci <- sample(rainbow(nrow(v20_fbsdaf_qual20snp_na3_nsb_nap)))
names(col_allloci) <- v20_fbsdaf_qual20snp_na3_nsb_nap$CHROMPOS

## create list objects to store time points and filtered allele
## frequency values. 
perrep_day <- perrep_fbsdaf_fltv20 <- lapply(levels(repl.name), (function(x1) NULL))
## extract the days, and the allele per replicate that reach at least 0.01
## frequency; don't want to clog the bottom of the figures
for(i in 1:nlevels(repl.name)){
  perrep_day[[i]] <- strsplit(
    grep(levels(repl.name)[i], names(v20_fbsdaf_qual20snp_na3_nsb_nap), value = TRUE), split = "_"
  ) %>% (function(obj1) sapply(obj1, "[[", 3)) %>% as.numeric
  perrep_fbsdaf_fltv20[[i]] <- v20_fbsdaf_qual20snp_na3_nsb_nap[
    apply(
      v20_fbsdaf_qual20snp_na3_nsb_nap[
        , grep(levels(repl.name)[i], names(v20_fbsdaf_qual20snp_na3_nsb_nap))
      ], 1, (function(x1) (any(x1 > 0.01, na.rm = TRUE)))
    ), 
    c(1:7, grep(levels(repl.name)[i], names(v20_fbsdaf_qual20snp_na3_nsb_nap)))
  ]
}

# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_afsovertime_sel_fbssnpsQ20_na3_nsb_nap.pdf"), width=6, height=6)
par(mfrow = c(3, 1))
for(i in 1:3){
  plotAfs(
    af_mat = as.matrix(perrep_fbsdaf_fltv20[[i]][, -(1:7)]), 
    tps = perrep_day[[i]], vline_sampletimes = F, 
    col_v = col_allloci[perrep_fbsdaf_fltv20[[i]]$CHROMPOS], 
    title_v = paste0("PBCV1 ", levels(repl.name)[i], "; fbs Q20 SNPs na3 nsb nap")
  )
  with(perrep_fbsdaf_fltv20[[i]], legend(
    "topleft", legend = CHROMPOS, fill = col_allloci[CHROMPOS], bty = 'n', 
    cex = .67, ncol = 2
  ))
}
# dev.off()
# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_afsovertime_dms_fbssnpsQ20_na3_nsb_nap.pdf"), width=6, height=6)
par(mfrow = c(3, 1))
for(i in 4:6){
  plotAfs(
    af_mat = as.matrix(perrep_fbsdaf_fltv20[[i]][, -(1:7)]), 
    tps = perrep_day[[i]], vline_sampletimes = F, 
    col_v = col_allloci[perrep_fbsdaf_fltv20[[i]]$CHROMPOS], 
    title_v = paste0("PBCV1 ", levels(repl.name)[i], "; fbs Q20 SNPs na3 nsb nap")
  )
  with(perrep_fbsdaf_fltv20[[i]], legend(
    "topleft", legend = CHROMPOS, fill = col_allloci[CHROMPOS], bty = 'n', 
    cex = .67, ncol = 2
  ))
}
# dev.off(); 
# system(paste0("open ", dir_v20, "figures/variantscalled_freebayes/fig_afsovertime_dms_fbssnpsQ20_na3_nsb_nap.pdf"))



v20_fbsvcf_qual20snp[
  v20_fbsdaf_qual20snp$CHROMPOS %in% c("PBCV1_scaffold_1::24378", "PBCV1_scaffold_1::324025"), 
]
v20_fbsdaf_qual20snp[
  v20_fbsdaf_qual20snp$CHROMPOS %in% c("PBCV1_scaffold_1::24378", "PBCV1_scaffold_1::324025"), 
]

v20_fbsvcf_qual20snp[c("10", "70", "123"), ]


## My filtering pipeline is not as stringent as in the SA analysis. 
## Removing variation at day 0 would not solve this. 
## Allele frequency trajectories over time. 
###################

###################
## continued analysis: save the filtered .vcf file, and also
## two subsets with all polymorphic variants at days 51 and 64


## the full filtered dataset: 
# writeLines(v20_fbsheader, con = outfile_v20filt)
# write.table(
#   matrix(gsub("CHROM", "#CHROM", x = names(v20_fbsvcf_qual20snp_na3_nsb_nap)), nrow = 1), 
#   file = outfile_v20filt, append = T, quote = F, sep = "\t", 
#   col.names = F, row.names = F
# )
# 
# write.table(
#   v20_fbsvcf_qual20snp_na3_nsb_nap, file = outfile_v20filt, 
#   append = T, quote = F, sep = "\t", col.names = F, row.names = F
# )

## check which alleles are polymorphic at days 51 and 58/64
idx_pm51 <- apply(
  v20_fbsdaf_qual20snp_na3_nsb_nap[, 7 + which(day.fbsord == 51)], 1, 
  (function(x1) any(x1 > 0.01, na.rm = T))
)

idx_pm61 <- apply(
  v20_fbsdaf_qual20snp_na3_nsb_nap[, 7 + which(day.fbsord > 51)], 1, 
  (function(x1) any(x1 > 0.01, na.rm = T))
)

## the set of SNPs polymorphic at day 51: 
# writeLines(v20_fbsheader, con = outfile_v20filt_d51)
# write.table(
#   matrix(gsub("CHROM", "#CHROM", x = names(v20_fbsvcf_qual20snp_na3_nsb_nap)), nrow = 1), 
#   file = outfile_v20filt_d51, append = T, quote = F, sep = "\t", 
#   col.names = F, row.names = F
# )
# write.table(
#   v20_fbsvcf_qual20snp_na3_nsb_nap[idx_pm51, ], file = outfile_v20filt_d51, 
#   append = T, quote = F, sep = "\t", col.names = F, row.names = F
# )

## the set of SNPs polymorphic at day 64: 
# writeLines(v20_fbsheader, con = outfile_v20filt_d61)
# write.table(
#   matrix(gsub("CHROM", "#CHROM", x = names(v20_fbsvcf_qual20snp_na3_nsb_nap)), nrow = 1), 
#   file = outfile_v20filt_d61, append = T, quote = F, sep = "\t", 
#   col.names = F, row.names = F
# )
# write.table(
#   v20_fbsvcf_qual20snp_na3_nsb_nap[idx_pm61, ], file = outfile_v20filt_d61, 
#   append = T, quote = F, sep = "\t", col.names = F, row.names = F
# )

## What's really more important are the derived allele frequency values. 
# write.table(
#   v20_fbsdaf_qual20snp_na3_nsb_nap, 
#   file = daffile_v20filt, quote = F, sep = "\t", col.names = T, row.names = F
# )
# write.table(
#   v20_fbsdaf_qual20snp_na3_nsb_nap[idx_pm51, ], 
#   file = daffile_v20filt_d51, quote = F, sep = "\t", col.names = T, row.names = F
# )
# write.table(
#   v20_fbsdaf_qual20snp_na3_nsb_nap[idx_pm61, ], 
#   file = daffile_v20filt_d61, quote = F, sep = "\t", col.names = T, row.names = F
# )

## continued analysis: save the filtered .vcf file
###################


