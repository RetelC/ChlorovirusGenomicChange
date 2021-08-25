#!/bin/R
#######################################
## Visualise viral polymorphism 
## Cas Retel, cas.retel@eawag.ch, 2020.02.04
#######################################
## After genotyping pooled sequencing samples (and filtering variant calls; 
## see virus2020-1_filter_vcfs.R). analyse genomic variation: 
## - create circus plots with location of SNPs (or maybe other 
## measures of diversity) for the six different replicates
## - visualise the distribution of SNP annotations
## - analyse repeatability between replicates (at day 61?)
## - Check variation in certain genomic regions over time. 
#######################################

rm(list=ls())
require(tidyverse)
pkgs <- list("RColorBrewer", "circlize")
sapply(pkgs, require, character.only=T)
source('~/Documents/HVInt/Project3_chlorovirusmuts/scripts/functions_vcf.R')
source('~/Documents/HVInt/scripts/plotAfs.R')

## some directories
dir_proj3 <- '~/Documents/HVInt/Project3_chlorovirusmuts/'
dir_v20 <- '~/Documents/HVInt/Project3_chlorovirusmuts/genom/Virus2020-1/'
dir_db <- '~/Dropbox/Talks_Algae_Virus/Project3_chlorovirusmuts/'

## some colours
source('~/Documents/Functions/gg_multiplot.R')
oldmai <- par("mai")
col_h <- "#009E73"
col_v <- "#E69F00"
col_ann <- c("#ffffff", "#bfbfbf", "#7b7d7b", "#000000")
col_trts <- rep(c("#EF5FA7", "#00A2FF"), each = 3)
col_brred <- "#FF2600"

## file names of .vcfs and .dafs 
vcffile.v20filt <- paste0(dir_v20, "fbs_dir/virus2020-1_fbs_filtered.vcf")

## miscellaneous; sample identifiers etc. 
pop.name <- readLines(paste0(dir_v20, "Virus2020-1.repls")) %>%
  (function(x1) strsplit(x1, split = "[.]")) %>%
  (function(x2) sapply(x2, "[[", 1))
n.pops <- length(pop.name)

treatment <- substr(pop.name, 1, 3) %>% 
  factor(, levels=c("SEL", "DMS"))
repl.name <- substr(pop.name, 1, 5) %>% 
  factor(, levels = c("SEL_1", "SEL_2", "SEL_3", "DMS_1", "DMS_2", "DMS_3"))
repl.num <- strsplit(pop.name, split = "_") %>%
  (function(obj1) sapply(obj1, "[[", 2)) %>% as.numeric
repl.day <- strsplit(pop.name, split = "_") %>%
  (function(obj1) sapply(obj1, "[[", 3)) %>% as.numeric

###################
## read in data. 
## Start with vcf file including annotations per mutation. 
## See how snpeff handles the freebayes multi-allelic mutations. 
v20.filename.vcf <- paste0(dir_v20, "sef_dir/virus2020-1_fbs_filtered.split.sef.vcf")
v20.header <- readLines(v20.filename.vcf) %>% 
  (function(x1) grep("^##", x = x1, value = T))
v20.vcf <- read.table(
  v20.filename.vcf, comment.char = "#", header = F
)
names(v20.vcf) <- readLines(v20.filename.vcf) %>% 
  (function(x1) grep("^#CHROM", x = x1, value = T)) %>%
  (function(x2) strsplit(x2, split = "\t")) %>% unlist %>%
  (function(x3) gsub("#CHROM", "CHROM", x = x3))

## extract annotations
v20.annots <- readVCFAnn(v20.filename.vcf)
## Now I realised that there are mutations in A144L and a147L, two ORF's 
## which overlap with A140/145R and are actually not transcribed. 
## For these four mutations, I'm manually going to check the mut.class 
## etc for the SNP in A140/145R. 
idx.140ol <- str_detect(pull(v20.annots, GENE.NAME), pattern = "(a144L)|(a147L)")

vcf.tmp <- read_delim(v20.filename.vcf, delim = "\t", col_names = F, comment = "#")
names(vcf.tmp) <- c(
  "CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
  "FILTER", "INFO", "FORMAT", str_c("SAMP.", 1:(ncol(vcf.tmp) - 9))
)
## check how much fields the INFO column has per annotation (per gene): 
fields.perannot <- pull(vcf.tmp, "INFO")[1] %>%
  (function(x1) str_split(x1, pattern = ";")[[1]][str_count(x1, pattern = ";") + 1]) %>%
  str_split(pattern = "\\|") %>% unlist %>%
  str_which(pattern = ",[ACTG]") %>% 
  (function(x2) x2[1] - 1)

## create an object to store the second annotation for idx.140ol
anns.second <- lapply(which(idx.140ol), function(x) NULL)
names(anns.second) <- which(idx.140ol)
for(il in 1:length(anns.second)){
  ## extract the second annotation
  anns.second[[il]] <- pull(vcf.tmp, "INFO")[which(idx.140ol)[il]] %>%
    (function(x1) str_split(x1, pattern = ";")[[1]][str_count(x1, pattern = ";") + 1]) %>%
    str_split(pattern = "\\|") %>% unlist %>%
    (function(x2) x2[(fields.perannot + 1):(2 * fields.perannot)])
  
  ## change it in v20.annots
  v20.annots[which(idx.140ol)[il], c("MUT.CLASS", "MUT.SEVERITY", "GENE.NAME")] <-
    anns.second[[il]][c(2, 3, 4)]
}

## generate some describing factors in the order freebayes put the 
## samples in
pop.name.fbsord <- colnames(v20.vcf)[-(1:9)]
treatment.fbsord <- substr(pop.name.fbsord, 1, 3) %>% 
  as_factor() %>% fct_relevel(c("SEL", "DMS"))
repl.name.fbsord <- substr(pop.name.fbsord, 1, 5) %>% 
  as_factor() %>% 
  fct_relevel(c("SEL_1", "SEL_2", "SEL_3", "DMS_1", "DMS_2", "DMS_3"))

repl.num.fbsord <- str_split(pop.name.fbsord, pattern = "_") %>% 
  (function(obj1) sapply(obj1, "[[", 2)) %>% as.numeric

repl.day.fbsord <- str_split(pop.name.fbsord, pattern = "_") %>%
  (function(obj1) sapply(obj1, "[[", 3)) %>% as.numeric

n.pops <- length(pop.name.fbsord)
## read in data. 
###################

###################
## calculate derived allele frequencies from the AO and RO fields in 
## the .vcf. 
RO <- AO <- matrix(
  numeric(n.pops * nrow(v20.vcf)), 
  nrow = nrow(v20.vcf)
)
colnames(RO) <- colnames(AO) <- 
  colnames(v20.vcf)[-(1:9)]

## figure out which element of the INFO field contains AO and RO
idx_ro <- str_split(v20.vcf$FORMAT, pattern = ":") %>%
  (function(obj) sapply(obj, function(x1) (which(grepl("^RO$", x = x1))))) %>%
  unique
idx_ao <- strsplit(v20.vcf$FORMAT, split = ":") %>%
  (function(obj) sapply(obj, function(x1) (which(grepl("^AO$", x = x1))))) %>%
  unique

## extract RO and AO's from the vcf: 
for(ipop in 1:n.pops){
  sample_split <- strsplit(v20.vcf[, 9 + ipop], split = ":")
  for(ir in 1:nrow(v20.vcf)){
    ## set AO2 to NA
    RO[ir, ipop] <- sample_split[[ir]][idx_ro] %>% as.numeric
    AO[ir, ipop] <- sample_split[[ir]][idx_ao] %>% as.numeric
  }
}


## Now calculate allele frequency
v20.dafs <- cbind(
  v20.annots[, c("CHROM", "POS", "MUT.SEVERITY")], 
  (AO / (AO + RO))
) %>% as_tibble

## set observations based on less than 10 measurements to NA
v20.dafs[, -(1:3)][((AO + RO) <= 10)] <- NA

## calculate derived allele frequencies
###################


###################
## find the loci that are polymorphic at days 51 and 61 (:= 58 or 64)
## Define polymorphic as: at least three reads supporting the alternative
## allele. 
AO <- as_tibble(AO)
idx.d61 <- AO %>% select(contains("_58") | contains("_64")) %>% 
  rowwise %>% 
  summarise(idx_d61 = any(c_across(everything()) > 2), .groups = "rowwise") %>%
  pull(idx_d61)
  

idx.d51 <- AO %>% select(contains("_51")) %>%
  rowwise %>% 
  summarise(idx_d51 = any(c_across(everything()) > 2), .groups = "rowwise") %>%
  pull(idx_d51)

## select and order derived allele frequency vectors in one tibble
v20.dafs.d61 <- v20.dafs %>%
  select(1:3, contains(str_c(levels(repl.name), rep(c("_58", "_64"), each = 3)))) %>%
  filter(idx.d61)

v20.annots.d61 <- v20.annots %>% filter(idx.d61)



## Figure out for every mutation, how many mutations are in its vicinity. 
## I need an object that stores the number of nearby mutations
v20.nmutsnearby.d61 <- cbind(
  v20.dafs.d61[, 1:3], 
  matrix(NA, nrow=nrow(v20.dafs.d61), ncol = ncol(v20.dafs.d61) - 3)
)
## count the number of polymorphic sites within 
dist.nb <- 500 ## maximum distance (in bases) to be counted as a neighboring
  ## site
for(ic in 4:ncol(v20.nmutsnearby.d61)){
  for(ir in 1:nrow(v20.nmutsnearby.d61)){
    v20.nmutsnearby.d61[ir, ic] <- (v20.dafs.d61[ir, ic] > 0) *
      sum(v20.dafs.d61[abs(v20.dafs.d61$POS - v20.dafs.d61$POS[ir]) < dist.nb, ic] > 0, 
          na.rm=T)
  }
}


## write to output file (file included in publication)
# write_tsv(
#   v20.dafs %>% filter(idx.d61),
#   file = paste0(dir_v20, "fbs_dir/virus2020-1_fbs_filtered_split.daf")
# )
# writeLines(
#   v20.header,
#   con = paste0(dir_v20, "fbs_dir/virus2020-1_fbs_filtered_split.vcf")
# )
# write_tsv(
#   names(v20.vcf) %>% str_replace("CHROM", "#CHROM") %>% t %>% as_tibble(), 
#   append = TRUE,
#   file = paste0(dir_v20, "fbs_dir/virus2020-1_fbs_filtered_split.vcf")
# )
# write_tsv(
#   v20.vcf %>% filter(idx.d61), append = TRUE,
#   file = paste0(dir_v20, "fbs_dir/virus2020-1_fbs_filtered_split.vcf")
# )

## find the loci that are polymorphic at days 51 and 61 (:= 58 or 64)
###################


###################
## Create circus plots (visualise spatial distribution of genomic variation 
## along the genome, per replicate)
# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_circplot_snps_d61.pdf"), width = 8, height = 4)

# pdf(paste0(dir_db, "manuscript/figures/fig_allsnps_circplot_d61.pdf"), width = 8, height = 4)
par(mfrow=c(1, 2), mai = c(.1, .1, .1, .1))
## initialise
circos.clear()
circos.par('track.height' = .06, start.degree = 88, gap.degree = 4)
## plot no. 1; size shows frequency
circos.initialize(factors = v20.dafs.d61$CHROM, xlim=c(0, 330611))
for(ipop in 0:(ncol(v20.dafs.d61) - 3)){
  if(ipop == 0) { 
    ## create an empty plot
    with(v20.dafs.d61, circos.trackPlotRegion(
      factors = CHROM, x = POS, ylim=c(0, 1), panel.fun = NULL, bg.border = 0
    ))
    ## 'horizontal' axis on the outside
    circos.axis(
      h = "bottom", major.at = (0:7)*.5e5, minor.ticks = 4, 
      labels = paste0(as.character((0:6)*50), c("", rep("k", 6)))
    )
  }else{
    ## create a track, no borders
    with(v20.dafs.d61, circos.trackPlotRegion(
      factors = CHROM, x = POS, ylim=c(0, 1), panel.fun = NULL, bg.border = 0
    ))
    ## add line with a colour for treatment
    circos.trackLines(
      factors = v20.dafs.d61$CHROM[1:2], x = c(0, 330611), y = rep(.5, 2), 
      col = col_trts[ipop], lwd=2
    )
    ## add points for mutations, size to allele frequency
    
    with(filter(v20.dafs.d61, v20.dafs.d61[, ipop + 3] > 0), circos.trackPoints(
      factors = CHROM, x = POS, y = rep(.5, length(POS)), pch = 16, 
      col = 1, 
      cex = 1.5 * pull(filter(v20.dafs.d61, v20.dafs.d61[, ipop + 3] > 0), ipop + 3)^.25
    ))
  }
}

## plot no. 2; size corresponds to the number of nearby mutations
circos.initialize(factors = v20.dafs.d61$CHROM, xlim=c(0, 330611))
for(ipop in 0:(ncol(v20.dafs.d61) - 3)){
  if(ipop == 0) { 
    ## create an empty plot
    with(v20.dafs.d61, circos.trackPlotRegion(
      factors = CHROM, x = POS, ylim=c(0, 1), panel.fun = NULL, bg.border = 0
    ))
    ## 'horizontal' axis on the outside
    circos.axis(
      h = "bottom", major.at = (0:7)*.5e5, minor.ticks = 4, 
      labels = paste0(as.character((0:6)*50), c("", rep("k", 6)))
    )
  }else{
    ## create a track, no borders
    with(v20.dafs.d61, circos.trackPlotRegion(
      factors = CHROM, x = POS, ylim=c(0, 1), panel.fun = NULL, bg.border = 0
    ))
    ## add line with a colour for treatment
    circos.trackLines(
      factors = v20.dafs.d61$CHROM[1:2], x = c(0, 330611), y = rep(.5, 2), 
      col = col_trts[ipop], lwd=2
    )
    ## add points for mutations, size to nmutsnearby
    with(v20.dafs.d61[which(v20.dafs.d61[, ipop + 3] > 0), ], circos.trackPoints(
      factors = CHROM, x = POS, y = rep(.5, length(POS)), pch = 16, 
      col = 1, 
      cex = .8 * (v20.nmutsnearby.d61[which(v20.dafs.d61[, ipop + 3] > 0), ipop + 3])^(1/3)
    ))
  }
}
# dev.off(); system(paste0("open ", dir_db, "manuscript/figures/fig_allsnps_circplot_d61.pdf"))
# dev.off(); system(paste0("open ", dir_v20, "figures/variantscalled_freebayes/fig_circplot_snps_d61.pdf"))

# pdf(paste0(dir_db, "manuscript/figures/fig_allsnps_legend.pdf"), width=2, height=4)
plot(0, 0, type='n', xaxt='n', yaxt='n',  ylab="", xlab="", bty='n')
legend("topleft", legend = c(.01, .1, 1), bty = 'n', pch=16, 
       pt.cex = 1.5 * (c(.01, .1, 1))^.25)
legend("bottomleft", legend = c(1, 4, 14), bty = 'n', pch=16, 
       pt.cex = .8 * (c(1, 4, 14))^(1/3))
# dev.off(); system(paste0("open ", dir_db, "manuscript/figures/fig_allsnps_legend.pdf"))

## Create circus plots
###################


###################
## visualise annotations 
head(v20.annots)
with(v20.annots, table(MUT.CLASS))

fig.annots.d61.uniq <- ggplot() + geom_bar(
  data=v20.annots.d61, aes(x=CHROM, fill=MUT.SEVERITY), colour=I('black')
) + scale_fill_manual(
  name = "Predicted phenotypic effect", drop = F, values=col_ann
) + labs(x="", y="") + ylim(0, 70) + theme_bw() + theme(
  plot.margin = unit(c(.2, .1, .2, 0), "cm"), 
  panel.grid=element_blank(), 
  plot.title=element_text(colour="#E69F00"), 
  axis.text.y=element_text(size=24), axis.ticks.y=element_line(size=1), 
  axis.text.x=element_blank(), axis.ticks.x=element_blank()
)


## Then, plot these again, but multiplied by how often they're observed

## create a vector of the predicted effects, multiplied by the number of
## replicates in which they are found
v20.annots.d61.perrep <- data.frame(
  CHROM = rep(
    v20.annots$CHROM, times = apply(
      AO[, grep("_[56][84]$", colnames(AO))], 1, 
      (function(x1) sum(x1 > 2))
    )
  ), 
  POS = rep(
    v20.annots$POS, times = apply(
      AO[, grep("_[56][84]$", colnames(AO))], 1, 
      (function(x1) sum(x1 > 2))
    )
  ), 
  MUT.SEVERITY = rep(
    v20.annots$MUT.SEVERITY, times = apply(
      AO[, grep("_[56][84]$", colnames(AO))], 1, 
      (function(x1) sum(x1 > 2))
    )
  ), 
  GENE.NAME = rep(
    v20.annots$GENE.NAME, times = apply(
      AO[, grep("_[56][84]$", colnames(AO))], 1, 
      (function(x1) sum(x1 > 2))
    )
  )
)


fig.annots.d61.perrep <- ggplot() + geom_bar(
  data=v20.annots.d61.perrep, aes(x=CHROM, fill=MUT.SEVERITY), colour=I('black')
) + scale_fill_manual(
  name = "Predicted phenotypic effect", drop=FALSE, values=col_ann
) + labs(x="", y="") + ylim(0, 190) + theme_bw() + theme(
  plot.margin = unit(c(.2, .1, .2, 0), "cm"), 
  panel.grid=element_blank(), 
  plot.title=element_text(colour="#E69F00"), 
  axis.text.y=element_text(size=24), axis.ticks.y=element_line(size=1), 
  axis.text.x=element_blank(), axis.ticks.x=element_blank()
)


# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_virus2020_annots_uniqsnps.pdf"), height=3, width=4.42)
# pdf(paste0(dir_db, "manuscript/figures/fig_annots_uniqsnps.pdf"), height = 3, width = 4.42)
fig.annots.d61.uniq
# dev.off(); 
# system(paste0("open ", dir_v20, "figures/variantscalled_freebayes/fig_virus2020_annots_uniqsnps.pdf"))

# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_virus2020_annots_allsnps.pdf"), height=3, width=4)
fig.annots.d61.perrep
# dev.off(); system(paste0("open ", dir_v20, "figures/variantscalled_freebayes/fig_virus2020_annots_allsnps.pdf"))
## the fact that these bars look exactly the same suggests that 
## annotation does not predict repeatability.

## visualise annotations 
###################


###################
## compare genome-wide distribution of PPF vs the empirical one. 

## ! code below requires reading in an object of 9.9e5 rows and takes 
## ~10 minutes to run, so I commented it out
# allmuts.annots <- readVCFAnn(paste0(dir_v18, "sef_dir/pbcv1_allpossiblemuts_sef.vcf"))

# fig_annots_allmuts <- ggplot() + geom_bar(
#   data=allmuts.annots, aes(x=CHROM, fill=MUT.SEVERITY), colour=I('black')
# ) + scale_fill_manual(
#   name = "Predicted phenotypic effect",
#   values=col_ann[levels(allmuts.annots$MUT.SEVERITY) %in% as.character(allmuts.annots$MUT.SEVERITY)]
# ) + labs(x="", y="") + theme_bw() +
#   scale_y_continuous(limits=c(0, 1100000), breaks = c(0, 5, 10)*1e5) +  theme(
#     plot.margin = unit(c(.2, .1, .2, 0), "cm"),
#     panel.grid=element_blank(),
#     plot.title=element_text(colour="#E69F00"),
#     axis.text.y=element_text(size=24), axis.ticks.y=element_line(size=1),
#     axis.text.x=element_blank(), axis.ticks.x=element_blank()
#   )
# 
# pdf(paste0(dir_v20, "figures/fig_virus_annots_allmuts.pdf"), height=3, width=5)
# fig_annots_allmuts
# dev.off(); system(paste0("open ", dir_v20, "figures/fig_virus_annots_allmuts.pdf"))

# allmuts.annotcounts <- with(allmuts.annots, table(MUT.SEVERITY))
v20.annotcounts.d61.uniq <- with(v20.annots.d61, table(MUT.SEVERITY))
v20.annotcounts.d61.perrep <- with(v20.annots.d61.perrep, table(MUT.SEVERITY))
# 
# chisq.test(
#   c(allmuts.annotcounts, v20.annotcounts.perrep),
#   y=rep(with(allmuts.annots, levels(MUT.SEVERITY)), 2)
# )
# ## don't think this is the right test with more than two levels...

# write.table(
#   rbind(
#     allmuts.annotcounts,
#     v20.annotcounts.d61.uniq,
#     v20.annotcounts.d61.perrep
#   ), file=paste0(dir_db, "genome_data/virus2020-1_annotationcounts.d61.txt"),
#   quote=FALSE, row.names=TRUE, col.names=TRUE
# )

# v20.annotcounts.table <- read.table(
#   paste0(dir_db, "genome_data/virus2020-1_annotationcounts.d61.txt"), 
#   stringsAsFactors = F, header = TRUE,
# )


## careful that this read.table object is correctly interpreted as a 
## 2 by 4 numeric matrix by R. 
## Do a chi-squared test of independence between genome-wide and empirical
## distribution of fitness effects. 
# chisq.test(
#   as.matrix(v20.annotcounts.table[1:2, 1:4]), simulate.p.value = T, B = 1e6
# )
## The difference is not significant. The empirical distribution of PPE's could
## be a random draw from the set of all possible SNPs at every possible
## position of the reference genome. 
## also multipled by the number of replicates in which they're found
# chisq.test(
#   as.matrix(v20.annotcounts.table[c(1, 3), 1:4]), simulate.p.value = T, B = 1e6
# )
## This difference is also not significant. 

## compare genome-wide distribution of PPF vs the empirical one. 
###################


###################
## differences between experimental treatments
## population genetic diversity; really not informative with a handful of 
## polymorphic loci. 
head(v20.dafs.d61)
source('~/Documents/HVInt/scripts/calcGenDiv.R')
n.sites = 3e5
v20.pigw.d61 <- apply(
  v20.dafs.d61[, -(1:3)],
  2, (function(x1) weirsD(afs = c(x1, rep(0, n.sites - nrow(v20.dafs.d61)))))
)
names(v20.pigw.d61) <-
  colnames(v20.dafs.d61)[-(1:3)]
plot(1:6, v20.pigw.d61, pch = 16, cex = 2, col = col_trts)
t.test(
  v20.pigw.d61 ~ factor(substr(colnames(v20.dafs.d61)[-(1:3)], 1, 3), levels = levels(treatment))
)
wilcox.test(
  v20.pigw.d61 ~ factor(substr(colnames(v20.dafs.d61)[-(1:3)], 1, 3), levels = levels(treatment))
)
## Nope, no difference despite the fact that diversity is consistenly lower
## in the DMS populations :-O. The same for the non-parametric equivalent 
## (wilcox.test() )

## number of variable sites
v20.nvarsites.d61 <- apply(
  AO[, grep("_[56][84]$", colnames(AO))], 2,
  (function(x1) sum(x1 > 2))
)[c(3, 1, 5, 2, 6, 4)]
v20.nvarsites.d61 <- apply(
  v20.dafs.d61[, -(1:3)], 2, (function(x1) sum(x1 > 0, na.rm = T))
)
 
plot(1:6, v20.nvarsites.d61, pch = 16, cex = 2, col = col_trts)
t.test(
  v20.nvarsites.d61 ~ factor(substr(colnames(v20.dafs.d61)[-(1:3)], 1, 3), levels = levels(treatment))
)
m1 <- glm(
  v20.nvarsites.d61 ~ factor(substr(colnames(v20.dafs.d61)[-(1:3)], 1, 3), levels = levels(treatment)), family = poisson(link = "log")
)
summary(m1)
## differences between experimental treatments
###################

###################
## Correlate predicted effect with repeatability and frequency change
## add a column to annots with the number of replicates in which a 
## SNP is found
v20.annots.d61 <- v20.annots.d61 %>%
  mutate(
    N.REPS.FOUND = v20.dafs.d61 %>%
      select(-(1:3)) %>% rowwise %>% 
      summarise(nrf = sum(across(everything()) > 0, na.rm = T)) %>% 
      pull(nrf)
  )

## also calculate selection coefficient. 
## selection coefficient would be a normalised transformation of 
## allele frequency change (assuming constant selection)
source('~/Documents/HVInt/scripts/calcSEffective.R')
v20.seff.d61 <- matrix(
  NA, nrow=nrow(v20.dafs.d61), ncol=ncol(v20.dafs.d61) - 3
)
for(ir in 1:nrow(v20.seff.d61)){
  for(ic in 1:ncol(v20.seff.d61)){
    v20.seff.d61[ir, ic] <- calcSEffective(
      afs=c(1e-4, v20.dafs.d61[ir, ic + 3] %>% pull), 
      tps=c(1, rep(c(58, 64), each = 3)[ic]), includezeroes=F
    )
  }
}
v20.annots.d61 <- v20.annots.d61 %>% 
  mutate(MEAN.SEFF = rowMeans(v20.seff.d61, na.rm = T))

## number of replicates in which a mutation reaches 1% frequency: 
v20.nrepsfound.1p <- apply(
  v20.dafs.d61[, -(1:3)], 1, (function(x1) sum(x1 > 0.01, na.rm=T))
)
# system(paste0("open ", dir_v20, "figures/variantscalled_freebayes/fig_*_muts.pdf"))

# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_nrepsvsannot.pdf"), width=4, height=3)
# pdf(paste0(dir_db, "manuscript/figures/fig_nrepsvsannots.pdf"), width = 5, height = 3)
par(mai=c(.22, .22, .22, .22))
with(v20.annots.d61, plot(
    N.REPS.FOUND ~ MUT.SEVERITY, xlim = c(.5, 4.5), bty = 'n'
))

set.seed(2307); with(v20.annots.d61, points(
  x = jitter(as.numeric(MUT.SEVERITY), factor=1), 
  y = N.REPS.FOUND, pch=16, col=1
))
# dev.off(); system(paste0("open ", dir_v20, "figures/variantscalled_freebayes/fig_nrepsvsannot.pdf"))

# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_nrepsvsannots_1perc.pdf"), width = 5, height = 3)
par(mai=c(.82, .82, .22, .22))
plot(v20.nrepsfound.1p ~ v20.annots.d61$MUT.SEVERITY, xlim=c(.5, 4.5), 
     xlab = "Predicted phenotypic effect", ylab = "Number of replicates", bty = 'n')
set.seed(2307); points(
  x = jitter(as.numeric(v20.annots.d61$MUT.SEVERITY), factor=1), 
  y = v20.nrepsfound.1p, pch=16, col=1
)
# dev.off(); system(paste0("open ", dir_v20, "figures/variantscalled_freebayes/fig_nrepsvsannots_1perc.pdf"))

## plot effective selection coefficient against predicted phenotypic effect
# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_seffvsannot.pdf"), width=4, height=3)
par(mai=c(.22, .22, .22, .22))
plot(
  as.numeric(v20.seff.d61) ~ 
    rep(v20.annots.d61$MUT.SEVERITY, ncol(v20.seff.d61)), 
  xlim=c(.5, 4.5), ylim=c(0, .28), bty = 'n'
)
set.seed(2307); points(
  x = jitter(as.numeric(
    rep(v20.annots.d61$MUT.SEVERITY, ncol(v20.seff.d61))
  ), factor=1), 
  y = as.numeric(v20.seff.d61), pch=16, col=alpha(1, .667), cex = .6
)
# dev.off(); system(paste0("open ", dir_v20, "figures/variantscalled_freebayes/fig_seffvsannot.pdf"))



## I want to see average selection coefficient per mutation, across
## six replicates. 
# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_meanseffvsannot.pdf"), width=4, height=3)
# pdf(paste0(dir_db, "manuscript/figures/fig_meanseffvsannot.pdf"), width = 5, height = 3)
par(mai=c(.22, .22, .22, .22))
with(v20.annots.d61, plot(
  MEAN.SEFF ~ MUT.SEVERITY, xlim = c(.5, 4.5), ylim = c(0, .19), bty = 'n'
))

set.seed(2307); with(v20.annots.d61, points(
  y = MEAN.SEFF, x = jitter(as.numeric(MUT.SEVERITY, factor = 1)), 
  pch=16, col=alpha(1, .667), cex = .8, 
))

# dev.off()
# dev.off(); system(paste0("open ", dir_v20, "figures/variantscalled_freebayes/fig_meanseffvsannot.pdf"))



## Stats: linear regression of no. of repls.found over mut.severity
## on the full dataset. 
lm(N.REPS.FOUND ~ as.numeric(MUT.SEVERITY), data = v20.annots.d61) %>% 
  summary
## this isn't right, because we expect mutations to have a phenotypic effect
## or not; do a Mann-Whitney U test with non-synonymous and synonymous SNPs only 
filter(v20.annots.d61, MUT.CLASS %in% c("synonymous_variant", "missense_variant")) %>%
  ggplot() + geom_jitter(aes(x = MUT.CLASS, y = N.REPS.FOUND), height = 0, width = .2)

wilcox.test(
  N.REPS.FOUND ~ MUT.SEVERITY, 
  data = filter(v20.annots.d61, MUT.CLASS %in% c("synonymous_variant", "missense_variant"))
)

## t-test to check for differences in average selection coefficient 
## (calculated across six replicates): 
filter(v20.annots.d61, MUT.CLASS %in% c("synonymous_variant", "missense_variant")) %>%
  ggplot() + geom_jitter(aes(x = MUT.CLASS, y = MEAN.SEFF), height = 0, width = .2)
t.test(
  MEAN.SEFF ~ MUT.SEVERITY, 
  data = filter(v20.annots.d61, MUT.CLASS %in% c("synonymous_variant", "missense_variant"))
)



## Maybe easier to comprehend than s_eff vs ppe: 
## plot allele frequency versus predicted phenotypic effect
## only include values larger than 0
idx_notzero_nva <- as.numeric(as.matrix(v20.dafs.d61[, -(1:3)])) > 0
# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_afvsannot.pdf"), width=4, height=3)
dev.off(); par(mai=c(.22, .22, .22, .22))
plot(
  as.numeric(as.matrix(v20.dafs.d61[, -(1:3)]))[idx_notzero_nva] ~ 
    rep(v20.annots.d61$MUT.SEVERITY, ncol(v20.dafs.d61) - 3)[idx_notzero_nva], 
    xlim=c(.5, 4.5), ylim=c(0, 1.02), bty = 'n'
)
set.seed(2307); points(
  x = jitter(as.numeric(
    rep(v20.annots.d61$MUT.SEVERITY, ncol(v20.dafs.d61) - 3)
  ), factor=1)[idx_notzero_nva], 
  y = as.numeric(as.matrix(v20.dafs.d61[, -(1:3)]))[idx_notzero_nva], 
  pch=16, col=alpha(1, .667)
)
# dev.off(); system(paste0("open ", dir_v20, "figures/variantscalled_freebayes/fig_afvsannot.pdf"))

## I want to see allele frequency versus the number of replicates
# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_nrepsvsseff_1perc.pdf"), width=4, height=3)
par(mai = c(.82, .82, .22, .22))
plot(
  v20.nrepsfound.1p ~ rowMeans(v20.seff.d61, na.rm = T), 
  ylab = "Number of replicates", xlab = "mean s_eff"
)
# dev.off(); system(paste0("open ", dir_v20, "figures/variantscalled_freebayes/fig_nrepsvsseff_1perc.pdf"))
## Correlate predicted effect to repeatability
###################


###################
## Assess if certain genes are hit more often than expected. 
with(v20.annots.d61, table(GENE.NAME))
with(v20.annots.d61.perrep, table(GENE.NAME))

## Those are the raw numbers. I plan to do Fisher's exact test, 
## which requires: 
## - genome length       - number of mutations observed
## - gene length         - number of mutations in gene
## Correct for multiple-testing by dividing the p-value by the number
## of genes in the viral genome. 

## generate a data frame to store this information per gene
v20.gff <- read.table(
  "~/Documents/HVInt/Chlorella/ref/snpeff/data/PBCV1/genes.gff", 
  sep="\t", stringsAsFactors = FALSE, comment.char="#"
)
## copied names from https://www.ensembl.org/info/website/upload/gff.html?redirect=no: 
names(v20.gff) <- c("CHROM", "SOURCE", "FEATURE", "START", "END", "SCORE", "STRAND", "FRAME", "ATTR")
v20.genes <- v20.gff[v20.gff$FEATURE == "gene", ]
v20.genes$gene.name <- v20.genes$ATTR %>%
  (function(x1) strsplit(x1, split=";")) %>%
  (function(x2) sapply(x2, "[[", 3)) %>%
  (function(x3) strsplit(x3, split="=")) %>%
  (function(x2) sapply(x2, "[[", 2))

## extract gene length 
v20.genes$gene.length <- with(v20.genes, 1 + (END - START))

## for every gene, count how often it is hit
timesInAnnots <- function(genename, annot.obj){
  ## return the number of times {genename} is present in 
  ## annot.obj$GENE.NAME. Function exists for vectorization purposes. 
  return(sum(genename == annot.obj$GENE.NAME))
}; timesInAnnots <- Vectorize(FUN = timesInAnnots, vectorize.args="genename")

v20.genes$hit.d61.uniq <- timesInAnnots(
  v20.genes$gene.name, annot.obj = v20.annots.d61
)
v20.genes$hit.d61.perrep <- timesInAnnots(
  v20.genes$gene.name, annot.obj = v20.annots.d61.perrep
)

## conservative (because not the full genome is evaluated) 
## estimate of genome length. 
v20.genomelength <- round(v20.gff[v20.gff$FEATURE == "region", "END"] * 0.99)

FET.pbcv1 <- function(genome.length, gene.length, nmuts.obs, nmuts.ingene){
  ## given input genome length, gene lengths, number of observed mutations
  ## and number of observed mutations per gene, 
  ## returns a (non-corrected) p-value of Fisher's exact test
  ## [...]
  fisher.test(
    x=matrix(c(
      nmuts.obs - nmuts.ingene, nmuts.ingene, 
      genome.length - gene.length, gene.length
    ), nrow=2)
  )$p.value
}
v20.genes$fet.pval.uniq <- numeric(nrow(v20.genes))
v20.genes$fet.pval.perrep <- numeric(nrow(v20.genes))

for(i in 1:nrow(v20.genes)){
  v20.genes$fet.pval.uniq[i] <- FET.pbcv1(
    genome.length = v20.genomelength, 
    gene.length = v20.genes[i, "gene.length"], 
    nmuts.obs = sum(v20.genes[, "hit.d61.uniq"]), 
    nmuts.ingene = v20.genes[i, "hit.d61.uniq"]
  )
  v20.genes$fet.pval.perrep[i] <- FET.pbcv1(
    genome.length = v20.genomelength * 6, 
    gene.length = v20.genes[i, "gene.length"] * 6, 
    nmuts.obs = sum(v20.genes[, "hit.d61.perrep"]), 
    nmuts.ingene = v20.genes[i, "hit.d61.perrep"]
  )
}
## done, I believe. 


with(v20.genes, (sort(fet.pval.uniq) * nrow(v20.genes)) %>% head)
with(v20.genes, (sort(fet.pval.perrep) * nrow(v20.genes)) %>% head)
## For the record, I don't think we should use pval.perrep in this case. 
## Evaluating a stretch of 10 positions in 6 evolutionary replicates 
## can hardly be considered 60 independent observations of 
## mutation occurrence. 

## control for FWER for this set of tests: 
v20.genes$fet.correctedpval.uniq <- p.adjust(
  v20.genes$fet.pval.uniq, method="holm"
)
sum(v20.genes$fet.correctedpval.uniq < .05)
mean(v20.genes$fet.correctedpval.uniq < .05)

v20.genes$fet.correctedpval.perrep <- p.adjust(
  v20.genes$fet.pval.perrep, method="holm"
)
sum(v20.genes$fet.correctedpval.perrep < .05)
mean(v20.genes$fet.correctedpval.perrep < .05)
## so as more or less expected, this confirms that observing the same 
## mutation in multiple replicates is highly unlikely. 

v20.genes[v20.genes$fet.correctedpval.uniq < .05, ]
# v20.genes[v20.genes$fet.correctedpval.perrep < .05, ]

## need to know how many PBCV-1 reference genome positions are 
## not in any protein-coding region. 
v.pos <- 1:v20.gff[v20.gff$FEATURE == "region", "END"]
for(i in 1:nrow(v20.genes)){
  v.pos <- v.pos[!between(v.pos, v20.genes[i, "START"], v20.genes[i, "END"])]
}
length(v.pos); length(v.pos) / v20.genomelength
## This matches what genome-wide snpeff investigation tells me
## I reported 20% in the Science Advances paper :$. 

# write.table(
#   v20.genes[(v20.genes$fet.correctedpval.uniq < .05), c(1, 4, 5, 9:17)],
#   file = paste0(dir_db, "genome_data/virus2020-1_geneshitoften.txt"),
#   quote = F, row.names=F, col.names = T
# )
# write.table(
#   v20.genes[(v20.genes$hit.d61.uniq > 0), c(1, 4, 5, 9:17)],
#   file = paste0(dir_db, "genome_data/virus2020-1_geneshit.txt"),
#   quote = F, row.names=F, col.names = T
# )

## Assess if certain genes are hit more often than expected. 
###################

###################
## zoom in on the mutations found in A122/123R, A140/145R and A540L
## first plot: gene A122R
A122R.coords <- c(62145, 66176)
A122R.xlim <- c(65080, 65420)
A140R.coords <- c(73107,	76490)
A140R.xlim <- c(75300, 76300)
A540L.coords <- c(257089, 260859)
A540L.xlim <- c(257200, 258800)

# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_A122R_muts.pdf"), height = 4, width = 4)
set.seed(2303)
par(mai = c(.22, .22, .22, .22))
plot(A122R.xlim, c(0, 7), type='n', xaxt='n', yaxt='n', 
     ylab="", xlab="", bty='n')
axis(side = 1, at = 100 * (600 : 700), las = 1, labels = F)
for(i in 4:ncol(v20.dafs.d61)){
  ## horizontal line through the plot, coloured by treatment
  lines(A122R.xlim, rep((ncol(v20.dafs.d61):1)[i], 2), col=col_trts[i - 3], lwd=3)
  ## beginning and ending of gene
  points(A122R.coords, rep((ncol(v20.dafs.d61):1)[i], 2), pch='|', col=1)
  ## dots where the mutations are; size reflecting frequency, colour
  ## giving predicted phenotypic effect
  posafs.ingene <- v20.dafs.d61[
    between(v20.dafs.d61$POS, A122R.coords[1], A122R.coords[2]) & (v20.dafs.d61[, i] > 0.01), 
    c(3, 2, i)
  ]
  points(
    x = pull(posafs.ingene, 2), 
    y = rep((ncol(v20.dafs.d61):1)[i], nrow(posafs.ingene)), 
    col = col_ann[pull(posafs.ingene, 1)], pch = 16, cex = 3 * (pull(posafs.ingene, 3)^.25)
  )
  rm(posafs.ingene)
}
# dev.off()

# pdf(paste0(dir_v20, "figures/fig_A122R_zoom.pdf"), height = 1, width = 4)
par(mai = c(.22, .22, .22, .22))
plot(
  A122R.coords, rep(.62, 2), type='p', xaxt='n', yaxt='n',ylab="", 
  xlab="", bty='n', pch=16, cex=.2, ylim=c(0, 1), 
  xlim = 1000 * c(61.8, 66.4)
)
axis(side = 1, at = 1000 * (61:67), las = 1, labels = F)
points(A122R.xlim, rep(.6, 2), col = col_brred, pch = '|', cex = 2)
# dev.off(); system(paste0("open ", dir_v20, "figures/fig_A122R_zoom.pdf"))


## second plot: A140/145R

# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_A140R_muts.pdf"), height = 4, width = 4)
set.seed(2306)
par(mai = c(.22, .22, .22, .22))
plot(A140R.xlim, c(0, 7), type='n', xaxt='n', yaxt='n', 
     ylab="", xlab="", bty='n')
axis(side = 1, at = 500 * (148 : 156), las = 1, labels = F)
for(i in 4:ncol(v20.dafs.d61)){
  ## horizontal line through the plot, coloured by treatment
  lines(A140R.xlim, rep((ncol(v20.dafs.d61):1)[i], 2), col=col_trts[i - 3], lwd=3)
  ## beginning and ending of gene
  points(A140R.coords, rep((ncol(v20.dafs.d61):1)[i], 2), pch='|', col=1)
  ## dots where the mutations are; size reflecting frequency, colour
  ## giving predicted phenotypic effect
  posafs.ingene <- v20.dafs.d61[
    between(v20.dafs.d61$POS, A140R.coords[1], A140R.coords[2]) & (v20.dafs.d61[, i] > 0.01), 
    c(3, 2, i)
  ]
  points(
    x = pull(posafs.ingene, 2), 
    y = jitter(rep((ncol(v20.dafs.d61):1)[i], nrow(posafs.ingene)), factor = 3), 
    col = col_ann[pull(posafs.ingene, 1)], pch = 16, cex = 3 * (pull(posafs.ingene, 3)^.25)
  )
  rm(posafs.ingene)
}
# dev.off()

# pdf(paste0(dir_v20, "figures/fig_A140R_zoom.pdf"), height = 1, width = 4)
par(mai = c(.22, .22, .22, .22))
plot(
  A140R.coords, rep(.62, 2), type='p', xaxt='n', yaxt='n',ylab="", 
  xlab="", bty='n', pch=16, cex=.2, ylim=c(0, 1), 
  xlim = 1000 * c(72.8, 76.6)
)
axis(side = 1, at = 1000 * (71:77), las = 1, labels = F)
points(A140R.xlim, rep(.6, 2), col = col_brred, pch = '|', cex = 2)
# dev.off(); system(paste0("open ", dir_v20, "figures/fig_A140R_zoom.pdf"))



## third plot: gene A540L
# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_A540R_muts.pdf"), height = 4, width = 4)
set.seed(2403)
par(mai = c(.22, .22, .22, .22))
plot(A540L.xlim, c(0, 7), type='n', xaxt='n', yaxt='n', 
     ylab="", xlab="", bty='n')
axis(side = 1, at = 500 * (510 : 520), las = 1, labels = F)
for(i in 4:ncol(v20.dafs.d61)){
  ## horizontal line through the plot, coloured by treatment
  lines(A540L.xlim, rep((ncol(v20.dafs.d61):1)[i], 2), col=col_trts[i - 3], lwd=3)
  ## beginning and ending of gene
  points(A540L.coords, rep((ncol(v20.dafs.d61):1)[i], 2), pch='|', col=1)
  ## dots where the mutations are; size reflecting frequency, colour
  ## giving predicted phenotypic effect
  posafs.ingene <- v20.dafs.d61[
    between(v20.dafs.d61$POS, A540L.coords[1], A540L.coords[2]) & (v20.dafs.d61[, i] > 0.01), 
    c(3, 2, i)
  ]
  points(
    x = pull(posafs.ingene, 2), 
    y = jitter(rep((ncol(v20.dafs.d61):1)[i], nrow(posafs.ingene)), factor = 2), 
    col = col_ann[pull(posafs.ingene, 1)], pch = 16, cex = 3 * (pull(posafs.ingene, 3)^.25)
  )
  rm(posafs.ingene)
}
# dev.off()
# system(paste0("open ", dir_v20, "figures/variantscalled_freebayes/fig_*_muts.pdf"))


# pdf(paste0(dir_v20, "figures/fig_A540L_zoom.pdf"), height = 1, width = 4)
par(mai = c(.22, .22, .22, .22))
plot(
  A540L.coords, rep(.62, 2), type='p', xaxt='n', yaxt='n',ylab="", 
  xlab="", bty='n', pch=16, cex=.2, ylim=c(0, 1), 
  xlim = 1000 * c(256.7, 261.3)
)
axis(side = 1, at = 1000 * (256:262), las = 1, labels = F)
points(A540L.xlim, rep(.6, 2), col = 2, pch = '|', cex = 2)
# dev.off(); system(paste0("open ", dir_v20, "figures/fig_A540L_zoom.pdf"))

# pdf(paste0(dir_v20, "figures/fig_genes_sizelegend.pdf"), height = 2.4, width = 2)
par(mai = c(.22, .22, .22, .22))
plot(
  rep(1, 3), 1:3, pch=16, cex = 3 * (10^(0:-2))^(.25), 
  bty = 'n', xaxt = 'n', yaxt = 'n', xlim=c(.5, 3), ylim=c(-1, 7)
)
text(rep(2, 3), 1:3, labels = as.character(10^-(0:2)))
# dev.off(); system(paste0("open ", dir_v20, "figures/fig_genes_sizelegend.pdf"))
## zoom in on the mutations found in A122/123R and A540L
###################



###################
## Plot mutations found in A122/123R, A140/145R and A540L over time
## Only include the ones that reach at least 1% frequency
head(v20.dafs)
idx.A122R <- between(v20.dafs$POS, A122R.coords[1], A122R.coords[2]) & 
  apply(v20.dafs[, -(1:3)], 1, (function(x1) any(x1 > 0.01, na.rm = TRUE)))
idx.A140R <- between(v20.dafs$POS, A140R.coords[1], A140R.coords[2]) & 
  apply(v20.dafs[, -(1:3)], 1, (function(x1) any(x1 > 0.01, na.rm = TRUE)))
idx.A540L <- between(v20.dafs$POS, A540L.coords[1], A540L.coords[2]) & 
  apply(v20.dafs[, -(1:3)], 1, (function(x1) any(x1 > 0.01, na.rm = TRUE)))
n.reps <- nlevels(repl.name)


reachesOnePerc <- function(afs){
  return(any(afs > 0.01, na.rm = TRUE))
}
reachesOnePercTail <- function(afs, remove_first = 1){
  if((!is.double(remove_first)) | (length(remove_first) > 1) | (remove_first < 1)){
    stop("argument remove_first is not a length-one positive integer")
  }
  return(any(afs[-(1:remove_first)] > 0.01, na.rm = TRUE))
}
firstLast <- function(x1){ return(c(x1[1], x1[length(x1)])) }


## 
v20.dafs_hotspots <- v20.dafs %>% 
  mutate(
    in_A122R = between(POS, A122R.coords[1], A122R.coords[2]), 
    in_A140R = between(POS, A140R.coords[1], A140R.coords[2]), 
    in_A540L = between(POS, A540L.coords[1], A540L.coords[2]), 
    POS_uniq = with(v20.dafs, str_c(POS, c("", "_2")[(1 + duplicated(POS))]))
  ) %>%
  rowwise() %>% 
  mutate(
    reaches_1p = reachesOnePerc(c_across(str_subset(names(v20.dafs), "(SEL)|(DMS)"))), 
  )

## generate random rainbow colour palettes
set.seed(2303)
col.A122R <- sample(rainbow(with(v20.dafs_hotspots, sum(reaches_1p & in_A122R))))
names(col.A122R) <- v20.dafs_hotspots %>% 
  filter(reaches_1p & in_A122R) %>% pull(POS_uniq)
col.A540L <- sample(rainbow(with(v20.dafs_hotspots, sum(reaches_1p & in_A540L))))
names(col.A540L) <- v20.dafs_hotspots %>% 
  filter(reaches_1p & in_A540L) %>% pull(POS_uniq)
col.A140R <- sample(rainbow(with(v20.dafs_hotspots, sum(reaches_1p & in_A140R))))
names(col.A140R) <- v20.dafs_hotspots %>% 
  filter(reaches_1p & in_A140R) %>% pull(POS_uniq)




# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_afsovertime_A122R.pdf"), width = 9, height = 6)
## generate plots for A122R: 
layout(t(matrix(c(1:6), nrow = 3)))
for(repli in levels(repl.name)){
  ## per replicate, filter the SNPs that reach 1% frequency at least once:  
  v20.dafs_repli_orf <- v20.dafs_hotspots %>% 
    filter(in_A122R) %>% 
    select(contains("POS"), contains(repli)) %>% 
    rowwise %>%
    mutate(
      reaches_1p_repli = 
        reachesOnePerc(
        c_across(str_subset(names(v20.dafs_hotspots), repli))
      )
    ) %>% 
    filter(reaches_1p_repli)
  ## extract time point info from variable names
  tps_repli <- str_subset(names(v20.dafs_repli_orf), repli) %>% 
    str_replace(pattern = "\\w+_\\d_(\\d+)", replacement = "\\1") %>%
    as.numeric()
  ## create plot, order by time point
  plotAfs(
    af_mat = select(v20.dafs_repli_orf, contains(repli))[, order(tps_repli)] %>% 
      as.matrix, 
    tps = tps_repli[order(tps_repli)], 
    col_v = col.A122R[pull(v20.dafs_repli_orf, POS_uniq)], 
    xlim_v = c(0, 67), vline_sampletimes = F, 
    title_v = str_c("SNPs inside A122R - ", repli)
  )
  
}
# dev.off()

# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_afsovertime_A140R.pdf"), width = 9, height = 6)
## generate plots for A140R: 
layout(t(matrix(c(1:6), nrow = 3)))
for(repli in levels(repl.name)){
  ## per replicate, filter the SNPs that reach 1% frequency at least once:  
  v20.dafs_repli_orf <- v20.dafs_hotspots %>% 
    filter(in_A140R) %>% 
    select(contains("POS"), contains(repli)) %>% 
    rowwise %>%
    mutate(
      reaches_1p_repli = reachesOnePerc(
        c_across(str_subset(names(v20.dafs_hotspots), repli))
      )
    ) %>% 
    filter(reaches_1p_repli)
  ## extract time point info from variable names
  tps_repli <- str_subset(names(v20.dafs_repli_orf), repli) %>% 
    str_replace(pattern = "\\w+_\\d_(\\d+)", replacement = "\\1") %>%
    as.numeric()
  ## create plot, order by time point
  plotAfs(
    af_mat = select(v20.dafs_repli_orf, contains(repli))[, order(tps_repli)] %>% 
      as.matrix, 
    tps = tps_repli[order(tps_repli)], 
    col_v = col.A140R[pull(v20.dafs_repli_orf, POS_uniq)], 
    xlim_v = c(0, 67), vline_sampletimes = F, 
    title_v = str_c("SNPs inside A140R - ", repli)
  )
  
}
# dev.off()

# pdf(paste0(dir_v20, "figures/variantscalled_freebayes/fig_afsovertime_A540L.pdf"), width = 9, height = 6)
## generate plots for A540L: 
layout(t(matrix(c(1:6), nrow = 3)))
for(repli in levels(repl.name)){
  ## per replicate, filter the SNPs that reach 1% frequency at least once:  
  v20.dafs_repli_orf <- v20.dafs_hotspots %>% 
    filter(in_A540L) %>% 
    select(contains("POS"), contains(repli)) %>% 
    rowwise %>%
    mutate(
      reaches_1p_repli = reachesOnePerc(
        c_across(str_subset(names(v20.dafs_hotspots), repli))
      )
    ) %>% 
    filter(reaches_1p_repli)
  ## extract time point info from variable names
  tps_repli <- str_subset(names(v20.dafs_repli_orf), repli) %>% 
    str_replace(pattern = "\\w+_\\d_(\\d+)", replacement = "\\1") %>%
    as.numeric()
  ## create plot, order by time point
  plotAfs(
    af_mat = select(v20.dafs_repli_orf, contains(repli))[, order(tps_repli)] %>% 
      as.matrix, 
    tps = tps_repli[order(tps_repli)], 
    col_v = col.A540L[pull(v20.dafs_repli_orf, POS_uniq)], 
    xlim_v = c(0, 67), vline_sampletimes = F, 
    title_v = str_c("SNPs inside A540L - ", repli)
  )
  
}
# dev.off()




## Create publication-quality figures: see 
## ~/Dropbox/Talks_Algae_Virus/Project3_chlorovirusmuts/manuscript/figures/plot_code.R 

## I'm curious to see where the second SNP at position 65130 has gone, 
## because I can't immediately locate it on the figures: 
v20.dafs[v20.dafs$POS == 65130, ]
## It's just very low-frequency. 

## for the publication, I'd like to know the frequency of the SNP that 
## reaches the highest overall frequency (it's in SEL_2): 
v20.dafs %>%
  filter(between(POS, A140R.coords[1], A140R.coords[2])) %>% 
  pivot_longer(
    cols = str_subset(names(v20.dafs), "(SEL)|(DMS)"), 
    names_to = "pop.name", values_to = "AF"
  ) %>% arrange(desc(AF)) %>% head
## I would also like to know the maximum frequency per replicate: 
v20.dafs %>% 
  filter(between(POS, A140R.coords[1], A140R.coords[2])) %>% 
  pivot_longer(
    cols = str_subset(names(v20.dafs), "(SEL)|(DMS)"), 
    names_to = "pop.name", values_to = "AF"
  ) %>% mutate(repl.name = str_sub(pop.name, 1, 5)) %>%
  group_by(repl.name) %>%
  summarise(
    max.AF = max(AF, na.rm = TRUE)
  )

## how much synonymous and non-synonymous SNPs are found in A140/145R? 
v20.annots.d61 %>%
  filter(between(POS, A140R.coords[1], A140R.coords[2])) %>%
  count(MUT.CLASS)
  

## Plot mutations found in A122/123R and A540L over time
###################