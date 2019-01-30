#### TCGA Pan Cancer Analysis
#### Immune
####
#### Mutant peptide and TCR co-occurrence
#### March 21, 2017
####
#### Scott Brown
#### sbrown@bcgsc.ca


## Read in data

setwd("/projects/sbrown_prj/TCGA/PanCanAtlas/analysis/figures/")

dat <- read.table("../mut_tcr_cooccurrence_170322.tsv", header=T, stringsAsFactors=F, sep="\t")

numSamp <- 7917
## from running and getting sample count:
#(python3) [sbrown@gphost03 /projects/sbrown_prj/TCGA/PanCanAtlas/analysis]
#λ python ../scripts/test_TCR_neo_cooccurrence.py TCGA_peptideMHC_binding_MC3_v0.2.8.CONTROLLED_predictedBindersOnly_withExpression_170317.tsv ../data/mitcr_complete/TCGA_mitcr_cdr3_result_161008.tsv ../data/mitcr_complete/mitcr_sampleStatistics_20160714.tsv mut_tcr_cooccurrence_170322.tsv

## get ratio of samples (how shared is it)
dat$mutRatio <- dat$found_together / dat$mutation_count
dat$tcrRatio <- dat$found_together / dat$TCR_count

## TCR-pMHC co-occurrence

pmhccdr3 <- subset(dat, type=="pmhc-cdr3")

pmhccdr3$fisher.p <- apply(pmhccdr3, 1, function(x){
  a <- as.numeric(x[["found_together"]])
  b <- as.numeric(x[["mutation_count"]]) - a
  c <- as.numeric(x[["TCR_count"]]) - a
  d <- numSamp - a - b - c
  return(fisher.test(matrix(c(a,b,c,d), nrow=2, byrow=T))$p.value)
})

pmhccdr3$chain <- apply(pmhccdr3, 1, function(x){
  if(grepl("alpha",x[["TCR"]])){
    return("alpha")
  }else{
    return("beta")
  }
})

#### Figure 5E ####
p <- ggplot(subset(pmhccdr3, fisher.p<0.05/nrow(pmhccdr3)), aes(x=mutRatio, y=tcrRatio, color=-log10(fisher.p*nrow(pmhccdr3)))) + geom_jitter(height=0.01, width=0.01, aes(size=found_together), shape=1, stroke=1.3, alpha=0.8) + scale_size(range=c(2,7), name="Number of Samples with pair") + scale_color_distiller(palette="Spectral", name="-log10(Adjusted Fisher's p)") + facet_wrap(~chain, strip.position = "left") + theme_bw() + xlab("(Samples with both) / (Samples with pMHC)") + ylab("(Samples with both) / (Samples with CDR3)") + ggtitle("Co-occurrence of CDR3 and pMHC")

ggsave("pmhc_cdr3_cooccur_withSize_byChain_sig_180117.pdf", p, width=10, height=4)


#### Supplemental Table - TCR-pMHC co-occurrence ####

write.table(pmhccdr3[order(pmhccdr3$found_together, decreasing = T),], "ST3_TCR-pMHC.tsv", quote = F, sep = "\t", row.names = F)


## TCR alpha-beta co-occurrence
alphabeta <- subset(dat, type=="cdr3a-cdr3b")
alphabeta$fisher.p <- apply(alphabeta, 1, function(x){
  a <- as.numeric(x[["found_together"]])
  b <- as.numeric(x[["mutation_count"]]) - a
  c <- as.numeric(x[["TCR_count"]]) - a
  d <- numSamp - a - b - c
  return(fisher.test(matrix(c(a,b,c,d), nrow=2, byrow=T))$p.value)
})

#### Figure 5D ####
p <- ggplot(subset(alphabeta, fisher.p<0.05/nrow(alphabeta)), aes(x=mutRatio, y=tcrRatio, color=-log10(fisher.p*nrow(alphabeta)))) + geom_jitter(height=0.01, width=0.01, aes(size=found_together), shape=1, stroke=1.3, alpha=0.8) + scale_size(range=c(2,12), name="Number of Samples with pair") + scale_color_distiller(palette="Spectral", name="-log10(Adjusted Fisher's p)") + theme_bw() + xlab("(Samples with both) / (Samples with CDR3a)") + ylab("(Samples with both) / (Samples with CDR3b)") + ggtitle("Co-occurrence of CDR3a and CDR3b")

ggsave("cdr3a_cdr3b_cooccur_withSize_byChain_sig_180117.pdf", p, width=7, height=5)


#### Supplemental Tabel - TCR alpha-beta co-occurrence ####

write.table(alphabeta[order(alphabeta$found_together, decreasing = T),], "ST3_TCRa-TCRb.tsv", quote = F, sep = "\t", row.names = F)

