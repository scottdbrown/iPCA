#### TCGA Pan Cancer Analysis
#### Immune
####
#### Recurrent neoantigens
#### October 13, 2016
####
#### Scott Brown
#### sbrown@bcgsc.ca


## Read in data

setwd("/projects/sbrown_prj/TCGA/PanCanAtlas/analysis/figures/")

pmhcs <- read.table("recurrentPMHC_170530.tsv", header=T, sep="\t", stringsAsFactors=F)

##pmhcs
pmhcs$x <- 1:nrow(pmhcs)

rpmhcs <- subset(pmhcs, numSub>1)

#### Figure 5A ####

pmhcs_nmut <- do.call(rbind,by(pmhcs, pmhcs$numMut, function(x){
  return(data.frame(numMut=as.numeric(x[["numMut"]][1]), count=nrow(x)))
}))

(p <- ggplot(pmhcs_nmut, aes(x=numMut, y=count)) + geom_point() + scale_y_log10() + xlab(expression(paste("pMHCs derived from ", italic("n"), " mutations"))) + ylab("Count (log scale)") + theme_bw())
ggsave("figure4/numMutPer_pMHC_180313.pdf",p,width=2.17,height=3.05, useDingbats=FALSE)



#### Figure 5B ####

## only recurrent
rpmhcs$label <- ""
rpmhcs$label[rpmhcs$pMHC=="KIGDFGLATEK_HLA-A*03:01"] <- "KIGDFGLATEK_HLA-A*03:01_BRAF-V600E"
rpmhcs$label[rpmhcs$pMHC=="KLVVVGADGV_HLA-A*02:01"] <- "KLVVVGADGV_HLA-A*02:01_KRAS/NRAS/HRAS-G12D"
rpmhcs$label[rpmhcs$pMHC=="LVVVGAVGV_HLA-A*02:01"] <- "LVVVGAVGV_HLA-A*02:01_KRAS/HRAS-G12V"
rpmhcs$label[rpmhcs$pMHC=="KLVVVGAVGV_HLA-A*02:01"] <- "KLVVVGAVGV_HLA-A*02:01_KRAS/HRAS-G12V"
rpmhcs$label[rpmhcs$pMHC=="KIGDFGLATEK_HLA-A*11:01"] <- "KIGDFGLATEK_HLA-A*11:01_BRAF-V600E"
rpmhcs$label[rpmhcs$pMHC=="SEITKQEKDF_HLA-B*44:02"] <- "SEITKQEKDF_HLA-B*44:02_PIK3CA-E545K"
rpmhcs$label[rpmhcs$pMHC=="KLVVVGACGV_HLA-A*02:01"] <- "KLVVVGACGV_HLA-A*02:01_KRAS/NRAS/HRAS-G12C"
rpmhcs$label[rpmhcs$pMHC=="VVVGAVGVGK_HLA-A*03:01"] <- "VVVGAVGVGK_HLA-A*03:01_KRAS-G12V"
rpmhcs$label[rpmhcs$pMHC=="VVGAVGVGK_HLA-A*03:01"] <- "VVGAVGVGK_HLA-A*03:01_KRAS-G12V"
rpmhcs$label[rpmhcs$pMHC=="KPIIIGHHAY_HLA-B*15:01"] <- "KPIIIGHHAY_HLA-B*15:01_IDH1-R132H"
rpmhcs$label[rpmhcs$pMHC=="AISTRDPLSK_HLA-A*03:01"] <- "AISTRDPLSK_HLA-A*03:01_PIK3CA-E542K"
rpmhcs$label[rpmhcs$pMHC=="KPIIIGHHAY_HLA-B*35:01"] <- "KPIIIGHHAY_HLA-B*35:01_IDH1-R132H"
rpmhcs$label[rpmhcs$pMHC=="VVGADGVGK_HLA-A*11:01"] <- "VVGADGVGK_HLA-A*11:01_KRAS-G12D"
rpmhcs$label[rpmhcs$pMHC=="STRDPLSEITK_HLA-A*11:01"] <- "STRDPLSEITK_HLA-A*11:01_PIK3CA-E545K"
rpmhcs$label[rpmhcs$pMHC=="VVVGADGVGK_HLA-A*11:01"] <- "VVVGADGVGK_HLA-A*11:01_KRAS-G12D"
rpmhcs$label[rpmhcs$pMHC=="KLVVVGAGDV_HLA-A*02:01"] <- "KLVVVGAGDV_HLA-A*02:01_KRAS-G13D"
rpmhcs$label[rpmhcs$pMHC=="GMNWRPILTI_HLA-A*02:01"] <- "GMNWRPILTI_HLA-A*02:01_TP53-R248W"
rpmhcs$label[rpmhcs$pMHC=="FGLATEKSRW_HLA-B*57:01"] <- "FGLATEKSRW_HLA-B*57:01_BRAF-V600E"
rpmhcs$label[rpmhcs$pMHC=="CMGGMNWRPI_HLA-A*02:01"] <- "CMGGMNWRPI_HLA-A*02:01_TP53-R248W"
rpmhcs$label[rpmhcs$pMHC=="ILDTAGREEY_HLA-A*01:01"] <- "ILDTAGREEY_HLA-A*01:01_KRAS/NRAS/HRAS-Q61R"
rpmhcs$label[rpmhcs$pMHC=="LVVVGAAGV_HLA-A*02:01"] <- "LVVVGAAGV_HLA-A*02:01_KRAS-G12A"
rpmhcs$label[rpmhcs$pMHC=="KLVVVGAAGV_HLA-A*02:01"] <- "KLVVVGAAGV_HLA-A*02:01_KRAS-G12A"



p <- ggplot(rpmhcs, aes(x=reorder(pMHC, -numSub), y=numSub, size=numMut)) + geom_point(shape=1) + geom_text_repel(data = subset(rpmhcs, numSub > 30), color="black", size=2, nudge_x=12000, segment.alpha=0.5, segment.color="grey", force=2.5, aes(label=label)) + scale_size_continuous(breaks=c(0,1,5,10,100), labels=c("0","1","2-5","6-10",">10"), name = "Number of \nmutations") + scale_x_discrete(expand=c(0.02,0)) + xlab("Peptide-MHCs ordered by frequency") + ylab("Number of Subjects") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave("figure4/recurrentPeptideMHCsOnly_hockeystick_170717.pdf", p, width=6, height=3.05, useDingbats=FALSE)




#### Supplemetal Table - Recurrent pMHCs ####

write.table(rpmhcs[order(rpmhcs$numSub, decreasing = T),], "ST2_RecurrentPMHCs.tsv", quote = F, sep = "\t", row.names = F)
