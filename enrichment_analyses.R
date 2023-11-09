####################### determining probe overlap for each platform
load("/dartfs/rc/lab/S/SalasLab/R/Annotation files/annotationEPICb5updated3.rda")
rm(annotDF)
annot$probeType <- substr(annot$Name, 1, 2)
EPICv1_probes <- rownames(annot)

HM450.hg38.manifest <- readRDS("/dartfs/rc/lab/S/SalasLab/Annotation files/HM450.hg38.manifest.rds")
HM450k_probes <- HM450.hg38.manifest@ranges@NAMES

EPICv2 <- read.csv("/dartfs/rc/lab/S/SalasLab/EPICv2/EPIC_v2_annotation.csv")
EPICv2_probes <- EPICv2$Name

setwd("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/")
for (tt in c("GBM","OLG","AST_LG","AST_HG")){
  load(paste0(tt,"_Libraries.RDATA"))
  lib_names <- ls()[grepl("Layer",ls())]
  tmp_libs <- list(Library_Layer0,Library_Layer1,Library_Layer2A,
                  Library_Layer2B,Library_Layer2C,Library_Layer2D,
                  Library_Layer3A,Library_Layer3B,Library_Layer4,
                  Library_Layer5A,Library_Layer5B,Library_Layer5C)
  names(tmp_libs) <- lib_names
  
  summ_df <- data.frame(matrix(NA, nrow=length(tmp_libs), ncol=5))
  colnames(summ_df) <- c("Layer","n_layer","n_450k","n_EPIC","n_EPICv2")
  summ_df$Layer <- lib_names
  rownames(summ_df) <- lib_names
  first <- TRUE
  for (lib_name in lib_names){
    tmp_lib <- tmp_libs[[lib_name]]
    summ_df[lib_name, "n_layer"] <- nrow(tmp_lib)
    summ_df[lib_name, "n_450k"] <- length(intersect(HM450k_probes,rownames(tmp_lib)))
    summ_df[lib_name, "n_EPIC"] <- length(intersect(EPICv1_probes,rownames(tmp_lib)))
    summ_df[lib_name, "n_EPICv2"] <- length(intersect(EPICv2_probes,rownames(tmp_lib)))
  
    if (first==TRUE){
      first <- FALSE
      intersect_cpgs <- intersect(EPICv2_probes,rownames(tmp_lib))
    } else {
      intersect_cpgs <- unique(intersect_cpgs,intersect(EPICv2_probes,rownames(tmp_lib)))
      
    }
  }
  
  write.csv(summ_df, file=paste0(tt,"_platform_overlap.csv"),row.names = FALSE)
  save(intersect_cpgs, file=paste0(tt,"_platformcommon_cpgs.RDA"))
}




####################### Cell-type specific DNAm
library(missMethyl)
library(IlluminaHumanMethylation450kmanifest)
library(lmerTest)
annot <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

#### Saving annotated CpG information
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/")
load("GBM_Libraries.RDATA")
rm(Library_Layer0)
tmp_libs <- list(Library_Layer1,Library_Layer2A,
                 Library_Layer2B,Library_Layer2C,Library_Layer2D,
                 Library_Layer3A,Library_Layer3B,Library_Layer4,
                 Library_Layer5A,Library_Layer5B,Library_Layer5C)
lib_names <- c("Library_Layer1","Library_Layer2A",
                 "Library_Layer2B","Library_Layer2C","Library_Layer2D",
                 "Library_Layer3A","Library_Layer3B","Library_Layer4",
                 "Library_Layer5A","Library_Layer5B","Library_Layer5C")
names(tmp_libs) <- lib_names
cell_cpgs <- c()
for (lib_name in lib_names){
  tmp_lib <- tmp_libs[[lib_name]]
  cell_cpgs <- c(cell_cpgs,rownames(tmp_lib))
}
cell_cpgs <- unique(cell_cpgs)

annot$DMC <- ifelse(annot$Name %in% cell_cpgs,1,0)
annot_DMC <- annot[annot$DMC == 1,]
write.csv(annot_DMC, file="CellLayers_DMCinfo.csv",row.names=FALSE)

#### Gene Ontology
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/")
tmp_results <- gometh(sig.cpg=cell_cpgs,
                        collection="GO",
                        array.type="450k")
save(tmp_results, file = "CellLayers_GO.RDA")

top10 <- top10[order(top10$FDR, decreaseing=FALSE),]
top10 <- top10[c(1:10),]

#### Chromosome
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/")
loci_test <- unique(annot$chr)
loci_test <- loci_test[-which(loci_test %in% c("chrY","chrX"))]
loci_summary <- data.frame(matrix(NA, nrow=length(loci_test), ncol=7))
colnames(loci_summary) <- c("Location","Coef","OR","SE","pval","OR_CI_hi","OR_CI_lo")
rownames(loci_summary) <- loci_test
loci_summary$Location <- loci_test
for (loci_name in loci_test){
  print(loci_name)
  annot$tmp <- ifelse(annot$chr == loci_name, 1, 0)
  tmp_model <- glm(tmp~DMC, data=annot, family = binomial)
  tmp_summary <- summary(tmp_model)
  tmp_summary <- tmp_summary$coefficients
  tmp_confint <- confint(tmp_model)
  loci_summary[loci_name, "Coef"] <- tmp_summary[2,1]
  loci_summary[loci_name, "SE"] <- tmp_summary[2,2]
  loci_summary[loci_name, "pval"] <- tmp_summary[2,4]
  loci_summary[loci_name, "OR"] <- exp(loci_summary[loci_name, "Coef"])
  loci_summary[loci_name, "OR_CI_hi"] <- exp(tmp_confint[2,2])
  loci_summary[loci_name, "OR_CI_lo"] <- exp(tmp_confint[2,1])
}

write.csv(loci_summary, file = "CellLayers_enrich_chr.csv")

#### Location
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/")
loci_test <- c("5'UTR", "3'UTR","1stExon","Body","TSS200","TSS1500")
loci_summary <- data.frame(matrix(NA, nrow=length(loci_test), ncol=7))
colnames(loci_summary) <- c("Location","Coef","OR","SE","pval","OR_CI_hi","OR_CI_lo")
rownames(loci_summary) <- loci_test
loci_summary$Location <- loci_test
for (loci_name in loci_test){
  print(loci_name)
  annot$tmp <- 0
  annot$tmp[grep(loci_name,annot$UCSC_RefGene_Group)] <- 1
  tmp_model <- glm(tmp~DMC, data=annot, family = binomial)
  tmp_summary <- summary(tmp_model)
  tmp_summary <- tmp_summary$coefficients
  tmp_confint <- confint(tmp_model)
  loci_summary[loci_name, "Coef"] <- tmp_summary[2,1]
  loci_summary[loci_name, "SE"] <- tmp_summary[2,2]
  loci_summary[loci_name, "pval"] <- tmp_summary[2,4]
  loci_summary[loci_name, "OR"] <- exp(loci_summary[loci_name, "Coef"])
  loci_summary[loci_name, "OR_CI_hi"] <- exp(tmp_confint[2,2])
  loci_summary[loci_name, "OR_CI_lo"] <- exp(tmp_confint[2,1])
}

write.csv(loci_summary, file = "CellLayers_enrich_location.csv")


#### Relation to CpG Isalnd
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/")
loci_test <- unique(annot$Relation_to_Island)
loci_summary <- data.frame(matrix(NA, nrow=length(loci_test), ncol=7))
colnames(loci_summary) <- c("Location","Coef","OR","SE","pval","OR_CI_hi","OR_CI_lo")
rownames(loci_summary) <- loci_test
loci_summary$Location <- loci_test
for (loci_name in loci_test){
  print(loci_name)
  annot$tmp <- ifelse(annot$Relation_to_Island == loci_name, 1, 0)
  tmp_model <- glm(tmp~DMC, data=annot, family = binomial)
  tmp_summary <- summary(tmp_model)
  tmp_summary <- tmp_summary$coefficients
  tmp_confint <- confint(tmp_model)
  loci_summary[loci_name, "Coef"] <- tmp_summary[2,1]
  loci_summary[loci_name, "SE"] <- tmp_summary[2,2]
  loci_summary[loci_name, "pval"] <- tmp_summary[2,4]
  loci_summary[loci_name, "OR"] <- exp(loci_summary[loci_name, "Coef"])
  loci_summary[loci_name, "OR_CI_hi"] <- exp(tmp_confint[2,2])
  loci_summary[loci_name, "OR_CI_lo"] <- exp(tmp_confint[2,1])
}

write.csv(loci_summary, file = "CellLayers_enrich_reltoisland.csv")

####################### L0 Libraries
#### Saving annotated CpG information
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/")
for (tt in c("GBM","OLG","AST_LG","AST_HG")){
  load(paste0(tt,"_Libraries.RDATA"))
  annot$DMC <- ifelse(annot$Name %in% rownames(Library_Layer0),1,0)
  annot_DMC <- annot[annot$DMC == 1,]
  write.csv(annot_DMC, file=paste0(tt,"_L0_DMCinfo.csv"),row.names=FALSE)
}

#### Gene Ontology
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/")
for (tt in c("GBM","OLG","AST_LG","AST_HG")){
  print(tt)
  load(paste0(tt,"_Libraries.RDATA"))
  tmp_results <- gometh(sig.cpg=rownames(Library_Layer0),
                        collection="GO",
                        array.type="450k")
  save(tmp_results, file = paste0(tt,"_L0_GO.RDA"))
}

#### Chromosome
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/")
for (tt in c("GBM","OLG","AST_LG","AST_HG")){
  print(tt)
  load(paste0(tt,"_Libraries.RDATA"))
  annot$DMC <- ifelse(annot$Name %in% rownames(Library_Layer0),1,0)
  loci_test <- unique(annot$chr)
  loci_test <- loci_test[-which(loci_test %in% c("chrY","chrX"))]
  loci_summary <- data.frame(matrix(NA, nrow=length(loci_test), ncol=7))
  colnames(loci_summary) <- c("Location","Coef","OR","SE","pval","OR_CI_hi","OR_CI_lo")
  rownames(loci_summary) <- loci_test
  loci_summary$Location <- loci_test
  for (loci_name in loci_test){
    print(loci_name)
    annot$tmp <- ifelse(annot$chr == loci_name, 1, 0)
    tmp_model <- glm(tmp~DMC, data=annot, family = binomial)
    tmp_summary <- summary(tmp_model)
    tmp_summary <- tmp_summary$coefficients
    tmp_confint <- confint(tmp_model)
    loci_summary[loci_name, "Coef"] <- tmp_summary[2,1]
    loci_summary[loci_name, "SE"] <- tmp_summary[2,2]
    loci_summary[loci_name, "pval"] <- tmp_summary[2,4]
    loci_summary[loci_name, "OR"] <- exp(loci_summary[loci_name, "Coef"])
    loci_summary[loci_name, "OR_CI_hi"] <- exp(tmp_confint[2,2])
    loci_summary[loci_name, "OR_CI_lo"] <- exp(tmp_confint[2,1])
  }
  write.csv(loci_summary, file = paste0(tt,"_L0_enrich_chr.csv"))
}

#### Location
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/")
loci_test <- c("5'UTR", "3'UTR","1stExon","Body","TSS200","TSS1500")
for (tt in c("GBM","OLG","AST_LG","AST_HG")){
  print(tt)
  load(paste0(tt,"_Libraries.RDATA"))
  annot$DMC <- ifelse(annot$Name %in% rownames(Library_Layer0),1,0)
  loci_summary <- data.frame(matrix(NA, nrow=length(loci_test), ncol=7))
  colnames(loci_summary) <- c("Location","Coef","OR","SE","pval","OR_CI_hi","OR_CI_lo")
  rownames(loci_summary) <- loci_test
  loci_summary$Location <- loci_test
  for (loci_name in loci_test){
    print(loci_name)
    annot$tmp <- 0
    annot$tmp[grep(loci_name,annot$UCSC_RefGene_Group)] <- 1
    tmp_model <- glm(tmp~DMC, data=annot, family = binomial)
    tmp_summary <- summary(tmp_model)
    tmp_summary <- tmp_summary$coefficients
    tmp_confint <- confint(tmp_model)
    loci_summary[loci_name, "Coef"] <- tmp_summary[2,1]
    loci_summary[loci_name, "SE"] <- tmp_summary[2,2]
    loci_summary[loci_name, "pval"] <- tmp_summary[2,4]
    loci_summary[loci_name, "OR"] <- exp(loci_summary[loci_name, "Coef"])
    loci_summary[loci_name, "OR_CI_hi"] <- exp(tmp_confint[2,2])
    loci_summary[loci_name, "OR_CI_lo"] <- exp(tmp_confint[2,1])
  }
  write.csv(loci_summary, file = paste0(tt,"_L0_enrich_location.csv"))
}

#### Relation to CpG Isalnd
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/")
for (tt in c("GBM","OLG","AST_LG","AST_HG")){
  load(paste0(tt,"_Libraries.RDATA"))
  print(tt)
  annot$DMC <- ifelse(annot$Name %in% rownames(Library_Layer0),1,0)
  loci_test <- unique(annot$Relation_to_Island)
  loci_summary <- data.frame(matrix(NA, nrow=length(loci_test), ncol=7))
  colnames(loci_summary) <- c("Location","Coef","OR","SE","pval","OR_CI_hi","OR_CI_lo")
  rownames(loci_summary) <- loci_test
  loci_summary$Location <- loci_test
  for (loci_name in loci_test){
    print(loci_name)
    annot$tmp <- ifelse(annot$Relation_to_Island == loci_name, 1, 0)
    tmp_model <- glm(tmp~DMC, data=annot, family = binomial)
    tmp_summary <- summary(tmp_model)
    tmp_summary <- tmp_summary$coefficients
    tmp_confint <- confint(tmp_model)
    loci_summary[loci_name, "Coef"] <- tmp_summary[2,1]
    loci_summary[loci_name, "SE"] <- tmp_summary[2,2]
    loci_summary[loci_name, "pval"] <- tmp_summary[2,4]
    loci_summary[loci_name, "OR"] <- exp(loci_summary[loci_name, "Coef"])
    loci_summary[loci_name, "OR_CI_hi"] <- exp(tmp_confint[2,2])
    loci_summary[loci_name, "OR_CI_lo"] <- exp(tmp_confint[2,1])
  }
  write.csv(loci_summary, file = paste0(tt,"_L0_enrich_reltoisland.csv"))
}


############### Plotting
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS")
df <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS/GBM_LO_DMCs_GeneRegion.csv")
library(ggplot2)
plot1 <- ggplot(df, aes(y = Location, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = OR_CI_lo, xmax = OR_CI_hi), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("Odds Ratio (95% CI)") + 
  ylab("Gene Region") + 
  theme_bw() 
pdf("GBM_L0_GeneRegion_Forestplot.pdf", height = 3, width = 4)
plot1
dev.off()

setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS")
df <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS/IDH_LO_DMCs_GeneRegion.csv")
library(ggplot2)
plot1 <- ggplot(df, aes(y = Location, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = OR_CI_lo, xmax = OR_CI_hi), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("Odds Ratio (95% CI)") + 
  ylab("Gene Region") + 
  theme_bw() 
pdf("IDH_L0_GeneRegion_Forestplot.pdf", height = 3, width = 4)
plot1
dev.off()


setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS")
df <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS/Layers_DMCs_GeneRegion.csv")
library(ggplot2)
plot1 <- ggplot(df, aes(y = Location, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = OR_CI_lo, xmax = OR_CI_hi), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("Odds Ratio (95% CI)") + 
  ylab("Gene Region") + 
  theme_bw() 
pdf("Layers_GeneRegion_Forestplot.pdf", height = 3, width = 4)
plot1
dev.off()


setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS")
df <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS/GBM_L0_DMCs_Relation_to_Island.csv")
library(ggplot2)
plot1 <- ggplot(df, aes(y = Location, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = OR_CI_lo, xmax = OR_CI_hi), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("Odds Ratio (95% CI)") + 
  ylab("Relaion to CpG Island") + 
  theme_bw() 
pdf("GBM_L0_RelationtoIsland_Forestplot.pdf", height = 3, width = 4)
plot1
dev.off()

setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS")
df <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS/IDH_L0_DMCs_Relation_to_Island.csv")
library(ggplot2)
plot1 <- ggplot(df, aes(y = Location, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = OR_CI_lo, xmax = OR_CI_hi), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("Odds Ratio (95% CI)") + 
  ylab("Relaion to CpG Island") + 
  theme_bw() 
pdf("IDH_L0_RelationtoIsland_Forestplot.pdf", height = 3, width = 4)
plot1
dev.off()



setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS")
df <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS/Layers_DMCs_Relation_to_Island.csv")
library(ggplot2)
plot1 <- ggplot(df, aes(y = Location, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = OR_CI_lo, xmax = OR_CI_hi), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("Odds Ratio (95% CI)") + 
  ylab("Relaion to CpG Island") + 
  theme_bw() 
pdf("Layers_RelationtoIsland_Forestplot.pdf", height = 3, width = 4)
plot1
dev.off()




setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS")
df <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS/GBM_hiinfilEWAS_DMCs_Relation_to_Island.csv")
library(ggplot2)
plot1 <- ggplot(df, aes(y = Location, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = OR_CI_lo, xmax = OR_CI_hi), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("Odds Ratio (95% CI)") + 
  ylab("Relaion to CpG Island") + 
  theme_bw() 
pdf("GBM_infilEWAS_RelationtoIsland_Forestplot.pdf", height = 3, width = 4)
plot1
dev.off()



setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS")
df <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS/GBM_hiinfilEWAS_DMCs_Relation_to_Island.csv")
library(ggplot2)
plot1 <- ggplot(df, aes(y = Location, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = OR_CI_lo, xmax = OR_CI_hi), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("Odds Ratio (95% CI)") + 
  ylab("Relaion to CpG Island") + 
  theme_bw() 
pdf("GBM_infilEWAS_RelationtoIsland_Forestplot.pdf", height = 3, width = 4)
plot1
dev.off()



setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS")
df <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS/IDH_hiinfilEWAS_DMCs_Relation_to_Island.csv")
library(ggplot2)
plot1 <- ggplot(df, aes(y = Location, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = OR_CI_lo, xmax = OR_CI_hi), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("Odds Ratio (95% CI)") + 
  ylab("Relaion to CpG Island") + 
  theme_bw() 
pdf("IDH_infilEWAS_RelationtoIsland_Forestplot.pdf", height = 3, width = 4)
plot1
dev.off()




setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS")
df <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS/GBM_hiinfilEWAS_DMCs_GeneRegion.csv")
library(ggplot2)
plot1 <- ggplot(df, aes(y = Location, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = OR_CI_lo, xmax = OR_CI_hi), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("Odds Ratio (95% CI)") + 
  ylab("Relaion to CpG Island") + 
  theme_bw() 
pdf("GBM_infilEWAS_GeneRegion_Forestplot.pdf", height = 3, width = 4)
plot1
dev.off()




setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS")
df <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/CapperEWAS/IDH_hiinfilEWAS_DMCs_GeneRegion.csv")
library(ggplot2)
plot1 <- ggplot(df, aes(y = Location, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = OR_CI_lo, xmax = OR_CI_hi), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("Odds Ratio (95% CI)") + 
  ylab("Relaion to CpG Island") + 
  theme_bw() 
pdf("IDH_infilEWAS_GeneRegion_Forestplot.pdf", height = 3, width = 4)
plot1
dev.off()