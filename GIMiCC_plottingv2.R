############ TCGA info
library(readxl)
library(dplyr)
library(readr)
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Brain_Tumor_Deconv_Ref/Processed_data/TCGA_ref/TCGAall_pheno.RDA")
TCGA_Supp_Info <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/TCGA_Supp_Info.xlsx")
TCGA_pheno <- left_join(TCGA_pheno, TCGA_Supp_Info, by="submitter_id")
TCGA_pheno$age_at_index <- as.numeric(TCGA_pheno$age_at_index)

df <- read_tsv(file = "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/TCGA_extra_info.tsv")
df <- df[c(5,6,7),]
df <- t(df)
colnames(df) <- df[1,]
df <- data.frame(df[-c(1,2),])
df$submitter_id <- rownames(df)
TCGA_pheno <- left_join(TCGA_pheno, df, by="submitter_id")

TCGA_pheno$Tumor_Type2 <- ifelse(TCGA_pheno$`IDH/codel subtype.y` == "IDHwt", "GBM",
                                 ifelse(TCGA_pheno$`IDH/codel subtype.y`== "IDHmut-codel", "OLG", 
                                        ifelse(TCGA_pheno$`IDH/codel subtype.y`== "IDHmut-non-codel", "AST", "N/A")))

TCGA_pheno$`TERT promoter status.y` <- ifelse(is.na(TCGA_pheno$`TERT promoter status.y`), "N/A", TCGA_pheno$`TERT promoter status.y`)
TCGA_pheno$`Chr 7 gain/Chr 10 loss` <- ifelse(is.na(TCGA_pheno$`Chr 7 gain/Chr 10 loss`), "N/A", TCGA_pheno$`TERT promoter status.y`)
TCGA_pheno$EGFR <- ifelse(is.na(TCGA_pheno$EGFR), "N/A", TCGA_pheno$EGFR)
TCGA_pheno$CDKN2A <- ifelse(is.na(TCGA_pheno$CDKN2A), "N/A", TCGA_pheno$CDKN2A)
TCGA_pheno$CDKN2B <- ifelse(is.na(TCGA_pheno$CDKN2B), "N/A", TCGA_pheno$CDKN2B)

TCGA_pheno$Tumor_Type2 <- ifelse(TCGA_pheno$Tumor_Type2 != "AST", TCGA_pheno$Tumor_Type2, 
                                 ifelse(TCGA_pheno$`TERT promoter status.y`== "Mutant" | 
                                          TCGA_pheno$`Chr 7 gain/Chr 10 loss`== "Gain chr 7 & loss chr 10" |
                                          TCGA_pheno$EGFR == "amp_rec" |
                                          TCGA_pheno$CDKN2A== "homdel_rec" |
                                          TCGA_pheno$CDKN2B== "homdel_rec" , "AST-HG", "AST-LG" ))

TCGA_pheno$EGFRv2 <- ifelse(TCGA_pheno$EGFR == "amp_rec", "Yes", "Not Detected")
TCGA_pheno$CDKN2Av2 <- ifelse(TCGA_pheno$CDKN2A == "homdel_rec", "Yes", "Not Detected")
TCGA_pheno$CDKN2Bv2 <- ifelse(TCGA_pheno$CDKN2B == "homdel_rec", "Yes", "Not Detected")

TCGA_pheno$IDH <- ifelse(TCGA_pheno$`IDH status.y` == "Mutant", "Yes", "No")
TCGA_pheno$co1p19q <- ifelse(TCGA_pheno$`1p/19q codeletion.y` == "codel", "Yes", "No")


library(tableone)
myVars <- c("age_at_index","gender","race","IDH","co1p19q","TERT promoter status.y",
            "Chr 7 gain/Chr 10 loss","EGFRv2","CDKN2Av2","CDKN2Bv2")
catVars <- c("gender","race","IDH","co1p19q","TERT promoter status.y",
             "Chr 7 gain/Chr 10 loss","EGFRv2","CDKN2Av2","CDKN2Bv2")

TCGA_pheno <- TCGA_pheno[TCGA_pheno$Tumor_Type2 != "N/A",]
tab2 <- CreateTableOne(vars = myVars, data = TCGA_pheno, factorVars = catVars, strata = "Tumor_Type2")

tab2 <- print(tab2)

setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/")
write.csv(tab2, file="TCGA_table1.csv")

####################### gene ontonolgy csv files
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/AST_HG_L0_GO.RDA")
tmp_results <- tmp_results[tmp_results$FDR <= 0.05,]
tmp_results$`GO Term` <- rownames(tmp_results)
tmp_results$Library <- "AST-HG L0"
final_results <- tmp_results

########## no results
# load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/AST_LG_L0_GO.RDA")
# tmp_results <- tmp_results[tmp_results$FDR <= 0.05,]
# tmp_results$`GO Term` <- rownames(tmp_results)
# tmp_results$Library <- "AST-LG L0"
# final_results <- rbind(final_results,tmp_results )

load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/OLG_L0_GO.RDA")
tmp_results <- tmp_results[tmp_results$FDR <= 0.05,]
tmp_results$`GO Term` <- rownames(tmp_results)
tmp_results$Library <- "OLG L0"
final_results <- rbind(final_results,tmp_results )


load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/GBM_L0_GO.RDA")
tmp_results <- tmp_results[tmp_results$FDR <= 0.05,]
tmp_results$`GO Term` <- rownames(tmp_results)
tmp_results$Library <- "GBM L0"
final_results <- rbind(final_results,tmp_results )


load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/CellLayers_GO.RDA")
tmp_results <- tmp_results[tmp_results$FDR <= 0.05,]
tmp_results$`GO Term` <- rownames(tmp_results)
tmp_results$Library <- "L1-L5"
final_results <- rbind(final_results,tmp_results )

write.csv(final_results,
          file = "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/GO_Summary.csv", row.names = FALSE)
###################################### Figure 2A
############### go into cluster and get betas 
library(limma)
library(dplyr)
library(ggplot2)
library(Metrics)
library(tibble)
library(openxlsx)
library(tidyr)
library(readr)
library(readxl)
library(dplyr)
load("/dartfs/rc/lab/S/SalasLab/Capper/Processed_data/Capper_pheno_agg2.RDA")
Capper_Supp_Info <- read_excel("/dartfs//rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Capper_Supp_Info.xlsx")
rownames(Capper_pheno) <- paste0(Capper_pheno$Slide,"_",Capper_pheno$Array)
Capper_pheno$ID <- substr(rownames(Capper_pheno), 12, 100)
pheno_df <- left_join(Capper_pheno, Capper_Supp_Info)
rownames(pheno_df) <- paste0(pheno_df$Slide,"_",pheno_df$Array)

load("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/GBM_Library0.RDA")
GBM_Layer0 <- Library_Layer0
load("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/OLI_Library0.RDA")
OLI_Layer0 <- Library_Layer0
load("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/AST_HG_Library0.RDA")
AST_HG_Layer0 <- Library_Layer0
load("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/AST_LG_Library0.RDA")
AST_LG_Layer0 <- Library_Layer0

tot_cpgs <- unique(c(rownames(AST_HG_Layer0), rownames(AST_LG_Layer0),
                     rownames(GBM_Layer0),rownames(OLI_Layer0)))


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

HC_samps <- rownames(pheno_df[pheno_df$Broad_Class == "Control_Healthy",])
GBM_samps <- rownames(pheno_df[pheno_df$Broad_Class == "Glioblastoma",])
OLG_samps <- rownames(pheno_df[pheno_df$methylation.class.ch1 == "O IDH",])
AST_LG_samps <- rownames(pheno_df[pheno_df$methylation.class.ch1 == "A IDH",])
AST_HG_samps <- rownames(pheno_df[pheno_df$methylation.class.ch1 == "A IDH, HG",])

tot_samps <- c(HC_samps, GBM_samps, OLG_samps, AST_LG_samps, AST_HG_samps)
Capper_GBM_subset <- pheno_df[which(rownames(pheno_df) %in% tot_samps),]

base_dir <- "/dartfs/rc/lab/S/SalasLab/Capper/Processed_data/"
first <- TRUE
for (i in unique(Capper_GBM_subset$processing_batch)){
  print(i)
  set_name <- paste0(base_dir,"Capper_batch",i,"/Capper_batch",i,"_betas_filtered.RDA")
  betas <- loadRData(set_name)
  betas <- betas[, which(colnames(betas) %in% rownames(Capper_GBM_subset))]
  if (first == TRUE){
    Capper_GBM_betas <- betas
    first <- FALSE
  } else {
    if (is.null(dim(betas))){
      CpGs <- tot_cpgs
      betas <- betas[CpGs]
      betas <- betas[order(names(betas))]
      Capper_GBM_betas <- Capper_GBM_betas[CpGs,]
      Capper_GBM_betas <- Capper_GBM_betas[order(rownames(Capper_GBM_betas)),]
      print(identical(names(betas), rownames(Capper_GBM_betas)))
    } else {
      CpGs <- tot_cpgs
      betas <- betas[CpGs,]
      betas <- betas[order(rownames(betas)),]
      Capper_GBM_betas <- Capper_GBM_betas[CpGs,]
      Capper_GBM_betas <- Capper_GBM_betas[order(rownames(Capper_GBM_betas)),]
      print(identical(rownames(betas), rownames(Capper_GBM_betas)))
    }
    Capper_GBM_betas <- cbind(Capper_GBM_betas, betas)
  }
}

save(Capper_GBM_betas, file="/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/tot_L0_cpgs.RDA")


# bring in local for plotting
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/tot_L0_cpgs.RDA")
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/GBM_Library0.RDA")
GBM_Layer0 <- Library_Layer0
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/OLI_Library0.RDA")
OLI_Layer0 <- Library_Layer0
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/AST_HG_Library0.RDA")
AST_HG_Layer0 <- Library_Layer0
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/AST_LG_Library0.RDA")
AST_LG_Layer0 <- Library_Layer0

table(GBM_Layer0$hyper)
table(OLI_Layer0$hyper)
table(AST_HG_Layer0$hyper)
table(AST_LG_Layer0$hyper)


library(readr)
library(readxl)
library(dplyr)
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Capper/Processed_data/Capper_pheno_agg2.RDA")
Capper_Supp_Info <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Capper_Supp_Info.xlsx")
rownames(Capper_pheno) <- paste0(Capper_pheno$Slide,"_",Capper_pheno$Array)
Capper_pheno$ID <- substr(rownames(Capper_pheno), 12, 100)
pheno_df <- left_join(Capper_pheno, Capper_Supp_Info)
rownames(pheno_df) <- paste0(pheno_df$Slide,"_",pheno_df$Array)

HC_samps <- rownames(pheno_df[pheno_df$Broad_Class == "Control_Healthy",])
GBM_samps <- rownames(pheno_df[pheno_df$Broad_Class == "Glioblastoma",])
OLG_samps <- rownames(pheno_df[pheno_df$methylation.class.ch1 == "O IDH",])
AST_LG_samps <- rownames(pheno_df[pheno_df$methylation.class.ch1 == "A IDH",])
AST_HG_samps <- rownames(pheno_df[pheno_df$methylation.class.ch1 == "A IDH, HG",])

train_samps <- rownames(pheno_df[pheno_df$Train_split == TRUE,])
test_samps <- rownames(pheno_df[pheno_df$Train_split == FALSE,])

HC_samps_test <- intersect(HC_samps, test_samps)
GBM_samps_test <- intersect(GBM_samps, test_samps)
OLG_samps_test <- intersect(OLG_samps, test_samps)
AST_LG_samps_test <- intersect(AST_LG_samps, test_samps)
AST_HG_samps_test <- intersect(AST_HG_samps, test_samps)

HC_samps_train <- intersect(HC_samps, train_samps)
GBM_samps_train <- intersect(GBM_samps, train_samps)
OLG_samps_train <- intersect(OLG_samps, train_samps)
AST_LG_samps_train <- intersect(AST_LG_samps, train_samps)
AST_HG_samps_train <- intersect(AST_HG_samps, train_samps)

summ_betas <- data.frame(matrix(nrow=nrow(Capper_GBM_betas), ncol=5))
colnames(summ_betas) <- c("HC", "GBM", "OLG", "AST-LG", "AST-HG")
rownames(summ_betas) <- rownames(Capper_GBM_betas)
summ_betas$GBM_L0 <- ifelse(rownames(summ_betas) %in% rownames(GBM_Layer0), "Yes", "No")
summ_betas$OLG_L0 <- ifelse(rownames(summ_betas) %in% rownames(OLI_Layer0),  "Yes", "No")
summ_betas$AST_HG_L0 <- ifelse(rownames(summ_betas) %in% rownames(AST_HG_Layer0),  "Yes", "No")
summ_betas$AST_LG_L0 <- ifelse(rownames(summ_betas) %in% rownames(AST_LG_Layer0),  "Yes", "No")
summ_betas_test <- summ_betas
summ_betas_train <- summ_betas
for (cpg in rownames(summ_betas)){
  summ_betas[cpg, "HC"] <- mean(Capper_GBM_betas[cpg, which(colnames(Capper_GBM_betas) %in% HC_samps)], na.rm=TRUE)
  summ_betas[cpg, "GBM"] <- mean(Capper_GBM_betas[cpg, which(colnames(Capper_GBM_betas) %in% GBM_samps)], na.rm=TRUE)
  summ_betas[cpg, "OLG"] <- mean(Capper_GBM_betas[cpg, which(colnames(Capper_GBM_betas) %in% OLG_samps)], na.rm=TRUE)
  summ_betas[cpg, "AST-LG"] <- mean(Capper_GBM_betas[cpg, which(colnames(Capper_GBM_betas) %in% AST_LG_samps)], na.rm=TRUE)
  summ_betas[cpg, "AST-HG"] <- mean(Capper_GBM_betas[cpg, which(colnames(Capper_GBM_betas) %in% AST_HG_samps)], na.rm=TRUE)

  summ_betas_test[cpg, "HC"] <- mean(Capper_GBM_betas[cpg, which(colnames(Capper_GBM_betas) %in% HC_samps_test)], na.rm=TRUE)
  summ_betas_test[cpg, "GBM"] <- mean(Capper_GBM_betas[cpg, which(colnames(Capper_GBM_betas) %in% GBM_samps_test)], na.rm=TRUE)
  summ_betas_test[cpg, "OLG"] <- mean(Capper_GBM_betas[cpg, which(colnames(Capper_GBM_betas) %in% OLG_samps_test)], na.rm=TRUE)
  summ_betas_test[cpg, "AST-LG"] <- mean(Capper_GBM_betas[cpg, which(colnames(Capper_GBM_betas) %in% AST_LG_samps_test)], na.rm=TRUE)
  summ_betas_test[cpg, "AST-HG"] <- mean(Capper_GBM_betas[cpg, which(colnames(Capper_GBM_betas) %in% AST_HG_samps_test)], na.rm=TRUE)
  
  summ_betas_train[cpg, "HC"] <- mean(Capper_GBM_betas[cpg, which(colnames(Capper_GBM_betas) %in% HC_samps_train)], na.rm=TRUE)
  summ_betas_train[cpg, "GBM"] <- mean(Capper_GBM_betas[cpg, which(colnames(Capper_GBM_betas) %in% GBM_samps_train)], na.rm=TRUE)
  summ_betas_train[cpg, "OLG"] <- mean(Capper_GBM_betas[cpg, which(colnames(Capper_GBM_betas) %in% OLG_samps_train)], na.rm=TRUE)
  summ_betas_train[cpg, "AST-LG"] <- mean(Capper_GBM_betas[cpg, which(colnames(Capper_GBM_betas) %in% AST_LG_samps_train)], na.rm=TRUE)
  summ_betas_train[cpg, "AST-HG"] <- mean(Capper_GBM_betas[cpg, which(colnames(Capper_GBM_betas) %in% AST_HG_samps_train)], na.rm=TRUE)
}


library(pheatmap)
heat_annot <- data.frame(summ_betas[,c(6:9)])
colnames(heat_annot) <- c("GBM L0", "OLG L0", "AST-HG L0", "AST-LG L0")
ann_colors <- list(`OLG L0` = c(Yes = "#3D7636", No = "white"),
                   `GBM L0` = c(Yes = "#BA6A78", No = "white"),
                   `AST-HG L0` = c(Yes = "#342186", No = "white"),
                   `AST-LG L0` = c(Yes = "#9ECAEC", No = "white"))
  


pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Manuscript/GIMiCC/L0_heatmap.pdf",
    onefile=TRUE, height=7, width=5)
pheatmap(
  mat=summ_betas[, c(1:5)],
  main="Total Dataset",
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = T, #samples
  #labels_col = rownames(TC_Library),
  annotation_row = heat_annot,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  cluster_cols=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10,
  annotation_legend=FALSE
)

pheatmap(
  mat=summ_betas_train[, c(1:5)],
  main="Training Set",
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = T, #samples
  #labels_col = rownames(TC_Library),
  annotation_row = heat_annot,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  cluster_cols=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10,
  annotation_legend=FALSE
)

pheatmap(
  mat=summ_betas_test[, c(1:5)],
  main="Test Set",
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = T, #samples
  #labels_col = rownames(TC_Library),
  annotation_row = heat_annot,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  cluster_cols=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10,
  annotation_legend=FALSE
)
dev.off()

###################################### Figure 2B,2C,2D,2E
library(dplyr)
library(readr)
library(readxl)
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Capper/Processed_data/Capper_pheno_agg2.RDA")
Capper_Supp_Info <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Capper_Supp_Info.xlsx")
rownames(Capper_pheno) <- paste0(Capper_pheno$Slide,"_",Capper_pheno$Array)
Capper_pheno$ID <- substr(rownames(Capper_pheno), 12, 100)
pheno_df <- left_join(Capper_pheno, Capper_Supp_Info)
rownames(pheno_df) <- paste0(pheno_df$Slide,"_",pheno_df$Array)
pheno_df$ID <- rownames(pheno_df)

load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/Capper_deconvo_results.RDA")
GBM_proj <- 100*GBM_proj/rowSums(GBM_proj)
OLG_proj <- 100*OLG_proj/rowSums(OLG_proj)
AST_HG_proj <- 100*AST_HG_proj/rowSums(AST_HG_proj)
AST_LG_proj <- 100*AST_LG_proj/rowSums(AST_LG_proj)

GBM_proj$ID <- rownames(GBM_proj)
OLG_proj$ID <- rownames(OLG_proj)
AST_HG_proj$ID <- rownames(AST_HG_proj)
AST_LG_proj$ID <- rownames(AST_LG_proj)

GBM_proj <- left_join(GBM_proj, pheno_df, by="ID")
OLG_proj <- left_join(OLG_proj, pheno_df, by="ID")
AST_HG_proj <- left_join(AST_HG_proj, pheno_df, by="ID")
AST_LG_proj <- left_join(AST_LG_proj, pheno_df, by="ID")

GBM_proj$Capper_purity <- 100*GBM_proj$`Estimated tumour purity according to TCGA / Ceccarelli et al. 2016`
OLG_proj$Capper_purity <- 100*OLG_proj$`Estimated tumour purity according to TCGA / Ceccarelli et al. 2016`
AST_HG_proj$Capper_purity <- 100*AST_HG_proj$`Estimated tumour purity according to TCGA / Ceccarelli et al. 2016`
AST_LG_proj$Capper_purity <- 100*AST_LG_proj$`Estimated tumour purity according to TCGA / Ceccarelli et al. 2016`

GBM_proj2 <- GBM_proj[GBM_proj$Broad_Class == "Glioblastoma",]
OLG_proj2 <- OLG_proj[OLG_proj$methylation.class.ch1 == "O IDH",]
AST_HG_proj2 <- AST_HG_proj[AST_HG_proj$methylation.class.ch1 == "A IDH, HG",]
AST_LG_proj2 <- AST_LG_proj[AST_LG_proj$methylation.class.ch1 == "A IDH",]

library(ggplot2)
library(ggpubr)
plot1 <- ggplot(GBM_proj2, aes(x=Capper_purity, y=Tumor)) + geom_point(aes(shape=Train_split),color="black",fill="#BA6A78")
plot1 <- plot1 + xlab("Capper Estimated Purity (%)")
plot1 <- plot1 + ylab("GIMiCC Estimated Purity (%)")
plot1 <- plot1 + labs(subtitle=paste0("GBM (n = ",sum(!is.na(GBM_proj2$Capper_purity)),")"))
plot1 <- plot1 + theme_classic()
plot1 <- plot1 + scale_shape_manual(values = c(21, 4))
plot1 <- plot1 + geom_abline(slope=1,intercept=0,lty=3)
plot1 <- plot1 + theme(legend.position = "none")

plot2 <- ggplot(OLG_proj2, aes(x=Capper_purity, y=Tumor)) + geom_point(aes(shape=Train_split),color="black", fill="#3D7636")
plot2 <- plot2 + xlab("Capper Estimated Purity (%)")
plot2 <- plot2 + ylab("GIMiCC Estimated Purity (%)")
plot2 <- plot2 + labs(subtitle=paste0("OLG (n = ",sum(!is.na(OLG_proj2$Capper_purity)),")"))
plot2 <- plot2 + theme_classic()
plot2 <- plot2 + scale_shape_manual(values = c(21, 4))
plot2 <- plot2 + geom_abline(slope=1,intercept=0,lty=3)
plot2 <- plot2 + theme(legend.position = "none")

plot3 <- ggplot(AST_HG_proj2, aes(x=Capper_purity, y=Tumor)) + geom_point(aes(shape=Train_split),color="black",fill="#342186")
plot3 <- plot3 + xlab("Capper Estimated Purity (%)")
plot3 <- plot3 + ylab("GIMiCC Estimated Purity (%)")
plot3 <- plot3 + labs(subtitle=paste0("AST-HG (n = ",sum(!is.na(AST_HG_proj2$Capper_purity)),")"))
plot3 <- plot3 + theme_classic()
plot3 <- plot3 + scale_shape_manual(values = c(21, 4))
plot3 <- plot3 + geom_abline(slope=1,intercept=0,lty=3)
plot3 <- plot3 + theme(legend.position = "none")

plot4 <- ggplot(AST_LG_proj2, aes(x=Capper_purity, y=Tumor)) + geom_point(aes(shape=Train_split),color="black",fill="#9ECAEC")
plot4 <- plot4 + xlab("Capper Estimated Purity (%)")
plot4 <- plot4 + ylab("GIMiCC Estimated Purity (%)")
plot4 <- plot4 + labs(subtitle=paste0("AST-LG (n = ",sum(!is.na(AST_LG_proj2$Capper_purity)),")"))
plot4 <- plot4 + theme_classic()
plot4 <- plot4 + scale_shape_manual(values = c(21, 4))
plot4 <- plot4 + geom_abline(slope=1,intercept=0,lty=3)
plot4 <- plot4 + theme(legend.position = "none")

pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Manuscript/GIMiCC/Capper_predicted_L0.pdf",
    onefile=TRUE, height=3, width=3)
plot1
plot2
plot3
plot4
dev.off()

################################## Figure 3 
library(readr)
library(ggplot2)
library(ggpubr)

setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets")
for (tt in c("GBM_L0","OLG_L0","AST_LG_L0","AST_HG_L0","CellLayers")){
  df <- read.csv(paste0(tt,"_enrich_chr.csv"))
  df$Location <- factor(df$Location, 
                        levels = c("chr22" , "chr21" , "chr20" , "chr19",
                                   "chr18", "chr17", "chr16", "chr15", "chr14",
                                   "chr13", "chr12", "chr11", "chr10", "chr9",
                                   "chr8", "chr7", "chr6", "chr5", "chr4",
                                   "chr3", "chr2", "chr1"))
  plot2 <- ggplot(df, aes(y = Location, x = OR)) 
  plot2 <- plot2 +  geom_point(shape = 18, size = 5)
  plot2 <- plot2 +  geom_errorbarh(aes(xmin = OR_CI_lo, 
                                       xmax = OR_CI_hi), 
                                   height = 0.25)
  plot2 <- plot2 +  geom_vline(xintercept = 1, color = "red",
                               linetype = "dashed", 
                               cex = 1, alpha = 0.5)
  plot2 <- plot2 +  xlab("Odds Ratio (95% CI)")
  plot2 <- plot2 +  ylab("Chromosome")
  plot2 <- plot2 +  theme_bw()
  
  df <- read.csv(paste0(tt,"_enrich_reltoisland.csv"))
  df$Location <- ifelse(df$Location == "N_Shelf", "N Shelf", 
                        ifelse(df$Location == "S_Shelf", "S Shelf", 
                               ifelse(df$Location == "N_Shore", "N Shore", 
                                      ifelse(df$Location == "S_Shore", "S Shore", 
                                             ifelse(df$Location == "OpenSea", "Open Sea", df$Location)))))
  df$Location <- factor(df$Location, 
                        levels = c("S Shelf" , "S Shore" , "Island" ,
                                  "N Shore", "N Shelf", "Open Sea"))
  
  plot2 <- ggplot(df, aes(y = Location, x = OR)) 
  plot2 <- plot2 +  geom_point(shape = 18, size = 5)
  plot2 <- plot2 +  geom_errorbarh(aes(xmin = OR_CI_lo, 
                                       xmax = OR_CI_hi), 
                                   height = 0.25)
  plot2 <- plot2 +  geom_vline(xintercept = 1, color = "red",
                               linetype = "dashed", 
                               cex = 1, alpha = 0.5)
  plot2 <- plot2 +  xlab("Odds Ratio (95% CI)")
  plot2 <- plot2 +  ylab("Proximity to CpG Islands")
  plot2 <- plot2 +  theme_bw()
  
  df <- read.csv(paste0(tt,"_enrich_location.csv"))
  df$Location <- factor(df$Location, 
                        levels = c("3'UTR" , "Body" , "1stExon" ,
                                   "5'UTR", "TSS200", "TSS1500"))
  plot3 <- ggplot(df, aes(y = Location, x = OR)) 
  plot3 <- plot3 +  geom_point(shape = 18, size = 5)
  plot3 <- plot3 +  geom_errorbarh(aes(xmin = OR_CI_lo, 
                                       xmax = OR_CI_hi), 
                                   height = 0.25)
  plot3 <- plot3 +  geom_vline(xintercept = 1, color = "red",
                               linetype = "dashed", 
                               cex = 1, alpha = 0.5)
  plot3 <- plot3 +  xlab("Odds Ratio (95% CI)")
  plot3 <- plot3 +  ylab("Genomic Region")
  plot3 <- plot3 +  theme_bw()
  
  pdf(paste0(tt,"_enrich_plots.pdf"), onefile=TRUE, width=4, height=3)
  print(plot2)
  print(plot2)
  print(plot3)
  dev.off()
}

################################ Figure 4
library(dplyr)
library(readr)
library(readxl)
library(ggplot2)
library(reshape2)
cell.types <- c("Endothelial","Stromal","Astrocyte","Microglia","Mono","Oligodendrocyte","GABA","GLU","Tumor",
                "Neu","Bmem","Bnv","CD4mem","CD4nv","CD8mem","CD8nv","Treg","NK")

load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Capper/Processed_data/Capper_pheno_agg2.RDA")
Capper_Supp_Info <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Capper_Supp_Info.xlsx")
rownames(Capper_pheno) <- paste0(Capper_pheno$Slide,"_",Capper_pheno$Array)
Capper_pheno$ID <- substr(rownames(Capper_pheno), 12, 100)
pheno_df <- left_join(Capper_pheno, Capper_Supp_Info)
rownames(pheno_df) <- paste0(pheno_df$Slide,"_",pheno_df$Array)
pheno_df$ID <- rownames(pheno_df)

load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/Capper_deconvo_results.RDA")
GBM_proj <- 100*GBM_proj/rowSums(GBM_proj)
OLG_proj <- 100*OLG_proj/rowSums(OLG_proj)
AST_HG_proj <- 100*AST_HG_proj/rowSums(AST_HG_proj)
AST_LG_proj <- 100*AST_LG_proj/rowSums(AST_LG_proj)

GBM_proj$ID <- rownames(GBM_proj)
OLG_proj$ID <- rownames(OLG_proj)
AST_HG_proj$ID <- rownames(AST_HG_proj)
AST_LG_proj$ID <- rownames(AST_LG_proj)

GBM_proj <- left_join(GBM_proj, pheno_df, by="ID")
OLG_proj <- left_join(OLG_proj, pheno_df, by="ID")
AST_HG_proj <- left_join(AST_HG_proj, pheno_df, by="ID")
AST_LG_proj <- left_join(AST_LG_proj, pheno_df, by="ID")

GBM_proj$Library <- "GBM"
OLG_proj$Library <- "OLG"
AST_HG_proj$Library <- "AST-HG"
AST_LG_proj$Library <- "AST-LG"

total_proj <- rbind(GBM_proj, OLG_proj, AST_HG_proj, AST_LG_proj)

####### Pannel 4A
total_proj <- rbind(GBM_proj, OLG_proj, AST_HG_proj, AST_LG_proj)
total_proj <- total_proj[total_proj$Broad_Class %in% c("Control_Healthy", "Glioma_IDH", "Glioblastoma"),]
total_proj <- total_proj[, c(cell.types, "Broad_Class","methylation.class.ch1","Train_split", "Library")]
total_proj <- melt(total_proj, id.vars = c("Broad_Class","methylation.class.ch1","Train_split", "Library"))
total_proj$TumorType <- ifelse(total_proj$Broad_Class == "Control_Healthy", "HC", 
                               ifelse(total_proj$Broad_Class == "Glioblastoma", "GBM", 
                                      ifelse(total_proj$methylation.class.ch1 == "O IDH", "OLG", 
                                             ifelse(total_proj$methylation.class.ch1 == "A IDH", "AST-LG", 
                                                    ifelse(total_proj$methylation.class.ch1 == "A IDH, HG", "AST-HG", NA
                                                    )))))
colnames(total_proj)[which(colnames(total_proj) == "variable")] <- "CellType"
colnames(total_proj)[which(colnames(total_proj) == "value")] <- "Prop"


total_proj$TumorType <- factor(total_proj$TumorType,
                                  levels = c("HC", "OLG", "AST-LG", "AST-HG", "GBM"))
total_proj$Library <- factor(total_proj$Library,
                               levels = c("OLG", "AST-LG", "AST-HG", "GBM"))

total_proj_test <- total_proj[total_proj$Train_split==FALSE,]
total_proj_train <- total_proj[total_proj$Train_split==TRUE,]

total_proj_test <- total_proj_test[total_proj_test$CellType == "Tumor",]
total_proj_train <- total_proj_train[total_proj_train$CellType == "Tumor",]

tmp_colors <- c(`OLG`="#3D7636",
                `AST-LG`="#9ECAEC",
                `AST-HG`="#342186",
                `GBM`="#BA6A78")
plot1 <- ggplot(total_proj_test, aes(x=TumorType, y=Prop, fill=Library)) + geom_boxplot(outlier.shape=NA)
plot1 <- plot1 + theme_classic()
plot1 <- plot1 + xlab("Tumor Sample Type")
plot1 <- plot1 + ylab("Tumor Purity")
plot1 <- plot1 + labs(fill="L0 Library Used")
plot1 <- plot1 + geom_point(position=position_jitterdodge(), shape=1, alpha = 0.3)
plot1 <- plot1 + scale_fill_manual(values = tmp_colors)
plot1 <- plot1 + theme(legend.position="top")

plot2 <- ggplot(total_proj_train, aes(x=TumorType, y=Prop, fill=Library)) + geom_boxplot(outlier.shape=NA)
plot2 <- plot2 + theme_classic()
plot2 <- plot2 + xlab("Tumor Sample Type")
plot2 <- plot2 + ylab("Tumor Purity")
plot2 <- plot2 + labs(fill="L0 Library Used")
plot2 <- plot2 + geom_point(position=position_jitterdodge(), shape=1, alpha = 0.3)
plot2 <- plot2 + scale_fill_manual(values = tmp_colors)
plot2 <- plot2 + theme(legend.position="top")

pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Manuscript/GIMiCC/Using_nonmatching_L0.pdf",
    onefile=TRUE, height=3, width=7)
plot1
plot2
dev.off()

############ Panel 4B.4C.4D
total_proj <- GBM_proj
total_proj <- total_proj[total_proj$Broad_Class %in% c("Control_Healthy", "Control_High_Infiltrate", "Control_Low_Yeild"),]
total_proj <- total_proj[, c(cell.types, "Broad_Class","methylation.class.ch1")]

cell.types2 <- c("Neuronal","Myeloid","Lymphoid","Angiogenic","Microglia","Astrocyte","Oligodendrocyte")
total_proj$Myeloid <- total_proj$Neu + total_proj$Mono
total_proj$Lymphoid <- total_proj$CD4nv + total_proj$CD4mem + total_proj$Treg + total_proj$CD8mem + total_proj$CD8nv+ total_proj$Bnv + total_proj$Bmem+ total_proj$NK
total_proj$Angiogenic <- total_proj$Endothelial + total_proj$Stromal
total_proj$Neuronal <- total_proj$GABA + total_proj$GLU
tmp_total_proj <- total_proj[, cell.types2]
tmp_total_proj <- 100*tmp_total_proj/rowSums(tmp_total_proj)
tmp_total_proj$Broad_Class <- total_proj$Broad_Class
tmp_total_proj$Region <- total_proj$methylation.class.ch1


library(RColorBrewer)
cell_colors <- c(Neuronal="#cc79a7",
                          Myeloid="#0072b2",
                          Lymphoid="#009e73",
                          Angiogenic="#d55e00",
                          Microglia="#56b4e9",
                          Astrocyte="#f0e442",
                          Oligodendrocyte="#e69f00")

tmp_pheno <- tmp_total_proj[tmp_total_proj$Broad_Class == "Control_Healthy",]
tmp_pheno <- tmp_pheno[order(tmp_pheno$Region, tmp_pheno$Oligodendrocyte, tmp_pheno$Neuronal),]
table(tmp_pheno$Region)
tmp_pheno$order1 <- c(1:nrow(tmp_pheno))
tmp_pheno <- melt(tmp_pheno, id.vars = c("Broad_Class","order1","Region"))
colnames(tmp_pheno)[which(colnames(tmp_pheno) == "variable")] <- "CellType"
colnames(tmp_pheno)[which(colnames(tmp_pheno) == "value")] <- "Prop"
plot1 <- ggplot(tmp_pheno, aes(fill=CellType, y=Prop, x=order1)) + geom_bar(position="stack", stat="identity")
plot1 <- plot1 + theme_classic()
plot1 <- plot1 + theme(axis.text.x=element_blank(),
                       legend.position = "none",
                       axis.ticks.x=element_blank())
plot1 <- plot1 + xlab("Healthy Control Samples")
plot1 <- plot1 + ylab("Proportion of Non-Tumor Fraction")
plot1 <- plot1 + scale_fill_manual(values=cell_colors)

tmp_pheno <- tmp_total_proj[tmp_total_proj$Broad_Class == "Control_High_Infiltrate",]
tmp_pheno <- tmp_pheno[order(tmp_pheno$Oligodendrocyte, tmp_pheno$Neuronal),]
tmp_pheno$order1 <- c(1:nrow(tmp_pheno))
tmp_pheno <- melt(tmp_pheno, id.vars = c("Broad_Class","order1","Region"))
colnames(tmp_pheno)[which(colnames(tmp_pheno) == "variable")] <- "CellType"
colnames(tmp_pheno)[which(colnames(tmp_pheno) == "value")] <- "Prop"
plot2 <- ggplot(tmp_pheno, aes(fill=CellType, y=Prop, x=order1)) + geom_bar(position="stack", stat="identity")
plot2 <- plot2 + theme_classic()
plot2 <- plot2 + theme(axis.text.x=element_blank(),
                       legend.position = "none",
                       axis.ticks.x=element_blank())
plot2 <- plot2 + xlab("High Infiltrate Samples")
plot2 <- plot2 + ylab("Proportion of Non-Tumor Fraction")
plot2 <- plot2 + scale_fill_manual(values=cell_colors)


tmp_pheno <- tmp_total_proj[tmp_total_proj$Broad_Class == "Control_Low_Yeild",]
tmp_pheno <- tmp_pheno[order(tmp_pheno$Oligodendrocyte, tmp_pheno$Neuronal),]
tmp_pheno$order1 <- c(1:nrow(tmp_pheno))
tmp_pheno <- melt(tmp_pheno, id.vars = c("Broad_Class","order1","Region"))
colnames(tmp_pheno)[which(colnames(tmp_pheno) == "variable")] <- "CellType"
colnames(tmp_pheno)[which(colnames(tmp_pheno) == "value")] <- "Prop"
plot3 <- ggplot(tmp_pheno, aes(fill=CellType, y=Prop, x=order1)) + geom_bar(position="stack", stat="identity")
plot3 <- plot3 + theme_classic()
plot3 <- plot3 + theme(axis.text.x=element_blank(),
                       legend.position = "top",
                       axis.ticks.x=element_blank())
plot3 <- plot3 + xlab("Low Yeild Samples")
plot3 <- plot3 + ylab("Proportion of Non-Tumor Fraction")
plot3 <- plot3 + scale_fill_manual(values=cell_colors)


pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Manuscript/GIMiCC/GIMiCC_October2023/Capper_Healthy_Barplots.pdf", width = 4, height = 4)
plot1
plot2
plot3
dev.off()


#################### Capper plotting for supplement
library(dplyr)
library(readr)
library(readxl)
library(ggplot2)
library(reshape2)
cell.types <- c("Endothelial","Stromal","Astrocyte","Microglia","Mono","Oligodendrocyte","GABA","GLU","Tumor",
                "Neu","Bmem","Bnv","CD4mem","CD4nv","CD8mem","CD8nv","Treg","NK")

load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Capper/Processed_data/Capper_pheno_agg2.RDA")
Capper_Supp_Info <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Capper_Supp_Info.xlsx")
rownames(Capper_pheno) <- paste0(Capper_pheno$Slide,"_",Capper_pheno$Array)
Capper_pheno$ID <- substr(rownames(Capper_pheno), 12, 100)
pheno_df <- left_join(Capper_pheno, Capper_Supp_Info)
rownames(pheno_df) <- paste0(pheno_df$Slide,"_",pheno_df$Array)
pheno_df$ID <- rownames(pheno_df)

pheno_df$TumorType <- ifelse(pheno_df$Broad_Class == "Control_Healthy", "HC", 
                               ifelse(pheno_df$Broad_Class == "Glioblastoma", "GBM", 
                                      ifelse(pheno_df$methylation.class.ch1 == "O IDH", "OLG", 
                                             ifelse(pheno_df$methylation.class.ch1 == "A IDH", "AST-LG", 
                                                    ifelse(pheno_df$methylation.class.ch1 == "A IDH, HG", "AST-HG", NA
                                                    )))))


load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/Capper_deconvo_results.RDA")
GBM_proj <- 100*GBM_proj/rowSums(GBM_proj)
OLG_proj <- 100*OLG_proj/rowSums(OLG_proj)
AST_HG_proj <- 100*AST_HG_proj/rowSums(AST_HG_proj)
AST_LG_proj <- 100*AST_LG_proj/rowSums(AST_LG_proj)

GBM_proj$ID <- rownames(GBM_proj)
OLG_proj$ID <- rownames(OLG_proj)
AST_HG_proj$ID <- rownames(AST_HG_proj)
AST_LG_proj$ID <- rownames(AST_LG_proj)

GBM_proj <- left_join(GBM_proj, pheno_df, by="ID")
OLG_proj <- left_join(OLG_proj, pheno_df, by="ID")
AST_HG_proj <- left_join(AST_HG_proj, pheno_df, by="ID")
AST_LG_proj <- left_join(AST_LG_proj, pheno_df, by="ID")

GBM_proj$Library <- "GBM"
OLG_proj$Library <- "OLG"
AST_HG_proj$Library <- "AST-HG"
AST_LG_proj$Library <- "AST-LG"

Capper_proj <- rbind(GBM_proj, OLG_proj, AST_HG_proj, AST_LG_proj)
Capper_proj <- Capper_proj[Capper_proj$Library == Capper_proj$TumorType,]
Capper_proj <- Capper_proj[!is.na(Capper_proj$TumorType),]


cell.types2 <- c("Tumor","Neuronal","Myeloid","Lymphoid","Angiogenic","Microglia","Astrocyte","Oligodendrocyte")
Capper_proj$Myeloid <- Capper_proj$Neu + Capper_proj$Mono
Capper_proj$Lymphoid <- Capper_proj$CD4nv + Capper_proj$CD4mem + Capper_proj$Treg + Capper_proj$CD8mem + Capper_proj$CD8nv+ Capper_proj$Bnv + Capper_proj$Bmem+ Capper_proj$NK
Capper_proj$Angiogenic <- Capper_proj$Endothelial + Capper_proj$Stromal
Capper_proj$Neuronal <- Capper_proj$GABA + Capper_proj$GLU
tmp_Capper_proj <- Capper_proj[, cell.types2]
tmp_Capper_proj <- 100*tmp_Capper_proj/rowSums(tmp_Capper_proj)
tmp_Capper_proj$ID <- Capper_proj$ID
tmp_Capper_proj$TumorType <- Capper_proj$TumorType

cell_colors <- c(Tumor="black",
                 Neuronal="#cc79a7",
                 Myeloid="#0072b2",
                 Lymphoid="#009e73",
                 Angiogenic="#d55e00",
                 Microglia="#56b4e9",
                 Astrocyte="#f0e442",
                 Oligodendrocyte="#e69f00")

library(dplyr)
library(ggplot2)
library(reshape2)
tmp_pheno <- tmp_Capper_proj[which(tmp_Capper_proj$TumorType == "GBM"),]
tmp_pheno <- tmp_pheno[order(tmp_pheno$Tumor,tmp_pheno$Oligodendrocyte,tmp_pheno$Myeloid,tmp_pheno$Lymphoid),]
tmp_pheno$order1 <- c(1:nrow(tmp_pheno))
tmp_pheno <- melt(tmp_pheno, id.vars = c("ID","order1","TumorType"))
colnames(tmp_pheno) <- c("ID","order1","TumorType","Celltype","Proportion")
plot1 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=order1)) + geom_bar(position="stack", stat="identity")
plot1 <- plot1 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot1 <- plot1 + scale_fill_manual(values=cell_colors)
plot1 <- plot1 + xlab("GBM") + ylab("Composition (%)")

tmp_pheno <- tmp_Capper_proj[which(tmp_Capper_proj$TumorType == "OLG"),]
tmp_pheno <- tmp_pheno[order(tmp_pheno$Tumor,tmp_pheno$Oligodendrocyte,tmp_pheno$Myeloid,tmp_pheno$Lymphoid),]
tmp_pheno$order1 <- c(1:nrow(tmp_pheno))
tmp_pheno <- melt(tmp_pheno, id.vars = c("ID","order1","TumorType"))
colnames(tmp_pheno) <- c("ID","order1","TumorType","Celltype","Proportion")
plot2 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=order1)) + geom_bar(position="stack", stat="identity")
plot2 <- plot2 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot2 <- plot2 + scale_fill_manual(values=cell_colors)
plot2 <- plot2 + xlab("OLG") + ylab("Composition (%)")

tmp_pheno <- tmp_Capper_proj[which(tmp_Capper_proj$TumorType == "AST-HG"),]
tmp_pheno <- tmp_pheno[order(tmp_pheno$Tumor,tmp_pheno$Oligodendrocyte,tmp_pheno$Myeloid,tmp_pheno$Lymphoid),]
tmp_pheno$order1 <- c(1:nrow(tmp_pheno))
tmp_pheno <- melt(tmp_pheno, id.vars = c("ID","order1","TumorType"))
colnames(tmp_pheno) <- c("ID","order1","TumorType","Celltype","Proportion")
plot3 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=order1)) + geom_bar(position="stack", stat="identity")
plot3 <- plot3 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot3 <- plot3 + scale_fill_manual(values=cell_colors)
plot3 <- plot3 + xlab("AST-HG") + ylab("Composition (%)")

tmp_pheno <- tmp_Capper_proj[which(tmp_Capper_proj$TumorType == "AST-LG"),]
tmp_pheno <- tmp_pheno[order(tmp_pheno$Tumor,tmp_pheno$Oligodendrocyte,tmp_pheno$Myeloid,tmp_pheno$Lymphoid),]
tmp_pheno$order1 <- c(1:nrow(tmp_pheno))
tmp_pheno <- melt(tmp_pheno, id.vars = c("ID","order1","TumorType"))
colnames(tmp_pheno) <- c("ID","order1","TumorType","Celltype","Proportion")
plot4 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=order1)) + geom_bar(position="stack", stat="identity")
plot4 <- plot4 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot4 <- plot4 + scale_fill_manual(values=cell_colors)
plot4 <- plot4 + xlab("AST-LG") + ylab("Composition (%)")

pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Manuscript/GIMiCC/GIMiCC_October2023/Capper_Barplots.pdf", width = 3, height = 3)
plot1
plot2
plot3
plot4
dev.off()


cell.types2 <- c("Microglia","NK","Mono","Neu","CD8nv","CD8mem",
                 "CD4nv","CD4mem","Treg","Bnv","Bmem")
cell_colors <- c(Microglia="#56b4e9",
                 NK="#009191",
                 Mono="#B56cfe",
                 Neu="#480091",
                 CD4nv="#DA6C00",
                 CD4mem="#924900",
                 CD8mem="#23FD23",
                 CD8nv="#FFFF6c",
                 Treg="grey4",
                 Bnv="#FFB5DA",
                 Bmem="#FE6CB6")

tmp_Capper_proj <- Capper_proj[, cell.types2]
#tmp_Capper_proj <- 100*tmp_Capper_proj/rowSums(tmp_Capper_proj)
tmp_Capper_proj$ID <- Capper_proj$ID
tmp_Capper_proj$TumorType <- Capper_proj$TumorType

tmp_pheno <- tmp_Capper_proj[which(tmp_Capper_proj$TumorType == "GBM"),]
tmp_pheno <- tmp_pheno[, c("ID",cell.types2)]
tmp_pheno <- melt(tmp_pheno, id.vars = c("ID"))
colnames(tmp_pheno) <- c("ID","Celltype","Proportion")
plot1 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=reorder(ID,Proportion,sum,decreasing=TRUE))) + geom_bar(position="stack", stat="identity")
plot1 <- plot1 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot1 <- plot1 + scale_fill_manual(values=cell_colors)
plot1 <- plot1 + xlab("GBM") + ylab("Composition (%)")
plot1 <- plot1 + ylim(0,60)

tmp_pheno <- tmp_Capper_proj[which(tmp_Capper_proj$TumorType == "OLG"),]
tmp_pheno <- tmp_pheno[, c("ID",cell.types2)]
tmp_pheno <- melt(tmp_pheno, id.vars = c("ID"))
colnames(tmp_pheno) <- c("ID","Celltype","Proportion")
plot2 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=reorder(ID,Proportion,sum,decreasing=TRUE))) + geom_bar(position="stack", stat="identity")
plot2 <- plot2 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot2 <- plot2 + scale_fill_manual(values=cell_colors)
plot2 <- plot2 + xlab("OLG") + ylab("Composition (%)")
plot2 <- plot2 + ylim(0,60)

tmp_pheno <- tmp_Capper_proj[which(tmp_Capper_proj$TumorType == "AST-HG"),]
tmp_pheno <- tmp_pheno[, c("ID",cell.types2)]
tmp_pheno <- melt(tmp_pheno, id.vars = c("ID"))
colnames(tmp_pheno) <- c("ID","Celltype","Proportion")
plot3 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=reorder(ID,Proportion,sum,decreasing=TRUE))) + geom_bar(position="stack", stat="identity")
plot3 <- plot3 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot3 <- plot3 + scale_fill_manual(values=cell_colors)
plot3 <- plot3 + xlab("AST-HG") + ylab("Composition (%)")
plot3 <- plot3 + ylim(0,60)

tmp_pheno <- tmp_Capper_proj[which(tmp_Capper_proj$TumorType == "AST-LG"),]
tmp_pheno <- tmp_pheno[, c("ID",cell.types2)]
tmp_pheno <- melt(tmp_pheno, id.vars = c("ID"))
colnames(tmp_pheno) <- c("ID","Celltype","Proportion")
plot4 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=reorder(ID,Proportion,sum,decreasing=TRUE))) + geom_bar(position="stack", stat="identity")
plot4 <- plot4 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot4 <- plot4 + scale_fill_manual(values=cell_colors)
plot4 <- plot4 + xlab("AST-LG") + ylab("Composition (%)")
plot4 <- plot4 + ylim(0,60)

pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Manuscript/GIMiCC/GIMiCC_October2023/Capper_Barplots_Immunemicro.pdf", width = 3, height = 3)
plot1
plot2
plot3
plot4
dev.off()

###################

################################ Figure 5
library(dplyr)
library(readr)
library(readxl)
library(ggplot2)
library(reshape2)
cell.types <- c("Endothelial","Stromal","Astrocyte","Microglia","Mono","Oligodendrocyte","GABA","GLU","Tumor",
                "Neu","Bmem","Bnv","CD4mem","CD4nv","CD8mem","CD8nv","Treg","NK")


load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Brain_Tumor_Deconv_Ref/Processed_data/TCGA_ref/TCGAall_pheno.RDA")
TCGA_Supp_Info <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/TCGA_Supp_Info.xlsx")
TCGA_pheno <- left_join(TCGA_pheno, TCGA_Supp_Info, by="submitter_id")

df <- read_tsv(file = "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/TCGA_extra_info.tsv")
df <- df[c(5,6,7),]
df <- t(df)
colnames(df) <- df[1,]
df <- data.frame(df[-c(1,2),])
df$submitter_id <- rownames(df)
TCGA_pheno <- left_join(TCGA_pheno, df, by="submitter_id")

TCGA_pheno$Tumor_Type2 <- ifelse(TCGA_pheno$`IDH/codel subtype.y` == "IDHwt", "GBM",
                                 ifelse(TCGA_pheno$`IDH/codel subtype.y`== "IDHmut-codel", "OLG", 
                                               ifelse(TCGA_pheno$`IDH/codel subtype.y`== "IDHmut-non-codel", "AST", "N/A")))

TCGA_pheno$`TERT promoter status.y` <- ifelse(is.na(TCGA_pheno$`TERT promoter status.y`), "N/A", TCGA_pheno$`TERT promoter status.y`)
TCGA_pheno$`Chr 7 gain/Chr 10 loss` <- ifelse(is.na(TCGA_pheno$`Chr 7 gain/Chr 10 loss`), "N/A", TCGA_pheno$`TERT promoter status.y`)
TCGA_pheno$EGFR <- ifelse(is.na(TCGA_pheno$EGFR), "N/A", TCGA_pheno$EGFR)
TCGA_pheno$CDKN2A <- ifelse(is.na(TCGA_pheno$CDKN2A), "N/A", TCGA_pheno$CDKN2A)
TCGA_pheno$CDKN2B <- ifelse(is.na(TCGA_pheno$CDKN2B), "N/A", TCGA_pheno$CDKN2B)

TCGA_pheno$Tumor_Type2 <- ifelse(TCGA_pheno$Tumor_Type2 != "AST", TCGA_pheno$Tumor_Type2, 
                                 ifelse(TCGA_pheno$`TERT promoter status.y`== "Mutant" | 
                                        TCGA_pheno$`Chr 7 gain/Chr 10 loss`== "Gain chr 7 & loss chr 10" |
                                        TCGA_pheno$EGFR == "amp_rec" |
                                        TCGA_pheno$CDKN2A== "homdel_rec" |
                                        TCGA_pheno$CDKN2B== "homdel_rec" , "AST-HG", "AST-LG" ))


load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/TCGA_deconvo_results.RDA")
GBM_proj <- 100*GBM_proj/rowSums(GBM_proj)
OLG_proj <- 100*OLG_proj/rowSums(OLG_proj)
AST_HG_proj <- 100*AST_HG_proj/rowSums(AST_HG_proj)
AST_LG_proj <- 100*AST_LG_proj/rowSums(AST_LG_proj)

GBM_proj$methylID <- rownames(GBM_proj)
OLG_proj$methylID <- rownames(OLG_proj)
AST_HG_proj$methylID <- rownames(AST_HG_proj)
AST_LG_proj$methylID <- rownames(AST_LG_proj)

GBM_proj <- left_join(GBM_proj, TCGA_pheno, by="methylID")
OLG_proj <- left_join(OLG_proj, TCGA_pheno, by="methylID")
AST_HG_proj <- left_join(AST_HG_proj, TCGA_pheno, by="methylID")
AST_LG_proj <- left_join(AST_LG_proj, TCGA_pheno, by="methylID")

GBM_proj$Library <- "GBM"
OLG_proj$Library <- "OLG"
AST_HG_proj$Library <- "AST-HG"
AST_LG_proj$Library <- "AST-LG"

TCGA_proj <- rbind(GBM_proj, OLG_proj, AST_HG_proj, AST_LG_proj)
TCGA_proj <- TCGA_proj[TCGA_proj$Library == TCGA_proj$Tumor_Type2,]

cell.types2 <- c("Tumor","Neuronal","Myeloid","Lymphoid","Angiogenic","Microglia","Astrocyte","Oligodendrocyte")
TCGA_proj$Myeloid <- TCGA_proj$Neu + TCGA_proj$Mono
TCGA_proj$Lymphoid <- TCGA_proj$CD4nv + TCGA_proj$CD4mem + TCGA_proj$Treg + TCGA_proj$CD8mem + TCGA_proj$CD8nv+ TCGA_proj$Bnv + TCGA_proj$Bmem+ TCGA_proj$NK
TCGA_proj$Angiogenic <- TCGA_proj$Endothelial + TCGA_proj$Stromal
TCGA_proj$Neuronal <- TCGA_proj$GABA + TCGA_proj$GLU
tmp_TCGA_proj <- TCGA_proj[, cell.types2]
tmp_TCGA_proj <- 100*tmp_TCGA_proj/rowSums(tmp_TCGA_proj)
tmp_TCGA_proj$methylID <- TCGA_proj$methylID
tmp_TCGA_proj$TumorType <- TCGA_proj$Tumor_Type2


library(RColorBrewer)
cell_colors <- c(Tumor="black",
                 Neuronal="#cc79a7",
                 Myeloid="#0072b2",
                 Lymphoid="#009e73",
                 Angiogenic="#d55e00",
                 Microglia="#56b4e9",
                 Astrocyte="#f0e442",
                 Oligodendrocyte="#e69f00")

library(dplyr)
library(ggplot2)
library(reshape2)
tmp_pheno <- tmp_TCGA_proj[which(tmp_TCGA_proj$TumorType == "GBM"),]
tmp_pheno <- tmp_pheno[order(tmp_pheno$Tumor,tmp_pheno$Oligodendrocyte,tmp_pheno$Myeloid,tmp_pheno$Lymphoid),]
tmp_pheno$order1 <- c(1:nrow(tmp_pheno))
tmp_pheno <- melt(tmp_pheno, id.vars = c("methylID","order1","TumorType"))
colnames(tmp_pheno) <- c("methylID","order1","TumorType","Celltype","Proportion")
plot1 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=order1)) + geom_bar(position="stack", stat="identity")
plot1 <- plot1 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot1 <- plot1 + scale_fill_manual(values=cell_colors)
plot1 <- plot1 + xlab("GBM") + ylab("Composition (%)")

tmp_pheno <- tmp_TCGA_proj[which(tmp_TCGA_proj$TumorType == "OLG"),]
tmp_pheno <- tmp_pheno[order(tmp_pheno$Tumor,tmp_pheno$Oligodendrocyte,tmp_pheno$Myeloid,tmp_pheno$Lymphoid),]
tmp_pheno$order1 <- c(1:nrow(tmp_pheno))
tmp_pheno <- melt(tmp_pheno, id.vars = c("methylID","order1","TumorType"))
colnames(tmp_pheno) <- c("methylID","order1","TumorType","Celltype","Proportion")
plot2 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=order1)) + geom_bar(position="stack", stat="identity")
plot2 <- plot2 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot2 <- plot2 + scale_fill_manual(values=cell_colors)
plot2 <- plot2 + xlab("OLG") + ylab("Composition (%)")

tmp_pheno <- tmp_TCGA_proj[which(tmp_TCGA_proj$TumorType == "AST-HG"),]
tmp_pheno <- tmp_pheno[order(tmp_pheno$Tumor,tmp_pheno$Oligodendrocyte,tmp_pheno$Myeloid,tmp_pheno$Lymphoid),]
tmp_pheno$order1 <- c(1:nrow(tmp_pheno))
tmp_pheno <- melt(tmp_pheno, id.vars = c("methylID","order1","TumorType"))
colnames(tmp_pheno) <- c("methylID","order1","TumorType","Celltype","Proportion")
plot3 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=order1)) + geom_bar(position="stack", stat="identity")
plot3 <- plot3 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot3 <- plot3 + scale_fill_manual(values=cell_colors)
plot3 <- plot3 + xlab("AST-HG") + ylab("Composition (%)")

tmp_pheno <- tmp_TCGA_proj[which(tmp_TCGA_proj$TumorType == "AST-LG"),]
tmp_pheno <- tmp_pheno[order(tmp_pheno$Tumor,tmp_pheno$Oligodendrocyte,tmp_pheno$Myeloid,tmp_pheno$Lymphoid),]
tmp_pheno$order1 <- c(1:nrow(tmp_pheno))
tmp_pheno <- melt(tmp_pheno, id.vars = c("methylID","order1","TumorType"))
colnames(tmp_pheno) <- c("methylID","order1","TumorType","Celltype","Proportion")
plot4 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=order1)) + geom_bar(position="stack", stat="identity")
plot4 <- plot4 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot4 <- plot4 + scale_fill_manual(values=cell_colors)
plot4 <- plot4 + xlab("AST-LG") + ylab("Composition (%)")

pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Manuscript/GIMiCC/GIMiCC_October2023/TCGA_Barplots.pdf", width = 3, height = 3)
plot1
plot2
plot3
plot4
dev.off()


cell.types2 <- c("Microglia","NK","Mono","Neu","CD8nv","CD8mem",
                 "CD4nv","CD4mem","Treg","Bnv","Bmem")
cell_colors <- c(Microglia="#56b4e9",
                 NK="#009191",
                 Mono="#B56cfe",
                 Neu="#480091",
                 CD4nv="#DA6C00",
                 CD4mem="#924900",
                 CD8mem="#23FD23",
                 CD8nv="#FFFF6c",
                 Treg="grey4",
                 Bnv="#FFB5DA",
                 Bmem="#FE6CB6")

tmp_TCGA_proj <- TCGA_proj[, cell.types2]
#tmp_TCGA_proj <- 100*tmp_TCGA_proj/rowSums(tmp_TCGA_proj)
tmp_TCGA_proj$methylID <- TCGA_proj$methylID
tmp_TCGA_proj$TumorType <- TCGA_proj$Tumor_Type2

tmp_pheno <- TCGA_pheno[which(TCGA_pheno$Tumor_Type2 == "GBM"),]
tmp_pheno <- left_join(tmp_pheno, tmp_TCGA_proj, by = "methylID")
tmp_pheno <- tmp_pheno[, c("methylID",cell.types2)]
tmp_pheno <- melt(tmp_pheno, id.vars = c("methylID"))
colnames(tmp_pheno) <- c("methylID","Celltype","Proportion")
plot1 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=reorder(methylID,Proportion,sum,decreasing=TRUE))) + geom_bar(position="stack", stat="identity")
plot1 <- plot1 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot1 <- plot1 + scale_fill_manual(values=cell_colors)
plot1 <- plot1 + xlab("GBM") + ylab("Composition (%)")
plot1 <- plot1 + ylim(0,60)

tmp_pheno <- TCGA_pheno[which(TCGA_pheno$Tumor_Type2 == "OLG"),]
tmp_pheno <- left_join(tmp_pheno, tmp_TCGA_proj, by = "methylID")
tmp_pheno <- tmp_pheno[, c("methylID",cell.types2)]
tmp_pheno <- melt(tmp_pheno, id.vars = c("methylID"))
colnames(tmp_pheno) <- c("methylID","Celltype","Proportion")
plot2 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=reorder(methylID,Proportion,sum,decreasing=TRUE))) + geom_bar(position="stack", stat="identity")
plot2 <- plot2 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot2 <- plot2 + scale_fill_manual(values=cell_colors)
plot2 <- plot2 + xlab("OLG") + ylab("Composition (%)")
plot2 <- plot2 + ylim(0,60)

tmp_pheno <- TCGA_pheno[which(TCGA_pheno$Tumor_Type2 == "AST-HG"),]
tmp_pheno <- left_join(tmp_pheno, tmp_TCGA_proj, by = "methylID")
tmp_pheno <- tmp_pheno[, c("methylID",cell.types2)]
tmp_pheno <- melt(tmp_pheno, id.vars = c("methylID"))
colnames(tmp_pheno) <- c("methylID","Celltype","Proportion")
plot3 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=reorder(methylID,Proportion,sum,decreasing=TRUE))) + geom_bar(position="stack", stat="identity")
plot3 <- plot3 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot3 <- plot3 + scale_fill_manual(values=cell_colors)
plot3 <- plot3 + xlab("AST-HG") + ylab("Composition (%)")
plot3 <- plot3 + ylim(0,60)

tmp_pheno <- TCGA_pheno[which(TCGA_pheno$Tumor_Type2 == "AST-LG"),]
tmp_pheno <- left_join(tmp_pheno, tmp_TCGA_proj, by = "methylID")
tmp_pheno <- tmp_pheno[, c("methylID",cell.types2)]
tmp_pheno <- melt(tmp_pheno, id.vars = c("methylID"))
colnames(tmp_pheno) <- c("methylID","Celltype","Proportion")
plot4 <- ggplot(tmp_pheno, aes(fill=Celltype, y=Proportion, x=reorder(methylID,Proportion,sum,decreasing=TRUE))) + geom_bar(position="stack", stat="identity")
plot4 <- plot4 + theme_classic() + theme(axis.text.x=element_blank(), legend.position = "none",axis.ticks.x=element_blank())
plot4 <- plot4 + scale_fill_manual(values=cell_colors)
plot4 <- plot4 + xlab("AST-LG") + ylab("Composition (%)")
plot4 <- plot4 + ylim(0,60)

pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Manuscript/GIMiCC/GIMiCC_October2023/TCGA_Barplots_Immunemicro.pdf", width = 3, height = 3)
plot1
plot2
plot3
plot4
dev.off()


########## Figure 6
library(dplyr)
library(readr)
library(readxl)
library(ggplot2)
library(reshape2)
cell.types <- c("Endothelial","Stromal","Astrocyte","Microglia","Mono","Oligodendrocyte","GABA","GLU","Tumor",
                "Neu","Bmem","Bnv","CD4mem","CD4nv","CD8mem","CD8nv","Treg","NK")


load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Brain_Tumor_Deconv_Ref/Processed_data/TCGA_ref/TCGAall_pheno.RDA")
TCGA_Supp_Info <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/TCGA_Supp_Info.xlsx")
TCGA_pheno <- left_join(TCGA_pheno, TCGA_Supp_Info, by="submitter_id")

df <- read_tsv(file = "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/TCGA_extra_info.tsv")
df <- df[c(5,6,7),]
df <- t(df)
colnames(df) <- df[1,]
df <- data.frame(df[-c(1,2),])
df$submitter_id <- rownames(df)
TCGA_pheno <- left_join(TCGA_pheno, df, by="submitter_id")

TCGA_pheno$Tumor_Type2 <- ifelse(TCGA_pheno$`IDH/codel subtype.y` == "IDHwt", "GBM",
                                 ifelse(TCGA_pheno$`IDH/codel subtype.y`== "IDHmut-codel", "OLG", 
                                        ifelse(TCGA_pheno$`IDH/codel subtype.y`== "IDHmut-non-codel", "AST", "N/A")))

TCGA_pheno$`TERT promoter status.y` <- ifelse(is.na(TCGA_pheno$`TERT promoter status.y`), "N/A", TCGA_pheno$`TERT promoter status.y`)
TCGA_pheno$`Chr 7 gain/Chr 10 loss` <- ifelse(is.na(TCGA_pheno$`Chr 7 gain/Chr 10 loss`), "N/A", TCGA_pheno$`TERT promoter status.y`)
TCGA_pheno$EGFR <- ifelse(is.na(TCGA_pheno$EGFR), "N/A", TCGA_pheno$EGFR)
TCGA_pheno$CDKN2A <- ifelse(is.na(TCGA_pheno$CDKN2A), "N/A", TCGA_pheno$CDKN2A)
TCGA_pheno$CDKN2B <- ifelse(is.na(TCGA_pheno$CDKN2B), "N/A", TCGA_pheno$CDKN2B)

TCGA_pheno$Tumor_Type2 <- ifelse(TCGA_pheno$Tumor_Type2 != "AST", TCGA_pheno$Tumor_Type2, 
                                 ifelse(TCGA_pheno$`TERT promoter status.y`== "Mutant" | 
                                          TCGA_pheno$`Chr 7 gain/Chr 10 loss`== "Gain chr 7 & loss chr 10" |
                                          TCGA_pheno$EGFR == "amp_rec" |
                                          TCGA_pheno$CDKN2A== "homdel_rec" |
                                          TCGA_pheno$CDKN2B== "homdel_rec" , "AST-HG", "AST-LG" ))


load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/TCGA_deconvo_results.RDA")
GBM_proj <- 100*GBM_proj/rowSums(GBM_proj)
OLG_proj <- 100*OLG_proj/rowSums(OLG_proj)
AST_HG_proj <- 100*AST_HG_proj/rowSums(AST_HG_proj)
AST_LG_proj <- 100*AST_LG_proj/rowSums(AST_LG_proj)

GBM_proj$methylID <- rownames(GBM_proj)
OLG_proj$methylID <- rownames(OLG_proj)
AST_HG_proj$methylID <- rownames(AST_HG_proj)
AST_LG_proj$methylID <- rownames(AST_LG_proj)

GBM_proj <- left_join(GBM_proj, TCGA_pheno, by="methylID")
OLG_proj <- left_join(OLG_proj, TCGA_pheno, by="methylID")
AST_HG_proj <- left_join(AST_HG_proj, TCGA_pheno, by="methylID")
AST_LG_proj <- left_join(AST_LG_proj, TCGA_pheno, by="methylID")

GBM_proj$Library <- "GBM"
OLG_proj$Library <- "OLG"
AST_HG_proj$Library <- "AST-HG"
AST_LG_proj$Library <- "AST-LG"

TCGA_proj <- rbind(GBM_proj, OLG_proj, AST_HG_proj, AST_LG_proj)
TCGA_proj <- TCGA_proj[TCGA_proj$Library == TCGA_proj$Tumor_Type2,]
TCGA_proj$high_infil <- ifelse(TCGA_proj$Mono + TCGA_proj$CD8nv + TCGA_proj$CD4mem + TCGA_proj$CD4nv + TCGA_proj$CD4mem + TCGA_proj$Treg + TCGA_proj$Bnv + TCGA_proj$Bmem +  TCGA_proj$Neu > 20, "Yes", "No")
TCGA_proj$high_infil_sum <- TCGA_proj$Mono + TCGA_proj$CD8nv + TCGA_proj$CD4mem + TCGA_proj$CD4nv + TCGA_proj$CD4mem + TCGA_proj$Treg + TCGA_proj$Bnv + TCGA_proj$Bmem +  TCGA_proj$Neu 
TCGA_proj$CPE <- TCGA_proj$CPE*100

tmp_pheno <- TCGA_proj[which(TCGA_proj$Tumor_Type2 == "GBM"),]
ggp1 <- ggplot(tmp_pheno, aes(CPE, Tumor, color=high_infil)) + geom_point()
ggp1 <- ggp1 + stat_smooth(method = "lm", geom = "smooth")
ggp1 <- ggp1 + geom_abline(slope=1,intercept=0,lty=2)
ggp1 <- ggp1 + theme_classic()+ theme(legend.position = "none")
ggp1 <- ggp1 + xlab("TCGA CPE (%)") + ylab("GIMiCC Tumor Purity (%)")
ggp1 <- ggp1 + scale_color_manual(values = c("Yes" = "red", "No" = "black"))
ggp1 <- ggp1 + labs(subtitle="GBM")
ggp1

tmp_pheno <- TCGA_proj[which(TCGA_proj$Tumor_Type2 == "OLG"),]
ggp2 <- ggplot(tmp_pheno, aes(CPE, Tumor, color=high_infil)) + geom_point()
ggp2 <- ggp2 + stat_smooth(method = "lm", geom = "smooth")
ggp2 <- ggp2 + geom_abline(slope=1,intercept=0,lty=2)
ggp2 <- ggp2 + theme_classic()+ theme(legend.position = "none")
ggp2 <- ggp2 + xlab("TCGA CPE (%)") + ylab("GIMiCC Tumor Purity (%)")
ggp2 <- ggp2 + scale_color_manual(values = c("Yes" = "red", "No" = "black"))
ggp2 <- ggp2 + labs(subtitle="OLG")
ggp2

tmp_pheno <- TCGA_proj[which(TCGA_proj$Tumor_Type2 == "AST-HG"),]
ggp3 <- ggplot(tmp_pheno, aes(CPE, Tumor, color=high_infil)) + geom_point()
ggp3 <- ggp3 + stat_smooth(method = "lm", geom = "smooth")
ggp3 <- ggp3 + geom_abline(slope=1,intercept=0,lty=2)
ggp3 <- ggp3 + theme_classic()+ theme(legend.position = "none")
ggp3 <- ggp3 + xlab("TCGA CPE (%)") + ylab("GIMiCC Tumor Purity (%)")
ggp3 <- ggp3 + scale_color_manual(values = c("Yes" = "red", "No" = "black"))
ggp3 <- ggp3 + labs(subtitle="AST-HG")
ggp3

tmp_pheno <- TCGA_proj[which(TCGA_proj$Tumor_Type2 == "AST-LG"),]
ggp4  <- ggplot(tmp_pheno, aes(CPE, Tumor, color=high_infil)) + geom_point()
ggp4  <- ggp4  + stat_smooth(method = "lm", geom = "smooth")
ggp4  <- ggp4  + geom_abline(slope=1,intercept=0,lty=2)
ggp4  <- ggp4  + theme_classic() + theme(legend.position = "none")
ggp4  <- ggp4  + xlab("TCGA CPE (%)") + ylab("GIMiCC Tumor Purity (%)")
ggp4  <- ggp4  + scale_color_manual(values = c("Yes" = "red", "No" = "black"))
ggp4 <- ggp4 + labs(subtitle="AST-LG")
ggp4

pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Manuscript/GIMiCC/GIMiCC_October2023/CPE_correlation.pdf", width = 3, height = 3)
ggp1
ggp2
ggp3
ggp4
dev.off()

TCGA_proj$very_diff <- ifelse(TCGA_proj$CPE - TCGA_proj$Tumor > 25, "Yes", "No")

tmp_pheno <- TCGA_proj[which(TCGA_proj$Tumor_Type2 == "GBM"),]
ggp1 <- ggplot(tmp_pheno, aes(CPE, Tumor, color=very_diff )) + geom_point()
ggp1 <- ggp1 + stat_smooth(method = "lm", geom = "smooth")
ggp1 <- ggp1 + geom_abline(slope=1,intercept=0,lty=2)
ggp1 <- ggp1 + theme_classic()+ theme(legend.position = "none")
ggp1 <- ggp1 + xlab("TCGA CPE (%)") + ylab("GIMiCC Tumor Purity (%)")
ggp1 <- ggp1 + scale_color_manual(values = c("Yes" = "red", "No" = "black"))
ggp1 <- ggp1 + labs(subtitle="GBM")


tmp_pheno <- TCGA_proj[which(TCGA_proj$Tumor_Type2 == "OLG"),]
ggp2 <- ggplot(tmp_pheno, aes(CPE, Tumor, color=very_diff )) + geom_point()
ggp2 <- ggp2 + stat_smooth(method = "lm", geom = "smooth")
ggp2 <- ggp2 + geom_abline(slope=1,intercept=0,lty=2)
ggp2 <- ggp2 + theme_classic()+ theme(legend.position = "none")
ggp2 <- ggp2 + xlab("TCGA CPE (%)") + ylab("GIMiCC Tumor Purity (%)")
ggp2 <- ggp2 + scale_color_manual(values = c("Yes" = "red", "No" = "black"))
ggp2 <- ggp2 + labs(subtitle="OLG")


tmp_pheno <- TCGA_proj[which(TCGA_proj$Tumor_Type2 == "AST-HG"),]
ggp3 <- ggplot(tmp_pheno, aes(CPE, Tumor, color=very_diff )) + geom_point()
ggp3 <- ggp3 + stat_smooth(method = "lm", geom = "smooth")
ggp3 <- ggp3 + geom_abline(slope=1,intercept=0,lty=2)
ggp3 <- ggp3 + theme_classic()+ theme(legend.position = "none")
ggp3 <- ggp3 + xlab("TCGA CPE (%)") + ylab("GIMiCC Tumor Purity (%)")
ggp3 <- ggp3 + scale_color_manual(values = c("Yes" = "red", "No" = "black"))
ggp3 <- ggp3 + labs(subtitle="AST-HG")


tmp_pheno <- TCGA_proj[which(TCGA_proj$Tumor_Type2 == "AST-LG"),]
ggp4  <- ggplot(tmp_pheno, aes(CPE, Tumor, color=very_diff )) + geom_point()
ggp4  <- ggp4  + stat_smooth(method = "lm", geom = "smooth")
ggp4  <- ggp4  + geom_abline(slope=1,intercept=0,lty=2)
ggp4  <- ggp4  + theme_classic() + theme(legend.position = "none")
ggp4  <- ggp4  + xlab("TCGA CPE (%)") + ylab("GIMiCC Tumor Purity (%)")
ggp4  <- ggp4  + scale_color_manual(values = c("Yes" = "red", "No" = "black"))
ggp4 <- ggp4 + labs(subtitle="AST-LG")

pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Manuscript/GIMiCC/GIMiCC_October2023/CPE_correlation_highlight_very_different.pdf", width = 3, height = 3)
ggp1
ggp2
ggp3
ggp4
dev.off()


very_different_samples <- TCGA_proj$methylID[TCGA_proj$high_infil == "Yes"]

TCGA_proj <- rbind(GBM_proj, OLG_proj, AST_HG_proj, AST_LG_proj)
TCGA_proj <- TCGA_proj[which(TCGA_proj$methylID %in% very_different_samples),]
TCGA_proj$high_infil <- ifelse(TCGA_proj$Mono + TCGA_proj$CD8nv + TCGA_proj$CD4mem + TCGA_proj$CD4nv + TCGA_proj$CD4mem + TCGA_proj$Treg + TCGA_proj$Bnv + TCGA_proj$Bmem +  TCGA_proj$Neu > 20, "Yes", "No")
TCGA_proj$high_infil_sum <- TCGA_proj$Mono + TCGA_proj$CD8nv + TCGA_proj$CD4mem + TCGA_proj$CD4nv + TCGA_proj$CD4mem + TCGA_proj$Treg + TCGA_proj$Bnv + TCGA_proj$Bmem +  TCGA_proj$Neu 
TCGA_proj$CPE <- TCGA_proj$CPE*100
TCGA_proj$CPE_diff <- TCGA_proj$CPE - TCGA_proj$Tumor 
TCGA_proj$Tumor_Type2 <- factor(TCGA_proj$Tumor_Type2,
                               levels = c("OLG", "AST-LG", "AST-HG", "GBM"))
TCGA_proj$Library <- factor(TCGA_proj$Library,
                             levels = c("OLG", "AST-LG", "AST-HG", "GBM"))


df <- TCGA_proj[, c("Tumor_Type2","Library","CPE_diff")]
tmp_colors <- c(`OLG`="#3D7636",
                `AST-LG`="#9ECAEC",
                `AST-HG`="#342186",
                `GBM`="#BA6A78")
plot1 <- ggplot(df, aes(x=Tumor_Type2, y=CPE_diff, fill=Library)) + geom_boxplot(outlier.shape=NA)
plot1 <- plot1 + theme_classic()
plot1 <- plot1 + xlab("Tumor Sample Type")
plot1 <- plot1 + ylab("CPE Purity - GIMiCC Purity")
plot1 <- plot1 + labs(fill="L0 Library Used")
plot1 <- plot1 + geom_point(position=position_jitterdodge(), shape=1, alpha = 0.3)
plot1 <- plot1 + scale_fill_manual(values = tmp_colors)
plot1 <- plot1 + theme(legend.position="top")

pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Manuscript/GIMiCC/GIMiCC_October2023/CPE_correlation_highlight_very_different.pdf", width = 3, height = 3)
ggp1
ggp2
ggp3
ggp4
dev.off()

pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Manuscript/GIMiCC/GIMiCC_October2023/CPE_correlation_changingL0_highinfil.pdf", width = 6, height = 3)
plot1
dev.off()

df <- TCGA_proj[TCGA_proj$Library == TCGA_proj$Tumor_Type2,]

table(df$Tumor_Type2, df$`Random Forest Sturm Cluster`)
table(df$Tumor_Type2, df$`RPPA cluster`)
table(df$Tumor_Type2, df$`Supervised DNA Methylation Cluster`)
table(df$Tumor_Type2, df$`Pan-Glioma DNA Methylation Cluster`)

############## pannel 6F
library(dplyr)
library(readr)
library(readxl)
library(ggplot2)
library(reshape2)
cell.types <- c("Endothelial","Stromal","Astrocyte","Microglia","Mono","Oligodendrocyte","GABA","GLU","Tumor",
                "Neu","Bmem","Bnv","CD4mem","CD4nv","CD8mem","CD8nv","Treg","NK")


load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Brain_Tumor_Deconv_Ref/Processed_data/TCGA_ref/TCGAall_pheno.RDA")
TCGA_Supp_Info <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/TCGA_Supp_Info.xlsx")
TCGA_pheno <- left_join(TCGA_pheno, TCGA_Supp_Info, by="submitter_id")

df <- read_tsv(file = "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/TCGA_extra_info.tsv")
df <- df[c(5,6,7),]
df <- t(df)
colnames(df) <- df[1,]
df <- data.frame(df[-c(1,2),])
df$submitter_id <- rownames(df)
TCGA_pheno <- left_join(TCGA_pheno, df, by="submitter_id")

TCGA_pheno$Tumor_Type2 <- ifelse(TCGA_pheno$`IDH/codel subtype.y` == "IDHwt", "GBM",
                                 ifelse(TCGA_pheno$`IDH/codel subtype.y`== "IDHmut-codel", "OLG", 
                                        ifelse(TCGA_pheno$`IDH/codel subtype.y`== "IDHmut-non-codel", "AST", "N/A")))

TCGA_pheno$`TERT promoter status.y` <- ifelse(is.na(TCGA_pheno$`TERT promoter status.y`), "N/A", TCGA_pheno$`TERT promoter status.y`)
TCGA_pheno$`Chr 7 gain/Chr 10 loss` <- ifelse(is.na(TCGA_pheno$`Chr 7 gain/Chr 10 loss`), "N/A", TCGA_pheno$`TERT promoter status.y`)
TCGA_pheno$EGFR <- ifelse(is.na(TCGA_pheno$EGFR), "N/A", TCGA_pheno$EGFR)
TCGA_pheno$CDKN2A <- ifelse(is.na(TCGA_pheno$CDKN2A), "N/A", TCGA_pheno$CDKN2A)
TCGA_pheno$CDKN2B <- ifelse(is.na(TCGA_pheno$CDKN2B), "N/A", TCGA_pheno$CDKN2B)

TCGA_pheno$Tumor_Type2 <- ifelse(TCGA_pheno$Tumor_Type2 != "AST", TCGA_pheno$Tumor_Type2, 
                                 ifelse(TCGA_pheno$`TERT promoter status.y`== "Mutant" | 
                                          TCGA_pheno$`Chr 7 gain/Chr 10 loss`== "Gain chr 7 & loss chr 10" |
                                          TCGA_pheno$EGFR == "amp_rec" |
                                          TCGA_pheno$CDKN2A== "homdel_rec" |
                                          TCGA_pheno$CDKN2B== "homdel_rec" , "AST-HG", "AST-LG" ))


load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/TCGA_deconvo_results.RDA")
GBM_proj <- 100*GBM_proj/rowSums(GBM_proj)
OLG_proj <- 100*OLG_proj/rowSums(OLG_proj)
AST_HG_proj <- 100*AST_HG_proj/rowSums(AST_HG_proj)
AST_LG_proj <- 100*AST_LG_proj/rowSums(AST_LG_proj)

GBM_proj$methylID <- rownames(GBM_proj)
OLG_proj$methylID <- rownames(OLG_proj)
AST_HG_proj$methylID <- rownames(AST_HG_proj)
AST_LG_proj$methylID <- rownames(AST_LG_proj)

GBM_proj <- left_join(GBM_proj, TCGA_pheno, by="methylID")
OLG_proj <- left_join(OLG_proj, TCGA_pheno, by="methylID")
AST_HG_proj <- left_join(AST_HG_proj, TCGA_pheno, by="methylID")
AST_LG_proj <- left_join(AST_LG_proj, TCGA_pheno, by="methylID")

GBM_proj$Library <- "GBM"
OLG_proj$Library <- "OLG"
AST_HG_proj$Library <- "AST-HG"
AST_LG_proj$Library <- "AST-LG"

TCGA_proj <- rbind(GBM_proj, OLG_proj, AST_HG_proj, AST_LG_proj)
TCGA_proj <- TCGA_proj[TCGA_proj$Library == TCGA_proj$Tumor_Type2,]
TCGA_proj$high_infil <- ifelse(TCGA_proj$Mono + TCGA_proj$CD8nv + TCGA_proj$CD4mem + TCGA_proj$CD4nv + TCGA_proj$CD4mem + TCGA_proj$Treg + TCGA_proj$Bnv + TCGA_proj$Bmem +  TCGA_proj$Neu > 20, "Yes", "No")
TCGA_proj$high_infil <- ifelse(is.na(TCGA_proj$high_infil), "No", TCGA_proj$high_infil)
TCGA_proj$high_infil_sum <- TCGA_proj$Mono + TCGA_proj$CD8nv + TCGA_proj$CD4mem + TCGA_proj$CD4nv + TCGA_proj$CD4mem + TCGA_proj$Treg + TCGA_proj$Bnv + TCGA_proj$Bmem +  TCGA_proj$Neu 
TCGA_proj$CPE <- TCGA_proj$CPE*100
TCGA_proj$very_diff <- ifelse(TCGA_proj$CPE - TCGA_proj$Tumor > 25, "Yes", "No")
TCGA_proj$Tumor_Type2 <- factor(TCGA_proj$Tumor_Type2,
                                levels = c("OLG", "AST-LG", "AST-HG", "GBM"))
tmp_colors <- c(`OLG`="#3D7636",
                `AST-LG`="#9ECAEC",
                `AST-HG`="#342186",
                `GBM`="#BA6A78")
tmp_colors2 <- c(`Yes`="brown2",
                 `No`="black")

TCGA_proj$EST_IMM <- as.numeric(TCGA_proj$`ESTIMATE immune score`)

plot1 <- ggplot(TCGA_proj, aes(x=Tumor_Type2, 
                               y=EST_IMM, 
                               color=high_infil,
                               fill=Tumor_Type2)) + geom_boxplot(outlier.shape = NA)
plot1 <- plot1 + theme_classic()
plot1 <- plot1 + scale_fill_manual(values = tmp_colors)
plot1 <- plot1 + scale_color_manual(values = tmp_colors2)
plot1 <- plot1 + ylab("ESTIMATE Immune Score")
plot1 <- plot1 + theme(axis.title.x = element_blank(),
                       legend.position = "none")
plot1 <- plot1 + geom_point(position=position_jitterdodge(), shape=1, alpha = 0.3)


pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Manuscript/GIMiCC/GIMiCC_October2023/Estimate_score.pdf", width = 3.5, height = 3)
plot1
dev.off()




############# Figure 7
### Overlap Table
summ_df <- data.frame(matrix(nrow=1,ncol=4))
colnames(summ_df) <- c("Model","Comparison","n_hyper","n_hypo")
counter <- 1
my_combos_to_test <- c("Unadj",
                       "TumorOnly",
                       "ManyBroad",
                       'ManyDetailedv1',
                       'ManyDetailedv2')

library(readxl)
setwd("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/CapperEWAS")
for (tt in c("ASTHG_vs_ASTLG",
             "ASTHG_vs_GBM",
             "ASTHG_vs_OLG")){
  for (tmodel in my_combos_to_test){
    df <- read_xlsx(paste0(tt,"_EWASResults.xlsx"),
                         sheet=tmodel)
    summ_df[counter, "Model"] <- tmodel
    summ_df[counter, "Comparison"] <- tt
    summ_df[counter, "n_hyper"] <- sum((df$delta_beta > 0.3) & (df$adj.P.Val < 0.05))
    summ_df[counter, "n_hypo"] <- sum((df$delta_beta < -0.3) & (df$adj.P.Val < 0.05))
    print(counter)
    counter <- counter+1
  }
}

write.csv(summ_df, 
     file="Tumor_Type_EWASresults.csv",
     row.names = FALSE)



summ_df <- data.frame(matrix(nrow=1,ncol=4))
colnames(summ_df) <- c("Model","Comparison","n_hyper","n_hypo")
counter <- 1
my_combos_to_test <- c("Unadj",
                       "TumorOnly",
                       "ManyBroad",
                       'ManyDetailedv1',
                       'ManyDetailedv2')
setwd("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/CapperEWAS")
for (tt in c("ASTLG_hivsloinfil",
             "ASTHG_hivsloinfil",
             "OLG_hivsloinfil",
             "GBM_hivsloinfil")){
  for (tmodel in my_combos_to_test){
    df <- read_xlsx(paste0(tt,"_EWASResults.xlsx"),
                    sheet=tmodel)
    summ_df[counter, "Model"] <- tmodel
    summ_df[counter, "Comparison"] <- tt
    summ_df[counter, "n_hyper"] <- sum((df$delta_beta > 0.3) & (df$adj.P.Val < 0.05), na.rm=TRUE)
    summ_df[counter, "n_hypo"] <- sum((df$delta_beta < -0.3) & (df$adj.P.Val < 0.05), na.rm=TRUE)
    print(counter)
    counter <- counter+1
  }
}

write.csv(summ_df, 
          file="Hi_Infil_EWASresults.csv",
          row.names = FALSE)

######## local
library(readxl)
library(ggplot2)
library(dplyr)
library(ggpubr)
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/CapperEWAS")


for (tt in c("ASTHG_vs_GBM", "ASTHG_vs_OLG", "ASTHG_vs_ASTLG")){
  print(tt)
  for (tmodel in c("Unadj","ManyDetailedv2")){
    print(tmodel)
    df <- read_xlsx(paste0(tt,"_EWASResults.xlsx"),
                    sheet=tmodel)
    df$tmpcat <- ifelse(df$adj.P.Val<0.05, "Sig","No Difference")
    df$tmpcat <- ifelse(df$tmpcat == "No Difference", "No Difference",
                                  ifelse(df$delta_beta  < -0.3, "Hypomethylated",
                                         ifelse(df$delta_beta  > 0.3, "Hypermethylated","No Difference")))
    df$log10pVal <- -log10(df$adj.P.Val)
    tmp <- max(df$log10pVal[df$log10pVal != Inf &
                              df$log10pVal != -Inf ], na.rm=TRUE)
    tmp2 <- min(df$log10pVal[df$log10pVal != Inf &
                              df$log10pVal != -Inf ], na.rm=TRUE)
    df$log10pVal[df$log10pVal == Inf] <- tmp
    df$log10pVal[df$log10pVal == -Inf] <- tmp2
    
    hyper_hypo_colors <- c(`No Difference` = "gray50", Hypomethylated = "yellow", Hypermethylated = "blue")
    
    p1 <- ggplot(df, aes(delta_beta, log10pVal)) +  geom_point(aes(fill = tmpcat),size = 9/5,shape=21, show.legend = FALSE) 
    p1 <- p1 + xlab(expression(paste(Delta,Beta)))  
    p1 <- p1 + ylab(expression("-log"[10]*"(Padj)")) + theme_classic() 
    p1 <- p1 + geom_hline(yintercept = -log10(0.05), color = "red", size = 0.8, linetype="dashed") 
    p1 <- p1 + geom_vline(xintercept = -0.3, color = "red", size = 0.8, linetype="dashed") 
    p1 <- p1 + geom_vline(xintercept = 0.3, color = "red", size = 0.8, linetype="dashed") 
    p1 <- p1 + scale_fill_manual(values = hyper_hypo_colors)
    
    pdf(paste0(tt,"_",tmodel,"_volcanoplot.pdf"), height = 3, width = 3)
    print(p1)
    dev.off()

  }
}

for (tt in c("ASTHG_hivsloinfil", "ASTLG_hivsloinfil",
             "GBM_hivsloinfil", "OLG_hivsloinfil")){
  print(tt)
  for (tmodel in c("Unadj","ManyDetailedv2")){
    print(tmodel)
    df <- read_xlsx(paste0(tt,"_EWASResults.xlsx"),
                    sheet=tmodel)
    df$tmpcat <- ifelse(df$adj.P.Val<0.05, "Sig","No Difference")
    df$tmpcat <- ifelse(df$tmpcat == "No Difference", "No Difference",
                        ifelse(df$delta_beta  < -0.3, "Hypomethylated",
                               ifelse(df$delta_beta  > 0.3, "Hypermethylated","No Difference")))
    df$log10pVal <- -log10(df$adj.P.Val)
    tmp <- max(df$log10pVal[df$log10pVal != Inf &
                              df$log10pVal != -Inf ], na.rm=TRUE)
    tmp2 <- min(df$log10pVal[df$log10pVal != Inf &
                               df$log10pVal != -Inf ], na.rm=TRUE)
    df$log10pVal[df$log10pVal == Inf] <- tmp
    df$log10pVal[df$log10pVal == -Inf] <- tmp2
    
    hyper_hypo_colors <- c(`No Difference` = "gray50", Hypomethylated = "yellow", Hypermethylated = "blue")
    
    p1 <- ggplot(df, aes(delta_beta, log10pVal)) +  geom_point(aes(fill = tmpcat),size = 9/5,shape=21, show.legend = FALSE) 
    p1 <- p1 + xlab(expression(paste(Delta,Beta)))  
    p1 <- p1 + ylab(expression("-log"[10]*"(Padj)")) + theme_classic() 
    p1 <- p1 + geom_hline(yintercept = -log10(0.05), color = "red", size = 0.8, linetype="dashed") 
    p1 <- p1 + geom_vline(xintercept = -0.3, color = "red", size = 0.8, linetype="dashed") 
    p1 <- p1 + geom_vline(xintercept = 0.3, color = "red", size = 0.8, linetype="dashed") 
    p1 <- p1 + scale_fill_manual(values = hyper_hypo_colors)
    
    pdf(paste0(tt,"_",tmodel,"_volcanoplot.pdf"), height = 3, width = 3)
    print(p1)
    dev.off()
    
  }
}