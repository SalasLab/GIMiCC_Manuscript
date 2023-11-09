# set.seed(11)
# library(caTools)
# 
# ############# splitting the Capper et al dataset
# load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Capper/Processed_data/Capper_pheno_agg.RDA")
# rownames(Capper_pheno) <- paste0(Capper_pheno$Slide,"_",Capper_pheno$Array)
# Capper_pheno$ID <- rownames(Capper_pheno)
# Capper_pheno$Train_split <- NA
# 
# tumor_types <- c("O IDH", "A IDH", "A IDH, HG",
#                  "GBM, G34", "GBM, MES", "GBM, MID",
#                  "GBM, MYCN", "GBM, RTK I", "GBM, RTK II",
#                  "GBM, RTK III", "DMG, K27",
#                  "CONTR, ADENOPIT", "CONTR, CEBM", "CONTR, HEMI",
#                  "CONTR, HYPTHAL", "CONTR, PINEAL", "CONTR, PONS", 
#                  "CONTR, WM")
# for (tumor_type in tumor_types){
#   tmp_subset <- Capper_pheno[Capper_pheno$methylation.class.ch1 == tumor_type,]
#   tmp_subset$Train_split <- sample.split(tmp_subset$ID, SplitRatio = 0.75)
#   for (i in 1:nrow(tmp_subset)){
#     tmp_samp <- tmp_subset$ID[i]
#     tmp_row <- which(Capper_pheno$ID == tmp_samp)
#     Capper_pheno$Train_split[tmp_row] <- tmp_subset$Train_split[i]
#   }
# }
# 
# save(Capper_pheno, 
#      file="//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Capper/Processed_data/Capper_pheno_agg2.RDA")



############# sensitivity analysis on the sample split on L0 dev
set.seed(45)
library(limma)
library(minfi)
library(dplyr)
library(ggplot2)
library(Metrics)
library(tibble)
library(openxlsx)
library(tidyr)
library(FlowSorted.Blood.EPIC)
library(InfiniumPurify)
library(caTools)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

tumor_types <- c("A IDH", "CONTR, ADENOPIT", "CONTR, CEBM", "CONTR, HEMI",
                 "CONTR, HYPTHAL", "CONTR, PINEAL", "CONTR, PONS", 
                 "CONTR, WM")

load("/dartfs/rc/lab/S/SalasLab/Capper/Processed_data/Capper_pheno_agg2.RDA")
Ctrl_pheno <- Capper_pheno[which(Capper_pheno$Broad_Class == "Control_Healthy"),]
Case_pheno <- Capper_pheno[which(Capper_pheno$methylation.class.ch1 == "A IDH"),]
total_pheno <- rbind(Ctrl_pheno, Case_pheno)

#### grab all the samples we will need
base_dir <- "/dartfs/rc/lab/S/SalasLab/Capper/Processed_data/"
for (i in unique(total_pheno$processing_batch)){
  print(i)
  set_name <- paste0(base_dir,"Capper_batch",i,"/Capper_batch",i,"_betas_filtered.RDA")
  betas <- loadRData(set_name)
  betas <- betas[, which(colnames(betas) %in% rownames(total_pheno))]
  if (i == 1){
    agg_betas <- betas
  } else {
    if (is.null(dim(betas))){
      CpGs <- intersect(rownames(agg_betas), names(betas))
      betas <- betas[CpGs]
      betas <- betas[order(names(betas))]
      agg_betas <- agg_betas[CpGs,]
      agg_betas <- agg_betas[order(rownames(agg_betas)),]
      print(identical(names(betas), rownames(agg_betas)))
    } else {
      CpGs <- intersect(rownames(agg_betas), rownames(betas))
      betas <- betas[CpGs,]
      betas <- betas[order(rownames(betas)),]
      agg_betas <- agg_betas[CpGs,]
      agg_betas <- agg_betas[order(rownames(agg_betas)),]
      print(identical(rownames(betas), rownames(agg_betas)))
    }
    agg_betas <- cbind(agg_betas, betas)
  }
}

samps <- intersect(rownames(total_pheno), colnames(agg_betas))
total_pheno <- total_pheno[samps,]
agg_betas <- agg_betas[,samps]
identical(rownames(total_pheno), colnames(agg_betas))


n_folds <- 10

for (n in 1:n_folds){
  message("Fold #: ", n)
  # split the data in a different way each fold
  total_pheno$Train_split <- NA
  for (tumor_type in tumor_types){
    tmp_subset <- total_pheno[total_pheno$methylation.class.ch1 == tumor_type,]
    tmp_subset$Train_split <- sample.split(tmp_subset$ID, SplitRatio = 0.75)
    for (i in 1:nrow(tmp_subset)){
      tmp_samp <- tmp_subset$ID[i]
      tmp_row <- which(total_pheno$ID == tmp_samp)
      total_pheno$Train_split[tmp_row] <- tmp_subset$Train_split[i]
    }
  }
  
  tmp_pheno <- total_pheno[total_pheno$Train_split == TRUE,]
  print(table(tmp_pheno$methylation.class.ch1))
  
  tumors <- rownames(tmp_pheno)[tmp_pheno$Broad_Class != "Control_Healthy"]
  ctrls <- rownames(tmp_pheno)[tmp_pheno$Broad_Class == "Control_Healthy"]
  
  tumor.data <- agg_betas[, tumors]
  normal.data <- agg_betas[, ctrls]

  iDMC2 <- InfiniumPurify:::get_iDMC(tumor.data, normal.data)
  idmc.dat <- data.frame(cbind(tumor.data[iDMC2, ],
                               normal.data[iDMC2, ]))
  colnames(idmc.dat) <- c(colnames(tumor.data), colnames(normal.data))
  idmc.dat$hyper <- rowMeans(idmc.dat[, tumors],
                             na.rm = TRUE) > rowMeans(idmc.dat[, ctrls],
                                                      na.rm = TRUE)
  
  Library_Layer0<-idmc.dat
  
  dir <- "/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/L0_sensitivity_analysis/"
  save(Library_Layer0, 
       file = paste0(dir, "L0_Fold",n,".RDA"))
}





