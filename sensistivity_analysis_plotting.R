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
dir <- "/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/L0_sensitivity_analysis/"

CpG_lists <- list()
Deconvo_results <- list()
used_in_construction <- list()

for (n in 1:n_folds){
  message("Fold #: ", n)
  
  load(file = paste0(dir, "L0_Fold",n,".RDA"))
  CpG_lists[[n]] <- rownames(Library_Layer0)
  used_in_construction[[n]] <- colnames(Library_Layer0)[1:(ncol(Library_Layer0)-1)]
  
  tumor_sample <- colnames(agg_betas)
  tumor_beta_iDMC <- agg_betas[rownames(agg_betas) %in% rownames(Library_Layer0), ]
  idmc_dat <- Library_Layer0[rownames(tumor_beta_iDMC), ]
  purity <- c()
  
  for (t in tumor_sample) {
    beta_adj <- c(
      tumor_beta_iDMC[idmc_dat$hyper == TRUE, t],
      1 - tumor_beta_iDMC[idmc_dat$hyper == FALSE, t]
    )
    pu <- InfiniumPurify:::.get_peak(beta_adj)
    purity[t] <- pu
  }
  purity_iDMC <- as.data.frame(purity)
  proj <- purity_iDMC
  proj$NonTumor <- 1 - proj$purity
  colnames(proj)[1] <- "Tumor"
  Deconvo_results[[n]] <- proj
  
  names(Deconvo_results)[n] <- paste0("L0_Fold",n)
  names(used_in_construction)[n] <- paste0("L0_Fold",n)
  names(CpG_lists)[n] <- paste0("L0_Fold",n)
  
}

setwd(dir)
save(Deconvo_results, used_in_construction, CpG_lists,
     file="L0_Sensitivity_Analysis.RDATA")



##################### below is local
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/L0_sensitivity_analysis/L0_Sensitivity_Analysis.RDATA")
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/L0_sensitivity_analysis")
library(UpSetR)

load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Capper/Processed_data/Capper_pheno_agg2.RDA")
Ctrl_pheno <- Capper_pheno[which(Capper_pheno$Broad_Class == "Control_Healthy"),]
Case_pheno <- Capper_pheno[which(Capper_pheno$methylation.class.ch1 == "A IDH"),]
total_pheno <- rbind(Ctrl_pheno, Case_pheno)

pdf("CpG_Selection_Overlapp.pdf", width=6, height=5.5)
upset(fromList(CpG_lists), nsets=10, order.by = "freq")
dev.off()

for (n in 1:10){
  if (n==1){
    df <- Deconvo_results[[n]]
    df$fold <- n
    df$Train_data <- ifelse(rownames(df) %in% used_in_construction[[n]], 1, 0)
    df$ID <- rownames(Deconvo_results[[n]])
    } else {
    tmpdf <- Deconvo_results[[n]]
    tmpdf$fold <- n
    tmpdf$Train_data <- ifelse(rownames(tmpdf) %in% used_in_construction[[n]], 1, 0)
    tmpdf$ID <- rownames(Deconvo_results[[n]])
    df <- rbind(df, tmpdf)
  }
}

df$Dx <- ifelse(df$ID %in% rownames(Ctrl_pheno), "HC", "Tumor")

df$fold <- as.factor(df$fold)
train <- df[df$Train_data == 1,]
test <- df[df$Train_data == 0,]

library(ggplot2)
library(ggpubr)

tmp_colors <- c("HC" = "#E5D1D6",
                "Tumor" = "#CD4022")

plot1 <- ggplot(test, aes(x=fold, y=Tumor, fill=Dx)) + geom_boxplot()
plot1 <- plot1 + geom_hline(yintercept=mean(test$Tumor[test$Dx=="Tumor"]))
plot1 <- plot1 + geom_hline(yintercept=mean(test$Tumor[test$Dx=="HC"]))
plot1 <- plot1 + geom_vline(xintercept=1.5, lty="dotted", alpha=0.5)
plot1 <- plot1 + geom_vline(xintercept=2.5, lty="dotted")
plot1 <- plot1 + geom_vline(xintercept=3.5, lty="dotted")
plot1 <- plot1 + geom_vline(xintercept=4.5, lty="dotted")
plot1 <- plot1 + geom_vline(xintercept=5.5, lty="dotted")
plot1 <- plot1 + geom_vline(xintercept=6.5, lty="dotted")
plot1 <- plot1 + geom_vline(xintercept=7.5, lty="dotted")
plot1 <- plot1 + geom_vline(xintercept=8.5, lty="dotted")
plot1 <- plot1 + geom_vline(xintercept=9.5, lty="dotted")
plot1 <- plot1 + theme_classic()
plot1 <- plot1 + scale_fill_manual(values=tmp_colors)
plot1 <- plot1 + xlab("Fold")
plot1 <- plot1 + ylab("Tumor Fraction")
plot1 <- plot1 + labs(subtitle="Test Data")

plot2 <- ggplot(train, aes(x=fold, y=Tumor, fill=Dx)) + geom_boxplot()
plot2 <- plot2 + geom_hline(yintercept=mean(train$Tumor[train$Dx=="Tumor"]))
plot2 <- plot2 + geom_hline(yintercept=mean(train$Tumor[train$Dx=="HC"]))
plot2 <- plot2 + geom_vline(xintercept=1.5, lty="dotted", alpha=0.5)
plot2 <- plot2 + geom_vline(xintercept=2.5, lty="dotted")
plot2 <- plot2 + geom_vline(xintercept=3.5, lty="dotted")
plot2 <- plot2 + geom_vline(xintercept=4.5, lty="dotted")
plot2 <- plot2 + geom_vline(xintercept=5.5, lty="dotted")
plot2 <- plot2 + geom_vline(xintercept=6.5, lty="dotted")
plot2 <- plot2 + geom_vline(xintercept=7.5, lty="dotted")
plot2 <- plot2 + geom_vline(xintercept=8.5, lty="dotted")
plot2 <- plot2 + geom_vline(xintercept=9.5, lty="dotted")
plot2 <- plot2 + theme_classic()
plot2 <- plot2 + scale_fill_manual(values=tmp_colors)
plot2 <- plot2 + xlab("Fold")
plot2 <- plot2 + ylab("Tumor Fraction")
plot2 <- plot2 + labs(subtitle="Training Data")



pdf("L0_prediction_by_Fold.pdf", width=6, height=3, onefile = TRUE)
plot1
plot2
dev.off()