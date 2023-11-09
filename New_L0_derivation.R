# Method 11
# this is a repeat of method 9 with full immune compartment
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

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

load("/dartfs/rc/lab/S/SalasLab/R/Annotation files/annotationEPICb5updated3.rda")
rm(annotDF)
annot2 <- annot
annot2 <- annot2[annot2$mask == FALSE,]
annot2 <- annot2[annot2$cg == "cg",]
annot2 <- annot2[annot2$CHR_hg38 != "chrX",]
annot2 <- annot2[annot2$CHR_hg38 != "chrY",]
keep_Cpgs <- annot2$Name
rm(annot,annot2)

########### IDH mut Grade 2/3
load("/dartfs/rc/lab/S/SalasLab/Capper/Processed_data/Capper_pheno_agg2.RDA")
rownames(Capper_pheno) <- paste0(Capper_pheno$Slide,"_",Capper_pheno$Array)
Capper_GBM_subset <- Capper_pheno[Capper_pheno$Broad_Class %in% c("Control_Healthy") |
                                  Capper_pheno$methylation.class.ch1 %in% c("A IDH"),]
Capper_GBM_subset <- Capper_GBM_subset[Capper_GBM_subset$Train_split==TRUE,]

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
      CpGs <- intersect(rownames(Capper_GBM_betas), names(betas))
      betas <- betas[CpGs]
      betas <- betas[order(names(betas))]
      Capper_GBM_betas <- Capper_GBM_betas[CpGs,]
      Capper_GBM_betas <- Capper_GBM_betas[order(rownames(Capper_GBM_betas)),]
      print(identical(names(betas), rownames(Capper_GBM_betas)))
    } else {
      CpGs <- intersect(rownames(Capper_GBM_betas), rownames(betas))
      betas <- betas[CpGs,]
      betas <- betas[order(rownames(betas)),]
      Capper_GBM_betas <- Capper_GBM_betas[CpGs,]
      Capper_GBM_betas <- Capper_GBM_betas[order(rownames(Capper_GBM_betas)),]
      print(identical(rownames(betas), rownames(Capper_GBM_betas)))
    }
    Capper_GBM_betas <- cbind(Capper_GBM_betas, betas)
  }
}

Capper_GBM_betas <- Capper_GBM_betas[which(rownames(Capper_GBM_betas) %in% keep_Cpgs),]

samps <- intersect(rownames(Capper_GBM_subset), colnames(Capper_GBM_betas))
Capper_GBM_subset <- Capper_GBM_subset[samps,]
Capper_GBM_betas <- Capper_GBM_betas[,samps]
identical(rownames(Capper_GBM_subset), colnames(Capper_GBM_betas))
tumors <- which(Capper_GBM_subset$Broad_Class != "Control_Healthy")
ctrls <- which(Capper_GBM_subset$Broad_Class == "Control_Healthy")

tumor.data <- Capper_GBM_betas[, tumors]
normal.data <- Capper_GBM_betas[, ctrls]
rm(Capper_GBM_betas)

iDMC2 <- InfiniumPurify:::get_iDMC(tumor.data, normal.data)
idmc.dat <- data.frame(cbind(tumor.data[iDMC2, ],
                             normal.data[iDMC2, ]))
colnames(idmc.dat) <- c(colnames(tumor.data), colnames(normal.data))
idmc.dat$hyper <- rowMeans(idmc.dat[, rownames(Capper_GBM_subset[tumors,])],
                           na.rm = TRUE) > rowMeans(idmc.dat[, rownames(Capper_GBM_subset[ctrls,])],
                                                    na.rm = TRUE)

Library_Layer0<-idmc.dat
save(Library_Layer0, file = "/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/AST_LG_Library0.RDA")


############ IDH mut Grade 4
load("/dartfs/rc/lab/S/SalasLab/Capper/Processed_data/Capper_pheno_agg2.RDA")
rownames(Capper_pheno) <- paste0(Capper_pheno$Slide,"_",Capper_pheno$Array)
Capper_GBM_subset <- Capper_pheno[Capper_pheno$Broad_Class %in% c("Control_Healthy") | 
                                    Capper_pheno$methylation.class.ch1 %in% c("A IDH, HG"),]
Capper_GBM_subset <- Capper_GBM_subset[Capper_GBM_subset$Train_split==TRUE,]

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
      CpGs <- intersect(rownames(Capper_GBM_betas), names(betas))
      betas <- betas[CpGs]
      betas <- betas[order(names(betas))]
      Capper_GBM_betas <- Capper_GBM_betas[CpGs,]
      Capper_GBM_betas <- Capper_GBM_betas[order(rownames(Capper_GBM_betas)),]
      print(identical(names(betas), rownames(Capper_GBM_betas)))
    } else {
      CpGs <- intersect(rownames(Capper_GBM_betas), rownames(betas))
      betas <- betas[CpGs,]
      betas <- betas[order(rownames(betas)),]
      Capper_GBM_betas <- Capper_GBM_betas[CpGs,]
      Capper_GBM_betas <- Capper_GBM_betas[order(rownames(Capper_GBM_betas)),]
      print(identical(rownames(betas), rownames(Capper_GBM_betas)))
    }
    Capper_GBM_betas <- cbind(Capper_GBM_betas, betas)
  }
}

Capper_GBM_betas <- Capper_GBM_betas[which(rownames(Capper_GBM_betas) %in% keep_Cpgs),]

samps <- intersect(rownames(Capper_GBM_subset), colnames(Capper_GBM_betas))
Capper_GBM_subset <- Capper_GBM_subset[samps,]
Capper_GBM_betas <- Capper_GBM_betas[,samps]
identical(rownames(Capper_GBM_subset), colnames(Capper_GBM_betas))
tumors <- which(Capper_GBM_subset$Broad_Class != "Control_Healthy")
ctrls <- which(Capper_GBM_subset$Broad_Class == "Control_Healthy")

tumor.data <- Capper_GBM_betas[, tumors]
normal.data <- Capper_GBM_betas[, ctrls]
rm(Capper_GBM_betas)

iDMC2 <- InfiniumPurify:::get_iDMC(tumor.data, normal.data)
idmc.dat <- data.frame(cbind(tumor.data[iDMC2, ],
                             normal.data[iDMC2, ]))
colnames(idmc.dat) <- c(colnames(tumor.data), colnames(normal.data))
idmc.dat$hyper <- rowMeans(idmc.dat[, rownames(Capper_GBM_subset[tumors,])],
                           na.rm = TRUE) > rowMeans(idmc.dat[, rownames(Capper_GBM_subset[ctrls,])],
                                                    na.rm = TRUE)

Library_Layer0<-idmc.dat
save(Library_Layer0, file = "/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/AST_HG_Library0.RDA")


############ Oligos
load("/dartfs/rc/lab/S/SalasLab/Capper/Processed_data/Capper_pheno_agg2.RDA")
rownames(Capper_pheno) <- paste0(Capper_pheno$Slide,"_",Capper_pheno$Array)
Capper_GBM_subset <- Capper_pheno[Capper_pheno$Broad_Class %in% c("Control_Healthy") | 
                                    Capper_pheno$methylation.class.ch1 %in% c("O IDH"),]
Capper_GBM_subset <- Capper_GBM_subset[Capper_GBM_subset$Train_split==TRUE,]

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
      CpGs <- intersect(rownames(Capper_GBM_betas), names(betas))
      betas <- betas[CpGs]
      betas <- betas[order(names(betas))]
      Capper_GBM_betas <- Capper_GBM_betas[CpGs,]
      Capper_GBM_betas <- Capper_GBM_betas[order(rownames(Capper_GBM_betas)),]
      print(identical(names(betas), rownames(Capper_GBM_betas)))
    } else {
      CpGs <- intersect(rownames(Capper_GBM_betas), rownames(betas))
      betas <- betas[CpGs,]
      betas <- betas[order(rownames(betas)),]
      Capper_GBM_betas <- Capper_GBM_betas[CpGs,]
      Capper_GBM_betas <- Capper_GBM_betas[order(rownames(Capper_GBM_betas)),]
      print(identical(rownames(betas), rownames(Capper_GBM_betas)))
    }
    Capper_GBM_betas <- cbind(Capper_GBM_betas, betas)
  }
}

Capper_GBM_betas <- Capper_GBM_betas[which(rownames(Capper_GBM_betas) %in% keep_Cpgs),]

samps <- intersect(rownames(Capper_GBM_subset), colnames(Capper_GBM_betas))
Capper_GBM_subset <- Capper_GBM_subset[samps,]
Capper_GBM_betas <- Capper_GBM_betas[,samps]
identical(rownames(Capper_GBM_subset), colnames(Capper_GBM_betas))
tumors <- which(Capper_GBM_subset$Broad_Class != "Control_Healthy")
ctrls <- which(Capper_GBM_subset$Broad_Class == "Control_Healthy")

tumor.data <- Capper_GBM_betas[, tumors]
normal.data <- Capper_GBM_betas[, ctrls]
rm(Capper_GBM_betas)

iDMC2 <- InfiniumPurify:::get_iDMC(tumor.data, normal.data)
idmc.dat <- data.frame(cbind(tumor.data[iDMC2, ],
                             normal.data[iDMC2, ]))
colnames(idmc.dat) <- c(colnames(tumor.data), colnames(normal.data))
idmc.dat$hyper <- rowMeans(idmc.dat[, rownames(Capper_GBM_subset[tumors,])],
                           na.rm = TRUE) > rowMeans(idmc.dat[, rownames(Capper_GBM_subset[ctrls,])],
                                                    na.rm = TRUE)

Library_Layer0<-idmc.dat
save(Library_Layer0, file = "/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/OLI_Library0.RDA")


############ GBMs
load("/dartfs/rc/lab/S/SalasLab/Capper/Processed_data/Capper_pheno_agg2.RDA")
rownames(Capper_pheno) <- paste0(Capper_pheno$Slide,"_",Capper_pheno$Array)
Capper_GBM_subset <- Capper_pheno[Capper_pheno$Broad_Class %in% c("Glioblastoma", "Control_Healthy"),]
Capper_GBM_subset <- Capper_GBM_subset[Capper_GBM_subset$Train_split==TRUE,]

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
      CpGs <- intersect(rownames(Capper_GBM_betas), names(betas))
      betas <- betas[CpGs]
      betas <- betas[order(names(betas))]
      Capper_GBM_betas <- Capper_GBM_betas[CpGs,]
      Capper_GBM_betas <- Capper_GBM_betas[order(rownames(Capper_GBM_betas)),]
      print(identical(names(betas), rownames(Capper_GBM_betas)))
    } else {
      CpGs <- intersect(rownames(Capper_GBM_betas), rownames(betas))
      betas <- betas[CpGs,]
      betas <- betas[order(rownames(betas)),]
      Capper_GBM_betas <- Capper_GBM_betas[CpGs,]
      Capper_GBM_betas <- Capper_GBM_betas[order(rownames(Capper_GBM_betas)),]
      print(identical(rownames(betas), rownames(Capper_GBM_betas)))
    }
    Capper_GBM_betas <- cbind(Capper_GBM_betas, betas)
  }
}

Capper_GBM_betas <- Capper_GBM_betas[which(rownames(Capper_GBM_betas) %in% keep_Cpgs),]

samps <- intersect(rownames(Capper_GBM_subset), colnames(Capper_GBM_betas))
Capper_GBM_subset <- Capper_GBM_subset[samps,]
Capper_GBM_betas <- Capper_GBM_betas[,samps]
identical(rownames(Capper_GBM_subset), colnames(Capper_GBM_betas))
tumors <- which(Capper_GBM_subset$Broad_Class != "Control_Healthy")
ctrls <- which(Capper_GBM_subset$Broad_Class == "Control_Healthy")

tumor.data <- Capper_GBM_betas[, tumors]
normal.data <- Capper_GBM_betas[, ctrls]
rm(Capper_GBM_betas)

iDMC2 <- InfiniumPurify:::get_iDMC(tumor.data, normal.data)
idmc.dat <- data.frame(cbind(tumor.data[iDMC2, ],
                             normal.data[iDMC2, ]))
colnames(idmc.dat) <- c(colnames(tumor.data), colnames(normal.data))
idmc.dat$hyper <- rowMeans(idmc.dat[, rownames(Capper_GBM_subset[tumors,])],
                           na.rm = TRUE) > rowMeans(idmc.dat[, rownames(Capper_GBM_subset[ctrls,])],
                                                    na.rm = TRUE)

Library_Layer0<-idmc.dat
save(Library_Layer0, file = "/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/GBM_Library0.RDA")

