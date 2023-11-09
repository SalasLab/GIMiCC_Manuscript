library(FlowSorted.Blood.EPIC)
library(minfi)
library(tibble)
library(dplyr)
library(InfiniumPurify)
library(ggplot2)
library(reshape2)
setwd("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/")

cell.types <- c("Endothelial","Stromal","Astrocyte","Microglia","Mono","Oligodendrocyte","GABA","GLU","Tumor",
                "Neu","Bmem","Bnv","CD4mem","CD4nv","CD8mem","CD8nv","Treg","NK")

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

load("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/DeconvoFunct.RDATA")

##### TCGA Data
setwd("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/")
betas <- readRDS("/dartfs/rc/lab/S/SalasLab/Brain_Tumor_Deconv_Ref/Processed_data/TCGA_ref/TCGAall_cbindfill.RDS")
load("/dartfs/rc/lab/S/SalasLab/Brain_Tumor_Deconv_Ref/Processed_data/TCGA_ref/TCGAall_pheno.RDA")
samps <- intersect(colnames(betas), TCGA_pheno$methylID)

betas <- betas[, which(colnames(betas) %in% samps)]

GBM_proj <- DeconvoFunct(betas, tumor.type = "GBM", h=5)
AST_HG_proj <- DeconvoFunct(betas, tumor.type = "AST-HG", h=5)
AST_LG_proj <- DeconvoFunct(betas, tumor.type = "AST-LG", h=5)
OLG_proj <- DeconvoFunct(betas, tumor.type = "OLG", h=5)

setwd("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/")
save(GBM_proj,AST_HG_proj,AST_LG_proj, OLG_proj, file="TCGA_deconvo_results.RDA")

##### Capper Data
rm(list=ls())
cell.types <- c("Endothelial","Stromal","Astrocyte","Microglia","Mono","Oligodendrocyte","GABA","GLU","Tumor",
                "Neu","Bmem","Bnv","CD4mem","CD4nv","CD8mem","CD8nv","Treg","NK")

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

load("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/DeconvoFunct.RDATA")

load("/dartfs/rc/lab/S/SalasLab/Capper/Processed_data/Capper_pheno_agg2.RDA")
rownames(Capper_pheno) <- paste0(Capper_pheno$Slide,"_",Capper_pheno$Array)

base_dir <- "/dartfs/rc/lab/S/SalasLab/Capper/Processed_data/"
first <- TRUE
for (i in 1:40){
  print(i)
  set_name <- paste0(base_dir,"Capper_batch",i,"/Capper_batch",i,"_betas_filtered.RDA")
  betas <- loadRData(set_name)
  tmpGBM_proj <- DeconvoFunct(betas, tumor.type = "GBM", h=5)
  tmpAST_HG_proj <- DeconvoFunct(betas, tumor.type = "AST-HG", h=5)
  tmpAST_LG_proj <- DeconvoFunct(betas, tumor.type = "AST-LG", h=5)
  tmpOLG_proj <- DeconvoFunct(betas, tumor.type = "OLG", h=5)
  if (first == TRUE){
    GBM_proj <- tmpGBM_proj
    AST_HG_proj <- tmpAST_HG_proj
    AST_LG_proj <- tmpAST_LG_proj
    OLG_proj <- tmpOLG_proj
    first <- FALSE
  } else {
    GBM_proj <- rbind(GBM_proj, tmpGBM_proj)
    AST_HG_proj <- rbind(AST_HG_proj, tmpAST_HG_proj)
    AST_LG_proj <- rbind(AST_LG_proj, tmpAST_LG_proj)
    OLG_proj <- rbind(OLG_proj, tmpOLG_proj)
  }
}


setwd("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/")
save(GBM_proj,AST_HG_proj,AST_LG_proj, OLG_proj, file="Capper_deconvo_results.RDA")

