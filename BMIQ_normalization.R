
#unzip---------------------------------------------- Steve did not rerun this chunk
library(R.utils)
# idatFiles <- list.files("Astrocyte_GSE166845_idats/", pattern = "idat.gz$", full = TRUE)
# sapply(idatFiles, gunzip, overwrite = TRUE)

setwd("/dartfs/rc/lab/S/SalasLab/Brain deconvolution/Raw data/Reboot_09072022")
library(tibble)
library(dplyr)
library(minfi)
#-----------------------------------------------
workdir <- "Astrocyte_GSE166845_idats/"
Astrocyte_RGset <- read.metharray.exp(workdir,force = TRUE)
Astrocyte_dp<-minfi::detectionP(Astrocyte_RGset)
sum(Astrocyte_dp>0.01)
Astrocyte_Raw_betas <- preprocessRaw(Astrocyte_RGset)
Astrocyte_Raw_betas<-getBeta(Astrocyte_Raw_betas)
identical(rownames(Astrocyte_Raw_betas),rownames(Astrocyte_dp))
identical(colnames(Astrocyte_Raw_betas),colnames(Astrocyte_dp))
for (i in 1:nrow(Astrocyte_Raw_betas)) {
  for (j in 1:ncol(Astrocyte_Raw_betas)) {
    Astrocyte_Raw_betas[i,j]<-ifelse(Astrocyte_dp[i,j]<0.01,Astrocyte_Raw_betas[i,j],NA)
  }
  
}
#----------------------------------------------------------------------------
workdir <- "Microglial_GSE191200_control_idats/"
Microglial_RGset <- read.metharray.exp(workdir,force = TRUE)
Microglial_dp<-minfi::detectionP(Microglial_RGset)
sum(Microglial_dp>0.01)
Microglial_Raw_betas <- preprocessRaw(Microglial_RGset)
Microglial_Raw_betas<-getBeta(Microglial_Raw_betas)
identical(rownames(Microglial_Raw_betas),rownames(Microglial_dp))
identical(colnames(Microglial_Raw_betas),colnames(Microglial_dp))
for (i in 1:nrow(Microglial_Raw_betas)) {
  for (j in 1:ncol(Microglial_Raw_betas)) {
    Microglial_Raw_betas[i,j]<-ifelse(Microglial_dp[i,j]<0.01,Microglial_Raw_betas[i,j],NA)
  }
  
}
#--------------------------------------------------------------------------------

load("/dartfs/rc/lab/S/SalasLab/Brain deconvolution/Raw data/CNS_reference.RData") #Steve had to change path
cns_pheno<-as.data.frame(cns_reference@colData)
cns_pheno <- rownames_to_column(as.data.frame(cns_pheno))
cns_keep_pheno<-cns_pheno %>% filter(CellType == "Endothelial" | CellType == "Stromal" | CellType == "GABA" | CellType == "GLU")

cns_reference_keep<-cns_reference[,cns_keep_pheno$rowname]
cns_Raw <- preprocessRaw(cns_reference_keep)
cns_raw_beta<-as.data.frame(getBeta(cns_Raw))
#------------------------------------------------------------------------------------------
load("/dartfs/rc/lab/S/SalasLab/Brain_Tumor_Deconv_Ref/Raw_data/Oligodendrocytes.RData") #Steve had to change path
Oligo_pheno<-phenoclean %>% filter(characteristics_ch1.2 == "cell type: Olig2")
Oligo_beta<-as.data.frame(Oligo[,colnames(Oligo) %in% Oligo_pheno$geo_accession])
for (i in 1:ncol(Oligo_beta)) {
    Oligo_beta[,i]<-ifelse(is.nan(Oligo_beta[,i]),NA,Oligo_beta[,i])
}
#---------------------------------------------------------


#####---------------------------------- Steve Added Leukocytes
load("/dartfs/rc/lab/S/SalasLab/Brain_Tumor_Deconv_Ref/Processed_data/Immune_cells/FlowSorted.BloodExtended.EPIC.RData")
immune_pheno<-as.data.frame(pData(FlowSorted.BloodExtended.EPIC))
samples_keep <- which(immune_pheno$Sample_Group %in% c("Bmem","Bnv","CD4subset","CD8subset","Tregnew","FlowSorted.Blood.EPIC"))
immune_RGSet_raw <- FlowSorted.BloodExtended.EPIC[,samples_keep]
immune_pheno<-as.data.frame(pData(immune_RGSet_raw))
rm(FlowSorted.BloodExtended.EPIC, samples_keep)
class(immune_RGSet_raw)

immune_dp<-minfi::detectionP(immune_RGSet_raw)
sum(immune_dp>0.01)
immune_raw_beta <- getBeta(immune_RGSet_raw)
identical(rownames(immune_raw_beta),rownames(immune_dp))
identical(colnames(immune_raw_beta),colnames(immune_dp))
for (i in 1:nrow(immune_raw_beta)) {
  for (j in 1:ncol(immune_raw_beta)) {
    immune_raw_beta[i,j]<-ifelse(immune_dp[i,j]<0.01,immune_raw_beta[i,j],NA)
  }
  
}
#####---------------------------------- 

#####---------------------------------- Steve Changed the manifest file because I did not have access to these 

# load("C:/Users/F004304/Dropbox (Dartmouth College)/Annotation/EPIC.hg19.manifest.RDATA")
# load("C:/Users/F004304/Dropbox (Dartmouth College)/Annotation/HM450.hg19.manifest.RDATA")
# 
# use1<-HM450_annotation %>% filter(MASK_general == "FALSE" & NonCpG == "FALSE" & SexProbe == "FALSE")
# use1<-use1$probeID
# 
# use2<-annotation %>% filter(MASK_general == "FALSE" & NonCpG == "FALSE" & SexProbe == "FALSE")
# use2<-use2$probeID
# use<-intersect(use1,use2)

use <- intersect(rownames(Astrocyte_Raw_betas), rownames(Microglial_Raw_betas))
use <- intersect(use, rownames(cns_raw_beta))
use <- intersect(use, rownames(Oligo_beta))
use <- intersect(use, rownames(immune_raw_beta))



load("/dartfs/rc/lab/S/SalasLab/R/Annotation files/annotationEPICb5updated3.rda")
rm(annotDF)
annot$probeType <- substr(annot$Name, 1, 2)

bad_rows <- rownames(annot[annot$CHR=="Y" | 
                             annot$CHR=="X" | 
                             annot$probeType=="rs" |
                             annot$probeType=="ch" |
                             annot$mask==TRUE,])
use <- use[-which(use %in% bad_rows)] # 345,063 CpGs
rm(annot)
#####---------------------------------- 


Astrocyte_Raw_betas<-Astrocyte_Raw_betas[use,]
Microglial_Raw_betas<-Microglial_Raw_betas[use,]
immune_raw_beta<-immune_raw_beta[use,] ##### Steve
cns_raw_beta<-cns_raw_beta[use,]

identical(rownames(immune_raw_beta), rownames(Microglial_Raw_betas)) ##### Steve double checking 

new_raw_beta<-cbind(Astrocyte_Raw_betas,Microglial_Raw_betas)
new_raw_beta<-cbind(new_raw_beta,immune_raw_beta) ##### Steve
new_raw_beta<-as.data.frame(cbind(new_raw_beta,cns_raw_beta))

Oligo_beta[Oligo_beta==1]<-max(new_raw_beta[new_raw_beta<0.99], na.rm = TRUE)
Oligo_beta[Oligo_beta==0]<-min(new_raw_beta[new_raw_beta>0.01], na.rm = TRUE)
Oligo_beta<-Oligo_beta[use,]
new_raw_beta<-cbind(new_raw_beta,Oligo_beta)

setwd("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/library_construction/")##### Steve

#BiocManager::install("ChAMP")
save(new_raw_beta,file = "/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/library_construction/Pre_BMIQ_Normalized2_betas.RDATA") ##### Steve saving this object before normalization steps

library(ChAMP)
new_raw_beta_update<- new_raw_beta[complete.cases(new_raw_beta), ]
new_raw_beta_normalized<- champ.norm(beta=new_raw_beta_update,method="BMIQ") #306,935

new_pheno<-as.data.frame(matrix(NA,ncol(new_raw_beta_normalized),2))
colnames(new_pheno)<-c("ID","CellType")
new_pheno$ID<-colnames(new_raw_beta_normalized)
new_pheno[1:6,2]<-"Astrocyte"
new_pheno[7:24,2]<-"Microglia"
new_pheno[25:70,2]<-immune_pheno$CellType
new_pheno[71:106,2]<-cns_keep_pheno$CellType
new_pheno[107:126,2]<-"Oligodendrocyte"

save(new_raw_beta_normalized,new_pheno, file = "/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/library_construction/BMIQ_Normalized2.RDATA")
