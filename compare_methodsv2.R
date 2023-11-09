library(minfi)
library(dplyr)

load("/dartfs/rc/lab/S/SalasLab/Brain_Tumor_Deconv_Ref/Processed_data/Immune_cells/FlowSorted.BloodExtended.EPIC.RData")
immune_pheno<-as.data.frame(pData(FlowSorted.BloodExtended.EPIC))
samples_keep <- which(immune_pheno$Sample_Group %in% c("MIX_base","MIX_Treg","MIX_basophil"))
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

immune_raw_beta <- immune_raw_beta[complete.cases(immune_raw_beta), ]

library(ChAMP)
betas <- champ.norm(beta=immune_raw_beta,method="BMIQ") #306,935

save(betas, file ="/dartfs/rc/lab/S/SalasLab/Brain_Tumor_Deconv_Ref/Processed_data/Immune_cells/MIX_betas.RDA")

MIX_pheno <- immune_pheno
save(MIX_pheno, file ="/dartfs/rc/lab/S/SalasLab/Brain_Tumor_Deconv_Ref/Processed_data/Immune_cells/MIX_pheno.RDA")

########################

load("/dartfs/rc/lab/S/SalasLab/Brain_Tumor_Deconv_Ref/Processed_data/Immune_cells/MIX_pheno.RDA")
load("/dartfs/rc/lab/S/SalasLab/Brain_Tumor_Deconv_Ref/Processed_data/Immune_cells/MIX_betas.RDA")

library(MethylCIBERSORT)
library(EpiDISH)

data("StromalMatrix_V2")
Int <- intersect(rownames(Mat), rownames(Stromal_v2))
Signature <- FeatureSelect.V4(CellLines.matrix = NULL, Heatmap = FALSE, export = TRUE, sigName = "ExampleType",
                              Stroma.matrix = Stromal_v2, deltaBeta = 0.2, FDR = 0.01, MaxDMRs = 100,
                              Phenotype.stroma = Stromal_v2.pheno)
data("V2_Signatures")

Signaturex <- Signatures$Glioma_v2_Signature.txt
Probes = rownames(Signature$SignatureMatrix)

brain_sig <- Signatures$Glioma_v2_Signature.txt
brain_sig1 <- brain_sig[,-1]
rownames(brain_sig1) <- brain_sig[,1]
brain_sig1<-brain_sig1/100
brain_sig1<-as.matrix(brain_sig1)

cpgs <- intersect(rownames(brain_sig1), rownames(betas))
betas1<-betas[cpgs,]
brain_sig1 <- brain_sig1[cpgs,]
betas1<-as.matrix(betas1)
betas2<-as.matrix(betas1[complete.cases(betas1), ])

Glioma_Prediction<-epidish(
  beta.m = betas2,
  ref.m = brain_sig1,
  method = "CBS",
)
Glioma_Prediction2<-as.data.frame(Glioma_Prediction$estF)
setwd("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/")
save(Glioma_Prediction2, file ="MIX_decon_CIBERSORT.rda")

########################
load("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets/DeconvoFunct.RDATA")


library(FlowSorted.Blood.EPIC)
gimmic_predictions_GBM <- DeconvoFunct(betas, h=5, tumor.type = "GBM")
gimmic_predictions_OLG <- DeconvoFunct(betas, h=5, tumor.type = "OLG")
gimmic_predictions_ASTHG <- DeconvoFunct(betas, h=5, tumor.type = "AST-HG")
gimmic_predictions_ASTLG <- DeconvoFunct(betas, h=5, tumor.type = "AST-LG")

setwd("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/")
save(gimmic_predictions_GBM,
     gimmic_predictions_OLG ,
     gimmic_predictions_ASTHG,
     gimmic_predictions_ASTLG, file ="MIX_decon_GIMiCC.rda")

########################

library(MethylResolver)
resolver_predictions <- MethylResolver(methylMix = betas, betaPrime = FALSE)
setwd("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/")
resolver_predictions <- read.delim("MethylResolver.txt")
save(resolver_predictions , file ="MIX_decon_MethylResolver.rda")


##################### load results locally for plot
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/")
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Brain_Tumor_Deconv_Ref/Processed_data/Immune_cells/MIX_pheno.RDA")
load("MIX_decon_MethylResolver.rda")
load("MIX_decon_CIBERSORT.rda")
load("MIX_decon_GIMiCC.rda")

df <- data.frame(matrix(NA,nrow=12,ncol=5))
rownames(df) <- rownames(MIX_pheno)
colnames(df) <- c("Gran","NK","Mono","Tcell","Bcell")

df_GIMMIC_OLG <- df
df_GIMMIC_ASTLG <- df
df_GIMMIC_ASTHG <- df
df_GIMMIC_GBM <- df
df_MethylResolver <- df
df_CIBERSORT <- df

df_GIMMIC_OLG$Gran <- (MIX_pheno$Gran*100) - (gimmic_predictions_OLG$Neu)
df_GIMMIC_ASTLG$Gran <- (MIX_pheno$Gran*100) - (gimmic_predictions_ASTLG$Neu)
df_GIMMIC_ASTHG$Gran <- (MIX_pheno$Gran*100) - (gimmic_predictions_ASTHG$Neu)
df_GIMMIC_GBM$Gran <- (MIX_pheno$Gran*100) - (gimmic_predictions_GBM$Neu)
df_MethylResolver$Gran <- (MIX_pheno$Gran*100) - (100*(resolver_predictions$Neu + resolver_predictions$Eos))
df_CIBERSORT$Gran <- (MIX_pheno$Gran*100) - (100*(Glioma_Prediction2$Neu + Glioma_Prediction2$Eos))

df_GIMMIC_OLG$NK <- (MIX_pheno$NK*100) - (gimmic_predictions_OLG$NK)
df_GIMMIC_ASTLG$NK <- (MIX_pheno$NK*100) - (gimmic_predictions_ASTLG$NK)
df_GIMMIC_ASTHG$NK <- (MIX_pheno$NK*100) - (gimmic_predictions_ASTHG$NK)
df_GIMMIC_GBM$NK <- (MIX_pheno$NK*100) - (gimmic_predictions_GBM$NK)
df_MethylResolver$NK <- (MIX_pheno$NK*100) - (100*resolver_predictions$NK)
df_CIBERSORT$NK <- (MIX_pheno$NK*100) - (100*Glioma_Prediction2$CD56)

df_GIMMIC_OLG$Mono <- (MIX_pheno$Mono*100) - (gimmic_predictions_OLG$Mono)
df_GIMMIC_ASTLG$Mono <- (MIX_pheno$Mono*100) - (gimmic_predictions_ASTLG$Mono)
df_GIMMIC_ASTHG$Mono <- (MIX_pheno$Mono*100) - (gimmic_predictions_ASTHG$Mono)
df_GIMMIC_GBM$Mono <- (MIX_pheno$Mono*100) - (gimmic_predictions_GBM$Mono)
df_MethylResolver$Mono <- (MIX_pheno$Mono*100) - (100*(resolver_predictions$Mon + resolver_predictions$Macro + resolver_predictions$Dendritic))
df_CIBERSORT$Mono <- (MIX_pheno$Mono*100) - (100*Glioma_Prediction2$CD14)

df_GIMMIC_OLG$Tcell <- ((MIX_pheno$CD4T+MIX_pheno$CD8T)*100) - ((gimmic_predictions_OLG$CD4nv + gimmic_predictions_OLG$CD8nv + gimmic_predictions_OLG$Treg + gimmic_predictions_OLG$CD8mem + gimmic_predictions_OLG$CD4mem))
df_GIMMIC_ASTLG$Tcell <- ((MIX_pheno$CD4T+MIX_pheno$CD8T)*100) - ((gimmic_predictions_ASTLG$CD4nv + gimmic_predictions_ASTLG$CD8nv + gimmic_predictions_ASTLG$Treg + gimmic_predictions_ASTLG$CD8mem + gimmic_predictions_ASTLG$CD4mem))
df_GIMMIC_ASTHG$Tcell <- ((MIX_pheno$CD4T+MIX_pheno$CD8T)*100) - ((gimmic_predictions_ASTHG$CD4nv + gimmic_predictions_ASTHG$CD8nv + gimmic_predictions_ASTHG$Treg + gimmic_predictions_ASTHG$CD8mem + gimmic_predictions_ASTHG$CD4mem))
df_GIMMIC_GBM$Tcell <- ((MIX_pheno$CD4T+MIX_pheno$CD8T)*100) - ((gimmic_predictions_GBM$CD4nv + gimmic_predictions_GBM$CD8nv + gimmic_predictions_GBM$Treg + gimmic_predictions_GBM$CD8mem + gimmic_predictions_GBM$CD4mem))
df_MethylResolver$Tcell <- ((MIX_pheno$CD4T+MIX_pheno$CD8T)*100) - (100*(resolver_predictions$Treg + resolver_predictions$Tnaive + resolver_predictions$Tmem + resolver_predictions$CD8))
df_CIBERSORT$Tcell <- ((MIX_pheno$CD4T+MIX_pheno$CD8T)*100) - (100*(Glioma_Prediction2$CD4_Eff + Glioma_Prediction2$CD8 + Glioma_Prediction2$Treg + Glioma_Prediction2$CD8))

df_GIMMIC_OLG$Bcell <- (MIX_pheno$Bcell*100) - ((gimmic_predictions_OLG$Bmem + gimmic_predictions_OLG$Bnv))
df_GIMMIC_ASTLG$Bcell <- (MIX_pheno$Bcell*100) - ((gimmic_predictions_ASTLG$Bmem + gimmic_predictions_ASTLG$Bnv))
df_GIMMIC_ASTHG$Bcell <- (MIX_pheno$Bcell*100) - ((gimmic_predictions_ASTHG$Bmem + gimmic_predictions_ASTHG$Bnv))
df_GIMMIC_GBM$Bcell <- (MIX_pheno$Bcell*100) - ((gimmic_predictions_GBM$Bmem + gimmic_predictions_GBM$Bnv))
df_MethylResolver$Bcell <- (MIX_pheno$Bcell*100) - (100*resolver_predictions$Bcell)
df_CIBERSORT$Bcell <- (MIX_pheno$Bcell*100) - (100*Glioma_Prediction2$CD19)

df_GIMMIC_OLG$Method <- "GIMiCC OLG"
df_GIMMIC_ASTHG$Method <- "GIMiCC AST-HG"
df_GIMMIC_ASTLG$Method <- "GIMiCC AST-LG"
df_GIMMIC_GBM$Method <- "GIMiCC GBM"
df_MethylResolver$Method <- "MethylResolver"
df_CIBERSORT$Method <- "methylCIBERSORT"

tmp_colors <- c(`GIMiCC OLG`="#3D7636",
                `GIMiCC AST-LG`="#9ECAEC",
                `GIMiCC AST-HG`="#342186",
                `GIMiCC GBM`="#BA6A78",
                MethylResolver="#EBE63C",
                methylCIBERSORT="#CA5F03")

df <- rbind(df_GIMMIC_OLG,df_GIMMIC_GBM,
            df_GIMMIC_ASTHG, df_GIMMIC_ASTLG,
            df_MethylResolver,df_CIBERSORT)
library(reshape2)
df2 <- melt(df, id.vars = "Method")

library(ggplot2)
p1 <- ggplot(df2, aes(x=variable, y=value, fill=Method)) + 
      geom_boxplot() +
      theme_classic() +
      geom_hline(yintercept = 0, lty=2) +
      xlab("Cell Type")+
      ylab("Real - Predicited") +
      scale_fill_manual(values=tmp_colors)

setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/")
pdf("compare_to_othermethods.pdf", width = 6, height = 3)
p1
dev.off()
