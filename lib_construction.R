# Method 11
# this is a repeat of method 9 with full immune compartment
library(ddpcr)
library(limma)
library(minfi)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(Metrics)
library(tibble)
library(openxlsx)
library(tidyr)
library(splitstackshape)
library(FlowSorted.Blood.EPIC)

############################################# 
rm(list=ls())
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Brain_deconvolution/Raw data/Reboot_09072022/newMeffiFunction.RDATA")

load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/April2023_Analysis/library_construction/BMIQ_Normalized2.RDATA")

load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Brain_deconvolution/Raw data/Reboot_09072022/Microglial_GSE191200_control_idats/microglia_pheno.RDATA")
microglia_pheno$ID <- paste0(microglia_pheno$Sample.ID,"_" , microglia_pheno$Sentrix_ID,"_" ,microglia_pheno$Sentrix_Position)
rmv <- microglia_pheno$ID[microglia_pheno$Age > 70]

dim(new_pheno)
dim(new_raw_beta_normalized)
new_pheno <- new_pheno[-which(new_pheno$ID %in% rmv),]
new_raw_beta_normalized <- new_raw_beta_normalized[, which(colnames(new_raw_beta_normalized) %in% new_pheno$ID)]
identical(colnames(new_raw_beta_normalized), new_pheno$ID)
dim(new_pheno)
dim(new_raw_beta_normalized)

#---------------------------------------------------------------------
betas1<-as.matrix(new_raw_beta_normalized)
new_pheno$Layer1<-ifelse(new_pheno$CellType %in% c("Astrocyte","Microglia","Oligodendrocyte"), "Glial", 
                         ifelse(new_pheno$CellType %in% c("Endothelial","Stromal"), "Endothelial and Stromal",
                            ifelse(new_pheno$CellType %in% c("GABA","GLU"),"Neuronal","Immune")))

table(new_pheno$CellType,new_pheno$Layer1)
cell.type<-new_pheno$Layer1

DMRt<-meffil.cell.type.specific.methylation.by.t(betas1,cell.type, number.sites = 25)
reference_library_two_50<-as.matrix(DMRt[[1]])
projection_two_50<-projectCellType_CP(betas1[rownames(reference_library_two_50),], 
                                      as.matrix(reference_library_two_50),lessThanOne =T)*100


Library_Layer1<-reference_library_two_50
#---------------------------------------------------------------------
new_pheno$Layer2A<-ifelse(new_pheno$CellType %in% c("Astrocyte","Microglia","Oligodendrocyte"), "Glial", 
                          ifelse(new_pheno$CellType %in% c("Endothelial","Stromal"), new_pheno$CellType,
                                 ifelse(new_pheno$CellType %in% c("GABA","GLU"),"Neuronal","Immune")))

table(new_pheno$CellType,new_pheno$Layer2A)
cell.type<-new_pheno$Layer2A

DMRt<-meffil.cell.type.specific.methylation.by.t(betas1,cell.type, number.sites = 25)
reference_library_two_50<-as.matrix(DMRt[[1]])
projection_two_50<-projectCellType_CP(betas1[rownames(reference_library_two_50),], 
                                      as.matrix(reference_library_two_50),lessThanOne =T)*100


Library_Layer2A<-reference_library_two_50

#---------------------------------------------------------------------
new_pheno$Layer2B<-ifelse(new_pheno$CellType %in% c("Astrocyte","Microglia","Oligodendrocyte"), new_pheno$CellType, 
                          ifelse(new_pheno$CellType %in% c("Endothelial","Stromal"), "Endothelial and Stromal",
                                 ifelse(new_pheno$CellType %in% c("GABA","GLU"),"Neuronal","Immune")))

table(new_pheno$CellType,new_pheno$Layer2B)
cell.type<-new_pheno$Layer2B

DMRt<-meffil.cell.type.specific.methylation.by.t(betas1,cell.type, number.sites = 25)
reference_library_two_50<-as.matrix(DMRt[[1]])
projection_two_50<-projectCellType_CP(betas1[rownames(reference_library_two_50),], 
                                      as.matrix(reference_library_two_50),lessThanOne =T)*100


Library_Layer2B<-reference_library_two_50


#---------------------------------------------------------------------
new_pheno$Layer2C<-ifelse(new_pheno$CellType %in% c("Astrocyte","Microglia","Oligodendrocyte"), "Glial", 
                          ifelse(new_pheno$CellType %in% c("Endothelial","Stromal"), "Endothelial and Stromal",
                                 ifelse(new_pheno$CellType %in% c("GABA","GLU"),new_pheno$CellType,"Immune")))

table(new_pheno$CellType,new_pheno$Layer2C)
cell.type<-new_pheno$Layer2C

DMRt<-meffil.cell.type.specific.methylation.by.t(betas1,cell.type, number.sites = 25)
reference_library_two_50<-as.matrix(DMRt[[1]])
projection_two_50<-projectCellType_CP(betas1[rownames(reference_library_two_50),], 
                                      as.matrix(reference_library_two_50),lessThanOne =T)*100


Library_Layer2C<-reference_library_two_50
#---------------------------------------------------------------------
new_pheno$Layer2D<-ifelse(new_pheno$CellType %in% c("Astrocyte","Microglia","Oligodendrocyte"), "Glial", 
                          ifelse(new_pheno$CellType %in% c("Endothelial","Stromal"), "Endothelial and Stromal",
                                 ifelse(new_pheno$CellType %in% c("GABA","GLU"),"Neuronal",
                                        ifelse(new_pheno$CellType %in% c("Neu","Mono"),"Myeloid","Lymphoid"))))

table(new_pheno$CellType,new_pheno$Layer2D)
cell.type<-new_pheno$Layer2D

DMRt<-meffil.cell.type.specific.methylation.by.t(betas1,cell.type, number.sites = 25)
reference_library_two_50<-as.matrix(DMRt[[1]])
projection_two_50<-projectCellType_CP(betas1[rownames(reference_library_two_50),], 
                                      as.matrix(reference_library_two_50),lessThanOne =T)*100


Library_Layer2D<-reference_library_two_50

#---------------------------------------------------------------------
new_pheno$Layer3A<-ifelse(new_pheno$CellType %in% c("Astrocyte","Microglia","Oligodendrocyte"), "Glial", 
                          ifelse(new_pheno$CellType %in% c("Endothelial","Stromal"), "Endothelial and Stromal",
                                 ifelse(new_pheno$CellType %in% c("GABA","GLU"),"Neuronal",
                                        ifelse(new_pheno$CellType %in% c("Neu","Mono"),"Myeloid",
                                               ifelse(new_pheno$CellType %in% c("Bnv","Bmem"),"Bcell",
                                                      ifelse(new_pheno$CellType %in% c("NK"),"NK","Tcell"))))))

table(new_pheno$CellType,new_pheno$Layer3A)
cell.type<-new_pheno$Layer3A

DMRt<-meffil.cell.type.specific.methylation.by.t(betas1,cell.type, number.sites = 25)
reference_library_two_50<-as.matrix(DMRt[[1]])
projection_two_50<-projectCellType_CP(betas1[rownames(reference_library_two_50),], 
                                      as.matrix(reference_library_two_50),lessThanOne =T)*100


Library_Layer3A<-reference_library_two_50

#---------------------------------------------------------------------
new_pheno$Layer3B<-ifelse(new_pheno$CellType %in% c("Astrocyte","Microglia","Oligodendrocyte"), "Glial", 
                          ifelse(new_pheno$CellType %in% c("Endothelial","Stromal"), "Endothelial and Stromal",
                                 ifelse(new_pheno$CellType %in% c("GABA","GLU"),"Neuronal",
                                        ifelse(new_pheno$CellType %in% c("Neu","Mono"),new_pheno$CellType,"Lymphoid"))))

table(new_pheno$CellType,new_pheno$Layer3B)
cell.type<-new_pheno$Layer3B

DMRt<-meffil.cell.type.specific.methylation.by.t(betas1,cell.type, number.sites = 25)
reference_library_two_50<-as.matrix(DMRt[[1]])
projection_two_50<-projectCellType_CP(betas1[rownames(reference_library_two_50),], 
                                      as.matrix(reference_library_two_50),lessThanOne =T)*100


Library_Layer3B<-reference_library_two_50

#---------------------------------------------------------------------
new_pheno$Layer4<-ifelse(new_pheno$CellType %in% c("Astrocyte","Microglia","Oligodendrocyte"), "Glial", 
                          ifelse(new_pheno$CellType %in% c("Endothelial","Stromal"), "Endothelial and Stromal",
                                 ifelse(new_pheno$CellType %in% c("GABA","GLU"),"Neuronal",
                                        ifelse(new_pheno$CellType %in% c("Neu","Mono"),"Myeloid",
                                               ifelse(new_pheno$CellType %in% c("Bnv","Bmem"),"Bcell",
                                                      ifelse(new_pheno$CellType %in% c("NK"),"NK",
                                                             ifelse(new_pheno$CellType %in% c("CD8mem","CD8nv"),"CD8Tcell","CD4Tcell")))))))

table(new_pheno$CellType,new_pheno$Layer4)
cell.type<-new_pheno$Layer4

DMRt<-meffil.cell.type.specific.methylation.by.t(betas1,cell.type, number.sites = 25)
reference_library_two_50<-as.matrix(DMRt[[1]])
projection_two_50<-projectCellType_CP(betas1[rownames(reference_library_two_50),], 
                                      as.matrix(reference_library_two_50),lessThanOne =T)*100


Library_Layer4<-reference_library_two_50

#---------------------------------------------------------------------
new_pheno$Layer5A<-ifelse(new_pheno$CellType %in% c("Astrocyte","Microglia","Oligodendrocyte"), "Glial", 
                         ifelse(new_pheno$CellType %in% c("Endothelial","Stromal"), "Endothelial and Stromal",
                                ifelse(new_pheno$CellType %in% c("GABA","GLU"),"Neuronal",
                                       ifelse(new_pheno$CellType %in% c("Neu","Mono"),"Myeloid",
                                              ifelse(new_pheno$CellType %in% c("Bnv","Bmem"),"Bcell",
                                                     ifelse(new_pheno$CellType %in% c("NK"),"NK",
                                                            ifelse(new_pheno$CellType %in% c("CD8mem","CD8nv"),"CD8Tcell",new_pheno$CellType)))))))

table(new_pheno$CellType,new_pheno$Layer5A)
cell.type<-new_pheno$Layer5A

DMRt<-meffil.cell.type.specific.methylation.by.t(betas1,cell.type, number.sites = 25)
reference_library_two_50<-as.matrix(DMRt[[1]])
projection_two_50<-projectCellType_CP(betas1[rownames(reference_library_two_50),], 
                                      as.matrix(reference_library_two_50),lessThanOne =T)*100


Library_Layer5A<-reference_library_two_50

#---------------------------------------------------------------------
new_pheno$Layer5B<-ifelse(new_pheno$CellType %in% c("Astrocyte","Microglia","Oligodendrocyte"), "Glial", 
                          ifelse(new_pheno$CellType %in% c("Endothelial","Stromal"), "Endothelial and Stromal",
                                 ifelse(new_pheno$CellType %in% c("GABA","GLU"),"Neuronal",
                                        ifelse(new_pheno$CellType %in% c("Neu","Mono"),"Myeloid",
                                               ifelse(new_pheno$CellType %in% c("Bnv","Bmem"),"Bcell",
                                                      ifelse(new_pheno$CellType %in% c("NK"),"NK",
                                                             ifelse(new_pheno$CellType %in% c("CD8mem","CD8nv"),new_pheno$CellType,"CD4Tcell")))))))

table(new_pheno$CellType,new_pheno$Layer5B)
cell.type<-new_pheno$Layer5B

DMRt<-meffil.cell.type.specific.methylation.by.t(betas1,cell.type, number.sites = 25)
reference_library_two_50<-as.matrix(DMRt[[1]])
projection_two_50<-projectCellType_CP(betas1[rownames(reference_library_two_50),], 
                                      as.matrix(reference_library_two_50),lessThanOne =T)*100


Library_Layer5B<-reference_library_two_50

#---------------------------------------------------------------------
new_pheno$Layer5C<-ifelse(new_pheno$CellType %in% c("Astrocyte","Microglia","Oligodendrocyte"), "Glial", 
                          ifelse(new_pheno$CellType %in% c("Endothelial","Stromal"), "Endothelial and Stromal",
                                 ifelse(new_pheno$CellType %in% c("GABA","GLU"),"Neuronal",
                                        ifelse(new_pheno$CellType %in% c("Neu","Mono"),"Myeloid",
                                               ifelse(new_pheno$CellType %in% c("Bnv","Bmem"),new_pheno$CellType,
                                                      ifelse(new_pheno$CellType %in% c("NK"),"NK",
                                                             ifelse(new_pheno$CellType %in% c("CD8mem","CD8nv"),"CD8Tcell","CD4Tcell")))))))

table(new_pheno$CellType,new_pheno$Layer5C)
cell.type<-new_pheno$Layer5C

DMRt<-meffil.cell.type.specific.methylation.by.t(betas1,cell.type, number.sites = 25)
reference_library_two_50<-as.matrix(DMRt[[1]])
projection_two_50<-projectCellType_CP(betas1[rownames(reference_library_two_50),], 
                                      as.matrix(reference_library_two_50),lessThanOne =T)*100


Library_Layer5C<-reference_library_two_50

load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/GBM_Library0.RDA")
GBM_Layer0 <- Library_Layer0
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/OLI_Library0.RDA")
OLI_Layer0 <- Library_Layer0
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/AST_HG_Library0.RDA")
AST_HG_Layer0 <- Library_Layer0
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction/AST_LG_Library0.RDA")
AST_LG_Layer0 <- Library_Layer0


#---------------------------------------------------------------------
library(UpSetR)
probe_list <- list(rownames(GBM_Layer0),
                   rownames(OLI_Layer0),
                   rownames(AST_HG_Layer0),
                   rownames(AST_LG_Layer0),
                   rownames(Library_Layer1),
                   rownames(Library_Layer2A),
                   rownames(Library_Layer2B),
                   rownames(Library_Layer2C),
                   rownames(Library_Layer2D),
                   rownames(Library_Layer3A),
                   rownames(Library_Layer3B),
                   rownames(Library_Layer4),
                   rownames(Library_Layer5A),
                   rownames(Library_Layer5B),
                   rownames(Library_Layer5C))
names(probe_list) <- c("Layer0_GBM", "Layer0_OLI", "Layer0_AST_HG", "Layer0_AST_LG",
                       "Layer1", "Layer2A", "Layer2B", "Layer2C","Layer2D",
                       "Layer3A","Layer3B", "Layer4", "Layer5A", "Layer5B","Layer5C")

setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_construction")

pdf("CpG_overlapp.pdf")
UpSetR::upset(UpSetR::fromList(probe_list), order.by = "freq",nsets=15)
dev.off()

setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/library_sets")
Library_Layer0 <- GBM_Layer0
save(Library_Layer0,Library_Layer1,Library_Layer2A,Library_Layer2B,Library_Layer2C, 
     Library_Layer2D,Library_Layer3A,Library_Layer3B,Library_Layer4,Library_Layer5A,
     Library_Layer5B,Library_Layer5C,file = "GBM_Libraries.RDATA")
Library_Layer0 <- OLI_Layer0
save(Library_Layer0,Library_Layer1,Library_Layer2A,Library_Layer2B,Library_Layer2C, 
     Library_Layer2D,Library_Layer3A,Library_Layer3B,Library_Layer4,Library_Layer5A,
     Library_Layer5B,Library_Layer5C,file = "OLI_Libraries.RDATA")
Library_Layer0 <- AST_HG_Layer0
save(Library_Layer0,Library_Layer1,Library_Layer2A,Library_Layer2B,Library_Layer2C, 
     Library_Layer2D,Library_Layer3A,Library_Layer3B,Library_Layer4,Library_Layer5A,
     Library_Layer5B,Library_Layer5C,file = "AST_HG_Libraries.RDATA")
Library_Layer0 <- AST_LG_Layer0
save(Library_Layer0,Library_Layer1,Library_Layer2A,Library_Layer2B,Library_Layer2C, 
     Library_Layer2D,Library_Layer3A,Library_Layer3B,Library_Layer4,Library_Layer5A,
     Library_Layer5B,Library_Layer5C,file = "AST_LG_Libraries.RDATA")


library(pheatmap)
pdf("Layer_heatmaps.pdf", onefile = TRUE)

pheatmap(
  GBM_Layer0[,1:(dim(GBM_Layer0)[2]-1)],
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = F, #samples
  #labels_col = rownames(TC_Library),
  #annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

pheatmap(
  OLI_Layer0[,1:(dim(OLI_Layer0)[2]-1)],
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = F, #samples
  #labels_col = rownames(TC_Library),
  #annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

pheatmap(
  AST_HG_Layer0[,1:(dim(AST_HG_Layer0)[2]-1)],
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = F, #samples
  #labels_col = rownames(TC_Library),
  #annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

pheatmap(
  AST_LG_Layer0[,1:(dim(AST_LG_Layer0)[2]-1)],
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = F, #samples
  #labels_col = rownames(TC_Library),
  #annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

pheatmap(
  Library_Layer1,
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = T, #samples
  #labels_col = rownames(TC_Library),
  #annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

pheatmap(
  Library_Layer2A,
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = T, #samples
  #labels_col = rownames(TC_Library),
  #annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

pheatmap(
  Library_Layer2B,
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = T, #samples
  #labels_col = rownames(TC_Library),
  #annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

pheatmap(
  Library_Layer2C,
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = T, #samples
  #labels_col = rownames(TC_Library),
  #annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

pheatmap(
  Library_Layer2D,
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = T, #samples
  #labels_col = rownames(TC_Library),
  #annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

pheatmap(
  Library_Layer3A,
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = T, #samples
  #labels_col = rownames(TC_Library),
  #annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

pheatmap(
  Library_Layer3B,
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = T, #samples
  #labels_col = rownames(TC_Library),
  #annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

pheatmap(
  Library_Layer4,
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = T, #samples
  #labels_col = rownames(TC_Library),
  #annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

pheatmap(
  Library_Layer5A,
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = T, #samples
  #labels_col = rownames(TC_Library),
  #annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

pheatmap(
  Library_Layer5B,
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = T, #samples
  #labels_col = rownames(TC_Library),
  #annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

pheatmap(
  Library_Layer5C,
  annotation_names_col = F,
  show_rownames = FALSE, #CpGs 
  show_colnames = T, #samples
  #labels_col = rownames(TC_Library),
  #annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

dev.off()





