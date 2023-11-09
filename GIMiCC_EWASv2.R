library(readxl)
library(readr)
library(dplyr)

load("/dartfs/rc/lab/S/SalasLab/PD/Scripts/Test_Deconvo_Func_Output/DeconvoEWAS2.RDATA")

load("/dartfs/rc/lab/S/SalasLab/Capper/Processed_data/Capper_pheno_agg2.RDA")
Capper_Supp_Info <- read_excel("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Capper_Supp_Info.xlsx")
rownames(Capper_pheno) <- paste0(Capper_pheno$Slide,"_",Capper_pheno$Array)
Capper_pheno$ID <- substr(rownames(Capper_pheno), 12, 100)
pheno_df <- left_join(Capper_pheno, Capper_Supp_Info)
rownames(pheno_df) <- paste0(pheno_df$Slide,"_",pheno_df$Array)
pheno_df$ID <- rownames(pheno_df)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


base_dir <- "/dartfs/rc/lab/S/SalasLab/Capper/Processed_data/"
first <- TRUE
for (i in unique(pheno_df$processing_batch)){
  print(i)
  set_name <- paste0(base_dir,"Capper_batch",i,"/Capper_batch",i,"_betas_filtered.RDA")
  betas <- loadRData(set_name)
  betas <- betas[, which(colnames(betas) %in% pheno_df$ID)]
  if (first == TRUE){
    agg_betas <- betas
    first <- FALSE
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


load("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/Capper_deconvo_results.RDA")
GBM_proj[is.na(GBM_proj)] <- 0
OLG_proj[is.na(OLG_proj)] <- 0
AST_HG_proj[is.na(AST_HG_proj)] <- 0
AST_LG_proj[is.na(AST_LG_proj)] <- 0

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

GBM_proj <- GBM_proj[GBM_proj$Broad_Class == "Glioblastoma",]
OLG_proj <- OLG_proj[OLG_proj$methylation.class.ch1 == "O IDH",]
AST_HG_proj <- AST_HG_proj[AST_HG_proj$methylation.class.ch1 == "A IDH, HG",]
AST_LG_proj <- AST_LG_proj[AST_LG_proj$methylation.class.ch1 == "A IDH",]

total_proj <- rbind(GBM_proj, OLG_proj, AST_HG_proj, AST_LG_proj)

library(dplyr)
total_proj$Tcell <- total_proj$CD8nv + total_proj$CD8mem + total_proj$CD4mem + total_proj$CD4nv + total_proj$Treg
total_proj$Bcell <- total_proj$Bnv + total_proj$Bmem
total_proj$Myeloid <- total_proj$Mono + total_proj$Neu
total_proj$Lymphoid <- total_proj$Tcell + total_proj$NK + total_proj$Bcell
total_proj$Immune <- total_proj$Lymphoid + total_proj$Myeloid 
total_proj$Glial <- total_proj$Astrocyte + total_proj$Microglia + total_proj$Oligodendrocyte
total_proj$Angiogenic <- total_proj$Stromal + total_proj$Endothelial 

cell.types <- c("Endothelial","Stromal","Astrocyte","Microglia","Mono","Oligodendrocyte","GABA","GLU","Tumor",
                "Neu","Bmem","Bnv","CD4mem","CD4nv","CD8mem","CD8nv","Treg","NK")

cell_types2 <- c(cell.types, "Angiogenic", "Glial", "Myeloid","Lymphoid","Tcell","Bcell","Immune")

for (cell in cell_types2){
  total_proj[,cell] <- total_proj[,cell] - mean(total_proj[,cell], na.rm=TRUE)
  total_proj[,cell] <- total_proj[,cell]/sd(total_proj[,cell], na.rm=TRUE)
}

rownames(total_proj) <- paste0(total_proj$Slide,"_",total_proj$Array)
both <- intersect(rownames(total_proj), colnames(agg_betas))
total_proj <- total_proj[both,]
agg_betas <- agg_betas[,both]

my_combos_to_test <- list(Unadj = NULL,
                          TumorOnly = "Tumor",
                          ManyBroad = c("Tumor","Angiogenic","Glial","Immune"),
                          ManyDetailedv1 = c("Tumor","Angiogenic","Astrocyte","Microglia" ,"Oligodendrocyte","Myeloid" ,"Lymphoid"),
                          ManyDetailedv2 = c("Tumor","Angiogenic","Astrocyte","Microglia" ,"Oligodendrocyte","Tcell" ,"Bcell" ,"NK" ,"Neu","Mono"))


# tmp_pheno <- total_proj[which(total_proj$Library %in% c("AST-LG", "AST-HG")),]
# tmp_betas <- agg_betas[, which(colnames(agg_betas) %in% rownames(tmp_pheno))]
# 
# DeconvoEWAS(df=tmp_pheno, # data frame with all data, colnames should match variable names within reg.model, cell.types, rownames should match colnames of betas 
#             betas=tmp_betas, # beta matrix for dataset
#             base_model="Library", # the right hand side of the regression model equation
#             user_p_adjust_method = "fdr", # adjustment method for multiple testing
#             cell_types=cell_types2, # list of cell types 
#             cell_combos_to_test = my_combos_to_test, # named list of vectors containing names of cell types to try, each entry of list will be one model, leave one entry as NULL to do unadj
#             model_type = "lm",
#             n_cpgs_to_test = nrow(betas), # how many CpGs to do EWAS with
#             cpg_filter_method = "deltaBeta", # method to filter top cpgs, either by total variance or by delta beta
#             cpg_user_names=NULL, # if cpg filter method is cpg list, should be names of cpgs user wants to test
#             group_var="Library", # what is the grouping variable
#             ref_group="AST-LG", # what is the reference group
#             case_group="AST-HG", # what is the group to compare 
#             save_output_csv = TRUE, # save the results in csv table
#             save_output_pdf = TRUE, # save a pdf version of the plot
#             file_output_prefix = "ASTHG_vs_ASTLG", # string to attach to saved file names
#             output_csv_dir = "/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/CapperEWAS", # path to folder to save files in
#             make_sum_equal_one = FALSE,
#             scale_cell_types = FALSE, # scale the cell types
#             confint.level=0.95)
# 
# 
# tmp_pheno <- total_proj[which(total_proj$Library %in% c("AST-LG", "OLG")),]
# tmp_betas <- agg_betas[, which(colnames(agg_betas) %in% rownames(tmp_pheno))]
# 
# DeconvoEWAS(df=tmp_pheno, # data frame with all data, colnames should match variable names within reg.model, cell.types, rownames should match colnames of betas 
#             betas=tmp_betas, # beta matrix for dataset
#             base_model="Library", # the right hand side of the regression model equation
#             user_p_adjust_method = "fdr", # adjustment method for multiple testing
#             cell_types=cell_types2, # list of cell types 
#             cell_combos_to_test = my_combos_to_test, # named list of vectors containing names of cell types to try, each entry of list will be one model, leave one entry as NULL to do unadj
#             model_type = "lm",
#             n_cpgs_to_test = nrow(betas), # how many CpGs to do EWAS with
#             cpg_filter_method = "deltaBeta", # method to filter top cpgs, either by total variance or by delta beta
#             cpg_user_names=NULL, # if cpg filter method is cpg list, should be names of cpgs user wants to test
#             group_var="Library", # what is the grouping variable
#             ref_group="AST-LG", # what is the reference group
#             case_group="OLG", # what is the group to compare 
#             save_output_csv = TRUE, # save the results in csv table
#             save_output_pdf = TRUE, # save a pdf version of the plot
#             file_output_prefix = "ASTHG_vs_OLG", # string to attach to saved file names
#             output_csv_dir = "/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/CapperEWAS", # path to folder to save files in
#             make_sum_equal_one = FALSE,
#             scale_cell_types = FALSE, # scale the cell types
#             confint.level=0.95)
# 
# 
# tmp_pheno <- total_proj[which(total_proj$Library %in% c("AST-LG", "GBM")),]
# tmp_betas <- agg_betas[, which(colnames(agg_betas) %in% rownames(tmp_pheno))]
# 
# DeconvoEWAS(df=tmp_pheno, # data frame with all data, colnames should match variable names within reg.model, cell.types, rownames should match colnames of betas 
#             betas=tmp_betas, # beta matrix for dataset
#             base_model="Library", # the right hand side of the regression model equation
#             user_p_adjust_method = "fdr", # adjustment method for multiple testing
#             cell_types=cell_types2, # list of cell types 
#             cell_combos_to_test = my_combos_to_test, # named list of vectors containing names of cell types to try, each entry of list will be one model, leave one entry as NULL to do unadj
#             model_type = "lm",
#             n_cpgs_to_test = nrow(betas), # how many CpGs to do EWAS with
#             cpg_filter_method = "deltaBeta", # method to filter top cpgs, either by total variance or by delta beta
#             cpg_user_names=NULL, # if cpg filter method is cpg list, should be names of cpgs user wants to test
#             group_var="Library", # what is the grouping variable
#             ref_group="AST-LG", # what is the reference group
#             case_group="GBM", # what is the group to compare 
#             save_output_csv = TRUE, # save the results in csv table
#             save_output_pdf = TRUE, # save a pdf version of the plot
#             file_output_prefix = "ASTHG_vs_GBM", # string to attach to saved file names
#             output_csv_dir = "/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/CapperEWAS", # path to folder to save files in
#             make_sum_equal_one = FALSE,
#             scale_cell_types = FALSE, # scale the cell types
#             confint.level=0.95)

tmp_pheno <- total_proj[which(total_proj$Library %in% c("GBM")),]
dec_vals <- quantile(tmp_pheno$Immune, probs = seq(0, 1, 1/10)) 

tmp_pheno$high_infil <- ifelse(tmp_pheno$Immune <= dec_vals[2], "low",
                                ifelse(tmp_pheno$Immune >= dec_vals[10], "high", "middle"
                                ))

tmp_pheno <- tmp_pheno[tmp_pheno$high_infil != "middle",]
tmp_betas <- agg_betas[, which(colnames(agg_betas) %in% rownames(tmp_pheno))]


DeconvoEWAS(df=tmp_pheno, # data frame with all data, colnames should match variable names within reg.model, cell.types, rownames should match colnames of betas 
            betas=tmp_betas, # beta matrix for dataset
            base_model="high_infil", # the right hand side of the regression model equation
            user_p_adjust_method = "fdr", # adjustment method for multiple testing
            cell_types=cell_types2, # list of cell types 
            cell_combos_to_test = my_combos_to_test, # named list of vectors containing names of cell types to try, each entry of list will be one model, leave one entry as NULL to do unadj
            model_type = "lm",
            n_cpgs_to_test = nrow(betas), # how many CpGs to do EWAS with
            cpg_filter_method = "deltaBeta", # method to filter top cpgs, either by total variance or by delta beta
            cpg_user_names=NULL, # if cpg filter method is cpg list, should be names of cpgs user wants to test
            group_var="high_infil", # what is the grouping variable
            ref_group="low", # what is the reference group
            case_group="high", # what is the group to compare 
            save_output_csv = TRUE, # save the results in csv table
            save_output_pdf = TRUE, # save a pdf version of the plot
            file_output_prefix = "GBM_hivsloinfil", # string to attach to saved file names
            output_csv_dir = "/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/CapperEWAS", # path to folder to save files in
            make_sum_equal_one = FALSE,
            scale_cell_types = FALSE, # scale the cell types
            confint.level=0.95)

tmp_pheno <- total_proj[which(total_proj$Library %in% c("AST-LG")),]
dec_vals <- quantile(tmp_pheno$Immune, probs = seq(0, 1, 1/10)) 

tmp_pheno$high_infil <- ifelse(tmp_pheno$Immune <= dec_vals[2], "low",
                               ifelse(tmp_pheno$Immune >= dec_vals[10], "high", "middle"
                               ))

tmp_pheno <- tmp_pheno[tmp_pheno$high_infil != "middle",]
tmp_betas <- agg_betas[, which(colnames(agg_betas) %in% rownames(tmp_pheno))]


DeconvoEWAS(df=tmp_pheno, # data frame with all data, colnames should match variable names within reg.model, cell.types, rownames should match colnames of betas 
            betas=tmp_betas, # beta matrix for dataset
            base_model="high_infil", # the right hand side of the regression model equation
            user_p_adjust_method = "fdr", # adjustment method for multiple testing
            cell_types=cell_types2, # list of cell types 
            cell_combos_to_test = my_combos_to_test, # named list of vectors containing names of cell types to try, each entry of list will be one model, leave one entry as NULL to do unadj
            model_type = "lm",
            n_cpgs_to_test = nrow(betas), # how many CpGs to do EWAS with
            cpg_filter_method = "deltaBeta", # method to filter top cpgs, either by total variance or by delta beta
            cpg_user_names=NULL, # if cpg filter method is cpg list, should be names of cpgs user wants to test
            group_var="high_infil", # what is the grouping variable
            ref_group="low", # what is the reference group
            case_group="high", # what is the group to compare 
            save_output_csv = TRUE, # save the results in csv table
            save_output_pdf = TRUE, # save a pdf version of the plot
            file_output_prefix = "ASTLG_hivsloinfil", # string to attach to saved file names
            output_csv_dir = "/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/CapperEWAS", # path to folder to save files in
            make_sum_equal_one = FALSE,
            scale_cell_types = FALSE, # scale the cell types
            confint.level=0.95)

tmp_pheno <- total_proj[which(total_proj$Library %in% c("OLG")),]
dec_vals <- quantile(tmp_pheno$Immune, probs = seq(0, 1, 1/10)) 

tmp_pheno$high_infil <- ifelse(tmp_pheno$Immune <= dec_vals[2], "low",
                               ifelse(tmp_pheno$Immune >= dec_vals[10], "high", "middle"
                               ))

tmp_pheno <- tmp_pheno[tmp_pheno$high_infil != "middle",]
tmp_betas <- agg_betas[, which(colnames(agg_betas) %in% rownames(tmp_pheno))]


DeconvoEWAS(df=tmp_pheno, # data frame with all data, colnames should match variable names within reg.model, cell.types, rownames should match colnames of betas 
            betas=tmp_betas, # beta matrix for dataset
            base_model="high_infil", # the right hand side of the regression model equation
            user_p_adjust_method = "fdr", # adjustment method for multiple testing
            cell_types=cell_types2, # list of cell types 
            cell_combos_to_test = my_combos_to_test, # named list of vectors containing names of cell types to try, each entry of list will be one model, leave one entry as NULL to do unadj
            model_type = "lm",
            n_cpgs_to_test = nrow(betas), # how many CpGs to do EWAS with
            cpg_filter_method = "deltaBeta", # method to filter top cpgs, either by total variance or by delta beta
            cpg_user_names=NULL, # if cpg filter method is cpg list, should be names of cpgs user wants to test
            group_var="high_infil", # what is the grouping variable
            ref_group="low", # what is the reference group
            case_group="high", # what is the group to compare 
            save_output_csv = TRUE, # save the results in csv table
            save_output_pdf = TRUE, # save a pdf version of the plot
            file_output_prefix = "OLG_hivsloinfil", # string to attach to saved file names
            output_csv_dir = "/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/CapperEWAS", # path to folder to save files in
            make_sum_equal_one = FALSE,
            scale_cell_types = FALSE, # scale the cell types
            confint.level=0.95)

tmp_pheno <- total_proj[which(total_proj$Library %in% c("AST-HG")),]
dec_vals <- quantile(tmp_pheno$Immune, probs = seq(0, 1, 1/10)) 

tmp_pheno$high_infil <- ifelse(tmp_pheno$Immune <= dec_vals[2], "low",
                                ifelse(tmp_pheno$Immune >= dec_vals[10], "high", "middle"
                                ))

tmp_pheno <- tmp_pheno[tmp_pheno$high_infil != "middle",]
tmp_betas <- agg_betas[, which(colnames(agg_betas) %in% rownames(tmp_pheno))]


DeconvoEWAS(df=tmp_pheno, # data frame with all data, colnames should match variable names within reg.model, cell.types, rownames should match colnames of betas 
            betas=tmp_betas, # beta matrix for dataset
            base_model="high_infil", # the right hand side of the regression model equation
            user_p_adjust_method = "fdr", # adjustment method for multiple testing
            cell_types=cell_types2, # list of cell types 
            cell_combos_to_test = my_combos_to_test, # named list of vectors containing names of cell types to try, each entry of list will be one model, leave one entry as NULL to do unadj
            model_type = "lm",
            n_cpgs_to_test = nrow(betas), # how many CpGs to do EWAS with
            cpg_filter_method = "deltaBeta", # method to filter top cpgs, either by total variance or by delta beta
            cpg_user_names=NULL, # if cpg filter method is cpg list, should be names of cpgs user wants to test
            group_var="high_infil", # what is the grouping variable
            ref_group="low", # what is the reference group
            case_group="high", # what is the group to compare 
            save_output_csv = TRUE, # save the results in csv table
            save_output_pdf = TRUE, # save a pdf version of the plot
            file_output_prefix = "ASTHG_hivsloinfil", # string to attach to saved file names
            output_csv_dir = "/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Processed_data/Oct2023_Analysis/CapperEWAS", # path to folder to save files in
            make_sum_equal_one = FALSE,
            scale_cell_types = FALSE, # scale the cell types
            confint.level=0.95)