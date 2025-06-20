library(shiny)
library(bs4Dash)
library(data.table)
library(tidyverse)
library(viridis)
library(mgcv)
library(patchwork)
library(ComplexHeatmap)
library(scHelper)
library(ArchR)

source('./custom_functions.R')

# Don't know what this does
options(scipen = 1)
options(digits = 2)

# prevent any attempt to load graphics within global.R set up
pdf(NULL)

####################################################################
##################       Aesthetic params      #####################

scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                              'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'node', 'streak', 
                              'PGC', 'BI', 'meso', 'endo', 'MIXED', 'Unmapped',
                              'Neural', 'Placodal', 'Non-neural', 'Contam')
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", 
                                "#53A651", "#6D8470", "#87638F", "#A5548D", "#C96555", "#ED761C", 
                                "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30", "#CC9F2C", "#AD6428", 
                                "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB",
                                "#A5718D", "#3F918C", "#ed5e5f", "#9792A3",
                                "#7C8483", "#EAEAEA")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 
                                       'MB','vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo',
                                       'Neural', 'Placodal', 'Non-neural', 'Contam',
                                       'MIXED', 'Unmapped')
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
names(stage_colours) <- stage_order

my_theme <- theme(axis.text=element_text(size=14),
                  axis.title=element_text(size=16))

########################################################################
####################       ArchR objects      ##########################

# load ArchR objects
HH5_ArchR <- loadArchRProject(path = "/home/ubuntu/RDS/sc_atac_data_EH/processed_data/Shiny_ArchR_data/scATAC_processed_files_ArchR/HH5_Save_ArchR_optimised", showLogo = FALSE)
HH6_ArchR <- loadArchRProject(path = "/home/ubuntu/RDS/sc_atac_data_EH/processed_data/Shiny_ArchR_data/scATAC_processed_files_ArchR/HH6_Save_ArchR_optimised", showLogo = FALSE)
HH7_ArchR <- loadArchRProject(path = "/home/ubuntu/RDS/sc_atac_data_EH/processed_data/Shiny_ArchR_data/scATAC_processed_files_ArchR/HH7_Save_ArchR_optimised", showLogo = FALSE)
ss4_ArchR <- loadArchRProject(path = "/home/ubuntu/RDS/sc_atac_data_EH/processed_data/Shiny_ArchR_data/scATAC_processed_files_ArchR/ss4_Save_ArchR_optimised", showLogo = FALSE)
ss8_ArchR <- loadArchRProject(path = "/home/ubuntu/RDS/sc_atac_data_EH/processed_data/Shiny_ArchR_data/scATAC_processed_files_ArchR/ss8_Save_ArchR_optimised", showLogo = FALSE)
FullData_ArchR <- loadArchRProject(path = "/home/ubuntu/RDS/sc_atac_data_EH/processed_data/Shiny_ArchR_data/scATAC_processed_files_ArchR/FullData_Save_ArchR_optimised", showLogo = FALSE)

# make the schelper cell type metadata of the fulldata named like the stages
FullData_ArchR <- addCellColData(FullData_ArchR, data = FullData_ArchR$scHelper_cell_type, cells = rownames(getCellColData(FullData_ArchR)), name = "transferred_scHelper_cell_type")
FullData_ArchR <- addCellColData(FullData_ArchR, data = FullData_ArchR$scHelper_cell_type_broad, cells = rownames(getCellColData(FullData_ArchR)), name = "transferred_scHelper_cell_type_broad")

# use MAGIC to impute weights
HH5_ArchR <- addImputeWeights(HH5_ArchR)
HH6_ArchR <- addImputeWeights(HH6_ArchR)
HH7_ArchR <- addImputeWeights(HH7_ArchR)
ss4_ArchR <- addImputeWeights(ss4_ArchR)
ss8_ArchR <- addImputeWeights(ss8_ArchR)
FullData_ArchR <- addImputeWeights(FullData_ArchR)

# make list of objects
ArchR_list <- list(FullData_ArchR, HH5_ArchR, HH6_ArchR, HH7_ArchR, ss4_ArchR, ss8_ArchR)
names(ArchR_list) <- c("Full Data", "HH5", "HH6", "HH7", "ss4", "ss8")

#####################################################################
####################       UI options      ##########################

# Data subsets including full data
data_subsets = c("Full Data", "HH5", "HH6", "HH7", "ss4", "ss8")

# Data subsets excluding full data
stages = c("HH5", "HH6", "HH7", "ss4", "ss8")

# Potential ways to group ArchR cells
ArchR_groupby_options <- c("Stage" = "stage", "Clusters" = "clusters", "Cell state" = "transferred_scHelper_cell_type", "Broad cell state" = "transferred_scHelper_cell_type_broad")

# Potential data modalities for a TF
TF_datatype_options <- c("Gene Score" = "GeneScoreMatrix", "Gene Expression (from scRNA-seq data)" = "GeneIntegrationMatrix", "TF activity" = "MotifMatrix")

# All potential TFs for which to plot featureplots. Should be present in all 3 datatypes. 
# Have checked that this final list of 384 TFs is the same across all stages, here have used ss8 but should be same for all
motifs <- unique(gsub("^(z:|deviations:)", "", getFeatures(ss8_ArchR, useMatrix = "MotifMatrix")))
genes <- getFeatures(ss8_ArchR, useMatrix = "GeneScoreMatrix")
matches_gene_score <- sapply(motifs, function(x) any(grepl(x, genes)))
TF_options_1 <- sort(motifs[matches_gene_score]) # 668
int_genes <- getFeatures(ss8_ArchR, useMatrix = "GeneIntegrationMatrix")
matches <- sapply(TF_options_1, function(x) any(grepl(x, int_genes)))
TF_options <- sort(TF_options_1[matches]) # 346

# Pre-loading UMAPs
# Initialize nested list
UMAP_Plots <- list()

# Loop through ArchR objects and generate plots
for (sample_name in names(ArchR_list)) {
  plots <- list()
  for (group in ArchR_groupby_options) {
    plots[[group]] <- ArchR_dimplot(ArchR_list[[sample_name]], group)
  }
  UMAP_Plots[[sample_name]] <- plots
}
