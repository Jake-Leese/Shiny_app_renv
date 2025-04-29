library(shiny)
library(shinythemes)
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

options(scipen = 1)
options(digits = 2)

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

# Paths to ArchR objects
archr_paths <- list(
  "Full Data" = "/home/ubuntu/RDS/sc_atac_data_EH/processed_data/FullData/FullData_Save-ArchR/FullData_Save-ArchR",
  "HH5" = "/home/ubuntu/RDS/sc_atac_data_EH/processed_data/HH5/HH5_Save-ArchR/HH5_Save-ArchR",
  "HH6" = "/home/ubuntu/RDS/sc_atac_data_EH/processed_data/HH6/HH6_Save-ArchR/HH6_Save-ArchR",
  "HH7" = "/home/ubuntu/RDS/sc_atac_data_EH/processed_data/HH7/HH7_Save-ArchR/HH7_Save-ArchR",
  "ss4" = "/home/ubuntu/RDS/sc_atac_data_EH/processed_data/ss4/ss4_Save-ArchR/ss4_Save-ArchR",
  "ss8" = "/home/ubuntu/RDS/sc_atac_data_EH/processed_data/ss8/ss8_Save-ArchR/ss8_Save-ArchR"
)

# Function to load ArchR object
dynamicLoadArchR <- function(dataset) {
  if (!is.null(dataset) && dataset %in% names(archr_paths)) {
    ArchRProj <- loadArchRProject(path = archr_paths[[dataset]], showLogo = FALSE)
    if (dataset == "Full Data") {
      ArchRProj <- addCellColData(ArchRProj, data = ArchRProj$scHelper_cell_type, cells = rownames(getCellColData(ArchRProj)), name = "transferred_scHelper_cell_type")
      ArchRProj <- addCellColData(ArchRProj, data = ArchRProj$scHelper_cell_type_broad, cells = rownames(getCellColData(ArchRProj)), name = "transferred_scHelper_cell_type_broad")
      ArchRProj <- addImputeWeights(ArchRProj)
      return(ArchRProj)
    } else {
      ArchRProj <- addImputeWeights(ArchRProj)
      return(ArchRProj)
    }
  }
  return(NULL)
}


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
ss8_ArchR <- loadArchRProject(path = "/home/ubuntu/RDS/sc_atac_data_EH/processed_data/ss8/ss8_Save-ArchR/ss8_Save-ArchR", showLogo = FALSE)
motifs <- unique(gsub("^(z:|deviations:)", "", getFeatures(ss8_ArchR, useMatrix = "MotifMatrix")))
genes <- getFeatures(ss8_ArchR, useMatrix = "GeneScoreMatrix")
matches_gene_score <- sapply(motifs, function(x) any(grepl(x, genes)))
TF_options_1 <- sort(motifs[matches_gene_score]) # 668
int_genes <- getFeatures(ss8_ArchR, useMatrix = "GeneIntegrationMatrix")
matches <- sapply(TF_options_1, function(x) any(grepl(x, int_genes)))
TF_options <- sort(TF_options_1[matches]) # 346
rm(ss8_ArchR)

#######################################################################################
################################    Server    #########################################

server <- function(input, output, session) {
  
  # Reactive value to store the loaded ArchR object
  archr_object <- reactiveVal(NULL)
  
  ####################################################################
  # Generate ArchR UMAP plots automatically when inputs change
  
  # Function to update the loaded ArchR object
  observeEvent(input$UMAP_subset_ArchR, {
    # Remove previous object from memory
    archr_object(NULL)
    gc()
    
    # Load new object
    archr_object(dynamicLoadArchR(input$UMAP_subset_ArchR))
  })
  
  observeEvent(input$UMAP_subset_ArchR, {
    output$dimplot_ArchR <- renderPlot({
      ArchR_dimplot(archr_object(), name = input$dimplot_ArchR_groupby) +
        theme_minimal()
    }, height = function() { session$clientData$output_dimplot_ArchR_width * 0.8 })
  })
  
  observeEvent(input$featureplot_ArchR_TF, {
    output$featureplot_ArchR <- renderPlot({
      ArchR_featureplot(archr_object(), TF = input$featureplot_ArchR_TF, datatype = input$featureplot_ArchR_datatype) +
        theme_minimal()
    }, height = function() { session$clientData$output_dimplot_ArchR_width * 0.8 })
  })
  
  ####################################################################
  # Generate ArchR Genome Browser 
  observeEvent(input$gbrowser_subset_ArchR, {
    # Remove previous object from memory
    archr_object(NULL)
    gc()
    
    # Load new object
    archr_object(dynamicLoadArchR(input$gbrowser_subset_ArchR))
  })
  
  observeEvent(input$run_gBrowser, {
    output$gbrowser_ArchR <- renderPlot(
      gBrowser_plot(ArchR = archr_object(), 
                    group_by = input$gbrowser_ArchR_groupby, 
                    region = input$gbrowser_ArchR_region,
                    extend_by = input$extend_by)
    )
  })
}

#######################################################################################
#####################     ArchR UMAPs and Feature plots    ############################

tab_ArchR_umap <- tabItem(tabName = "ArchR_UMAP",
                          
                          fluidRow(
                            # Left column for dimplot_ArchR
                            column(6,
                                   fluidRow(
                                     column(12,
                                            radioButtons("UMAP_subset_ArchR", "Data subset to plot", data_subsets, inline = TRUE, width = '800')
                                     )
                                   ),
                                   fluidRow(
                                     column(12,
                                            radioButtons("dimplot_ArchR_groupby", "UMAP coloured by", ArchR_groupby_options, inline = TRUE, width = '800')
                                     )
                                   ),
                                   fluidRow(
                                     column(12,
                                            box(plotOutput("dimplot_ArchR"), width = 12, height = "35vw")  # Adjusted height to maintain aspect ratio
                                     )
                                   )
                            ),
                            # Right column for featureplot_ArchR
                            column(6,
                                   fluidRow(
                                     column(12,
                                            radioButtons("featureplot_ArchR_datatype", "Data to plot", TF_datatype_options, inline = TRUE, width = '800')
                                     )
                                   ),
                                   fluidRow(
                                     column(12,
                                            selectInput("featureplot_ArchR_TF", "Transcription factor to plot", choices = TF_options, width = "250")
                                     )
                                   ),
                                   fluidRow(
                                     column(12,
                                            box(plotOutput("featureplot_ArchR"), width = 12, height = "35vw")  # Adjusted height to maintain aspect ratio
                                     )
                                   )
                            )
                          )
)


##############################################################################
#####################     ArchR Genome Browser    ############################

tab_ArchR_gbrowser <- tabItem(tabName = "ArchR_Genome_Browser",
                              fluidRow(
                                column(12,
                                       radioButtons("gbrowser_subset_ArchR", "Data subset to plot", data_subsets, inline = TRUE, selected = 'ss8', width = '800')
                                )
                              ),
                              fluidRow(
                                column(12,
                                       radioButtons("gbrowser_ArchR_groupby", "How to group cells", ArchR_groupby_options, inline = TRUE, selected = 'clusters', width = '800')
                                )
                              ),
                              fluidRow(
                                #column(12, selectInput("gbrowser_ArchR_region", "Select region or gene to visualise", selected = 'SOX8', choices = NULL, multiple = FALSE, width = "250"))
                                column(12, textInput("gbrowser_ArchR_region", "Select region or gene to visualise", width = "250"))
                              ),
                              fluidRow(
                                column(12, sliderInput("extend_by", "Size around ROI in base pairs:", min = 0, max = 50000, value = 5000))
                              ),
                              fluidRow(column(6, actionButton("run_gBrowser", "Generate plot!"))),
                              fluidRow(
                                column(12,
                                       box(
                                         plotOutput("gbrowser_ArchR"),
                                         width = 12,
                                         height = "35vw"
                                       )
                                )
                              )
)

ui <- dashboardPage(
  header = dashboardHeader(
    title = dashboardBrand(
      title = "10x ATAC Neural Plate Border",
      href = "https://github.com/evaham1/atac_neural_plate_border"
    )
  ),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("ArchR_UMAP", tabName = "ArchR_UMAP", icon = icon("arrows-alt")),
      menuItem("ArchR_Genome_Browser", tabName = "ArchR_Genome_Browser", icon = icon("arrows-alt"))
    )
  ),
  
  dashboardBody(
    tabItems(
      tab_ArchR_umap,
      tab_ArchR_gbrowser
    )
  )
)

shinyApp(ui, server)

