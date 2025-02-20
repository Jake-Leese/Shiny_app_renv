library(shinythemes)

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




