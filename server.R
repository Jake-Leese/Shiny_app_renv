
server <- function(input, output, session){
  
  ####################################################################
  # Generate ArchR UMAP plots automatically when inputs change
  observeEvent(input$UMAP_subset_ArchR, {
    output$dimplot_ArchR <- renderPlot({
      ArchR_dimplot(ArchR_list[[input$UMAP_subset_ArchR]], name = input$dimplot_ArchR_groupby) +
        my_theme
    }, height = function() { session$clientData$output_dimplot_ArchR_width * 0.8 })
  })
  
  observeEvent(input$featureplot_ArchR_TF, {
    output$featureplot_ArchR <- renderPlot({
      ArchR_featureplot(ArchR_list[[input$UMAP_subset_ArchR]], TF = input$featureplot_ArchR_TF, datatype = input$featureplot_ArchR_datatype) +
        my_theme
    }, height = function() { session$clientData$output_dimplot_ArchR_width * 0.8 })
  })
  
  ####################################################################
  # Generate ArchR Genome Browser 
  observeEvent(input$run_gBrowser, {
    output$gbrowser_ArchR <- renderPlot(
      gBrowser_plot(ArchR = ArchR_list[[input$gbrowser_subset_ArchR]], 
                    group_by = input$gbrowser_ArchR_groupby, 
                    region = input$gbrowser_ArchR_region,
                    extend_by = input$extend_by)
    )
  }
  )
    
}

