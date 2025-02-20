## Function to extract cell types avaliable for a given stage
check_cell_types <- function(stage){
  metadata <- SEACells_metadata
  if (!stage == "Full Data"){metadata <- metadata %>% filter(stage == !!stage)}
  cell_types <- unique(metadata$scHelper_cell_type)
  return(cell_types)
}

## Function to plot genome browser plot
gBrowser_plot <- function(ArchR, group_by, region, extend_by){
  # set colours:
  if (grepl("scHelper_cell_type", group_by) ){
    available_cell_types <- as.vector(unique(getCellColData(ArchR, select = group_by)[,1]))
    available_cell_types <- replace(available_cell_types, is.na(available_cell_types), "Unmapped")
    pal <- scHelper_cell_type_colours[available_cell_types]
  } 
  if (group_by == "stage"){
    pal = stage_colours
  }
  if (group_by == "clusters"){
    pal = NULL
  }
  
  # REGION input:
  if (grepl("chr", region)){
    print("region detected!")
    df <- region_string_to_df(region)
    GRanges_highlight <- makeGRangesFromDataFrame(df)
    df_extended <- df %>% 
      mutate(start = start - extend_by) %>%
      mutate(end = end + extend_by)
    GRanges_plot <- makeGRangesFromDataFrame(df_extended)
    
    p <- plotBrowserTrack(
      ArchRProj = ArchR, 
      groupBy = group_by, 
      region = GRanges_plot,
      highlight = GRanges_highlight,
      features = getPeakSet(ArchR),
      plotSummary = c("bulkTrack", "geneTrack", "featureTrack"),
      facetbaseSize = 12,
      baseSize = 12,
      title = "",
      pal = pal
    )
    grid::grid.draw(p)
  } else {
    # GENE SYMBOL input:
    print("gene symbol detected!")
    p <- plotBrowserTrack(
      ArchRProj = ArchR, 
      groupBy = group_by, 
      geneSymbol = region,
      upstream = extend_by,
      downstream = extend_by,
      features = getPeakSet(ArchR),
      plotSummary = c("bulkTrack", "geneTrack", "featureTrack"),
      facetbaseSize = 12,
      baseSize = 12,
      title = "",
      pal = pal
    )
    grid::grid.draw(p[[1]])
    }
  
}

##### function to take a string (eg "chr9:17029120-17029653") and turn into a dataframe that can be converted to granges
region_string_to_df <- function(region){
  coords <- strsplit(region, ":")[[1]][2]
  df <- data.frame(c(strsplit(region, ":")[[1]][1]), 
                   c(strsplit(coords, "-")[[1]][1]),
                   strsplit(coords, "-")[[1]][2])
  colnames(df) <- c("chromosome", "start", "end")
  df$start <- as.numeric(df$start)
  df$end <- as.numeric(df$end)
  return(df)
}

### function to make ArchR DimPlot with correct colours depending on 'name' which is the chosen 'dimplot_ArchR_groupby'
ArchR_dimplot <- function(ArchR, name){
  
  # if plotting stage
  if(name == "stage"){
    print(plotEmbedding(ArchR, name = name,
                        plotAs = "points", size = 1.8, randomize = TRUE,
                        labelSize = 8, labelAsFactors = FALSE, baseSize = 20, legendSize = 20,
                        pal = stage_colours) + 
          theme_ArchR(legendTextSize = 17, baseSize = 17, plotMarginCm = 0.5) + 
      guides(colour = guide_legend(override.aes = list(size=7))))
  }
  
  # if plotting cell states
  if (grepl("scHelper_cell_type", name)){
    available_cell_types <- as.vector(unique(getCellColData(ArchR, select = name)[,1]))
    available_cell_types <- replace(available_cell_types, is.na(available_cell_types), "Unmapped")
    scHelper_cols <- scHelper_cell_type_colours[available_cell_types]
    print(plotEmbedding(ArchR, name = name,
                  plotAs = "points", size = 1.8, randomize = TRUE,
                  labelSize = 0, labelAsFactors = FALSE, baseSize = 0, legendSize = 20,
                  pal = scHelper_cols) + 
            theme_ArchR(legendTextSize = 17, baseSize = 17, plotMarginCm = 0.5) + 
      guides(colour = guide_legend(override.aes = list(size=7))))
  }

    # if plotting clusters
  if (name == "clusters"){
    print(plotEmbedding(ArchR, name = name,
                  plotAs = "points", size = 1.8, randomize = TRUE,
                  labelSize = 6, labelAsFactors = FALSE, baseSize = 20, legendSize = 20) + 
      theme_ArchR(legendTextSize = 17, baseSize = 17, plotMarginCm = 0.5) + 
      guides(colour = guide_legend(override.aes = list(size=7))))
  }
  
}

## function to make ArchR feature plot with input of TF and what kind of data to plot of it: gene score, gex or activity from ChromVar
ArchR_featureplot <- function(ArchR, TF, datatype){
    if(datatype == "MotifMatrix"){
      cols = "solarExtra"
      TF = paste0("z:", TF)}
    if(datatype == "GeneIntegrationMatrix"){cols = "blueYellow"}
  if(datatype == "GeneScoreMatrix"){cols = "horizon"}
    print(plotEmbedding(ArchR, colorBy = datatype, name = TF, continuousSet = cols,
                        plotAs = "points", size = 1.8, randomize = TRUE,
                        labelSize = 10, labelAsFactors = FALSE, baseSize = 20, legendSize = 20) + 
            theme_ArchR(legendTextSize = 17, baseSize = 17, plotMarginCm = 0.5)) + theme(legend.key.width=unit(2,"cm"))
  }

