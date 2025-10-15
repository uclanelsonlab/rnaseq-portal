# ===================================================================================================
# ALL EXPERIMENTS PCA PLOT MODULE
# ===================================================================================================
# This module contains the plot generation function for the "All Experiments" view
# Shows overview of all samples colored by FPE number and shaped by TNFa treatment

generate_all_plot <- function(cts, coldata) {
  tryCatch({
    # Try to use precomputed PCA if available
    pca_file <- '../data/fpe_all_pca.rds'
    if (file.exists(pca_file)) {
      pre <- readRDS(pca_file)
      pcaData <- pre$pcaData
      percentVar <- round(100 * pre$percentVar)
    } else {
      # Create DESeq2 dataset with all samples
      dds <- DESeqDataSetFromMatrix(
        countData = cts, 
        colData = coldata, 
        design = ~ treatment + treatment.time)
      
      # Apply variance stabilizing transformation
      vsd <- vst(dds, blind = FALSE)
      
      # Generate PCA data
      pcaData <- plotPCA(
        object = vsd, 
        intgroup = c('FPE.num', 'participant_id', 'treatment', 'TNFa.positive', 'treatment.time', 'replicate.num'), 
        returnData = TRUE)
      
      # Calculate percentage variance explained
      percentVar <- round(100 * attr(pcaData, 'percentVar'))
    }
    
    # Create the plot
    ggplot(pcaData, aes(x = PC1, y = PC2, color = FPE.num, shape = TNFa.positive)) +
      geom_point(size = 3) +
      xlab(paste0('PC1: ', percentVar[1], '% variance')) +
      ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
      coord_fixed() +
      scale_color_brewer(palette = 'Set1', name = "Experiment") +
      scale_shape_discrete(name = "TNF?", labels = c("NO", "YES")) +
      ggtitle(label = 'All Experiments') +
      theme_gray(base_size = 20) 
    
  }, error = function(e) {
    # Return error plot if generation fails
    ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
      geom_text(label = paste("Error in All Experiments plot:", e$message), 
                size = 4, color = "red") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "gray95", color = NA)) +
      labs(title = "Plot Generation Error")
  })
}
