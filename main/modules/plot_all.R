# ===================================================================================================
# ALL EXPERIMENTS PCA PLOT MODULE
# ===================================================================================================
# This module contains the plot generation function for the "All Experiments" view
# Shows overview of all samples colored by FPE number and shaped by TNFa treatment

generate_all_plot <- function(cts, coldata) {
  tryCatch({
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
    
    # Create the plot
    ggplot(pcaData, aes(x = PC1, y = PC2, color = FPE.num, shape = TNFa.positive)) +
      geom_point(size = 3) +
      xlab(paste0('PC1: ', percentVar[1], '% variance')) +
      ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
      coord_fixed() +
      scale_color_brewer(palette = 'Set1') +
      ggtitle(label = 'All Experiments - PCA Plot') +
      theme_minimal() +
      theme(
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      )
    
  }, error = function(e) {
    # Return error plot if generation fails
    ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
      geom_text(label = paste("Error in All Experiments plot:", e$message), 
                size = 4, color = "red") +
      theme_void() +
      labs(title = "Plot Generation Error")
  })
}
