# ===================================================================================================
# EXPERIMENT 4 (FPE4) PCA PLOT MODULE
# ===================================================================================================
# This module contains the plot generation function for Fibroblast Priming Experiment #4
# Shows samples colored by source and shaped by treatment

generate_fpe4_plot <- function(cts, coldata) {
  tryCatch({
    # Filter samples for FPE4 only
    selected.samples <- coldata$FPE.num == 'FPE4'
    
    # Check if we have FPE4 samples
    if(sum(selected.samples) == 0) {
      stop("No FPE4 samples found in the data")
    }
    
    # Use precomputed PCA if available
    pca_file <- '../data/fpe4_pca.rds'
    if (file.exists(pca_file)) {
      pre <- readRDS(pca_file)
      pcaData <- pre$pcaData
      percentVar <- round(100 * pre$percentVar)
    } else {
      # Create DESeq2 dataset with FPE4 samples only
      dds <- DESeqDataSetFromMatrix(
        countData = cts[, selected.samples], 
        colData = coldata[selected.samples, ], 
        design = ~ treatment)
      
      # Apply variance stabilizing transformation
      vsd <- vst(dds, blind = FALSE)
      
      # Generate PCA data
      pcaData <- plotPCA(
        object = vsd, 
        intgroup = c('FPE.num', 'participant_id', 'treatment', 'treatment.time', 'replicate.num'), 
        returnData = TRUE)
      percentVar <- round(100 * attr(pcaData, 'percentVar'))
    }
    
    
    # Create the plot
    ggplot(data = pcaData, aes(x = PC1, y = PC2, shape = treatment, color = source)) +
      geom_point(size = 3) +
      xlab(paste0('PC1: ', percentVar[1], '% variance')) +
      ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
      coord_fixed() +
      scale_color_brewer(palette = 'Dark2') +
      scale_shape_discrete(labels = c("None", "TNFa")) +
      ggtitle(label = 'Fibroblast Priming Experiment 4') +
      theme_gray(base_size = 20) +
      guides(
        color = guide_legend(title = "Source"),
        shape = guide_legend(title = "Treatment")
      )
    
  }, error = function(e) {
    # Return error plot if generation fails
    ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
      geom_text(label = paste("Error in FPE4 plot:", e$message), 
                size = 4, color = "red") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "gray95", color = NA)) +
      labs(title = "Plot Generation Error")
  })
}
