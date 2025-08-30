# ===================================================================================================
# EXPERIMENT 6 (FPE6) PCA PLOT MODULE
# ===================================================================================================
# This module contains the plot generation function for Fibroblast Priming Experiment #6
# Shows samples colored by sub-treatment and shaped by TNFa status, with special highlighting

generate_fpe6_plot <- function(cts, coldata, coldata_fpe6) {
  tryCatch({
    # Filter samples for FPE6 only
    selected.samples <- coldata$FPE.num == 'FPE6'
    
    # Check if we have FPE6 samples
    if(sum(selected.samples) == 0) {
      stop("No FPE6 samples found in the data")
    }
    
    # Use the special FPE6 processed coldata
    coldata_subset <- coldata_fpe6[selected.samples, ]
    
    # Create DESeq2 dataset with FPE6 samples only
    dds <- DESeqDataSetFromMatrix(
      countData = cts[, selected.samples], 
      colData = coldata_subset, 
      design = ~ treatment)
    
    # Apply variance stabilizing transformation
    vsd <- vst(dds, blind = FALSE)
    
    # Generate PCA data
    pcaData <- plotPCA(
      object = vsd, 
      intgroup = c('FPE.num', 'participant_id', 'treatment', 'treatment.time', 'replicate.num'), 
      returnData = TRUE)
    
    # Calculate percentage variance explained
    percentVar <- round(100 * attr(pcaData, 'percentVar'))
    
    # Create the plot with special highlighting for 'none' samples
    ggplot(data = pcaData, aes(x = PC1, y = PC2, color = sub.treatment, shape = TNFa.positive)) +
      geom_point(size = 3) +
      geom_point(data = pcaData[pcaData$none, ], size = 3, color = 'black') +
      xlab(paste0('PC1: ', percentVar[1], '% variance')) +
      ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
      coord_fixed() +
      scale_color_brewer(palette = 'Set3') +
      ggtitle(label = 'Fibroblast Priming #6') +
      theme_minimal() +
      theme(
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      guides(
        color = guide_legend(title = "Sub-treatment"),
        shape = guide_legend(title = "TNFa Positive")
      )
    
  }, error = function(e) {
    # Return error plot if generation fails
    ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
      geom_text(label = paste("Error in FPE6 plot:", e$message), 
                size = 4, color = "red") +
      theme_void() +
      labs(title = "Plot Generation Error")
  })
}
