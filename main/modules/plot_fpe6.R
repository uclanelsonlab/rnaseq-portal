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
    
    # Use precomputed PCA if available
    pca_file <- '../data/fpe6_pca.rds'
    if (file.exists(pca_file)) {
      pre <- readRDS(pca_file)
      pcaData <- pre$pcaData
      percentVar <- round(100 * pre$percentVar)
    } else {
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
        intgroup = c('FPE.num', 'participant_id', 'treatment', 'treatment.time', 'replicate.num', 
                     'sub.treatment', 'TNFa.positive', 'none'), 
        returnData = TRUE)
      percentVar <- round(100 * attr(pcaData, 'percentVar'))
    }
    
    # percentVar already set
    
    # Reorder sub.treatment factor levels
    pcaData$sub.treatment <- factor(pcaData$sub.treatment, 
                                     levels = c('TNFa', 'rapa', 'TGFb+rapa', 'TGFb', 'TGFb+SB', 
                                               'SB', 'rapa+SB', 'KC7F2', 'TGFb+KC7F2', 
                                               'EGFRi', 'EGFRi+EGF', 'EGF'))
    
    # Create the plot with special highlighting for 'none' samples
    ggplot(data = pcaData, aes(x = PC1, y = PC2, color = sub.treatment, shape = TNFa.positive)) +
      geom_point(size = 3) +
      geom_point(data = pcaData[pcaData$none, ], size = 3, color = 'black') +
      xlab(paste0('PC1: ', percentVar[1], '% variance')) +
      ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
      coord_fixed() +
      scale_shape_discrete(labels = c("NO", "YES")) +
      ggtitle(label = 'Fibroblast Priming 6') +
      theme_gray(base_size = 20)+
      guides(
        shape = guide_legend(title = "TNF?", order = 1),
        color = guide_legend(title = "Sub-treatment", order = 2)
      )
    
  }, error = function(e) {
    # Return error plot if generation fails
    ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
      geom_text(label = paste("Error in FPE6 plot:", e$message), 
                size = 4, color = "red") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "gray95", color = NA)) +
      labs(title = "Plot Generation Error")
  })
}
