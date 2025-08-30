# ===================================================================================================
# EXPERIMENT 7 (FPE7) PCA PLOT MODULE
# ===================================================================================================
# This module contains the plot generation function for Fibroblast Priming Experiment #7
# Shows samples colored by participant ID and shaped by treatment, with affected samples highlighted

generate_fpe7_plot <- function(cts, coldata) {
  tryCatch({
    # Filter samples for FPE7 only
    selected.samples <- coldata$FPE.num == 'FPE7'
    
    # Check if we have FPE7 samples
    if(sum(selected.samples) == 0) {
      stop("No FPE7 samples found in the data")
    }
    
    # Create DESeq2 dataset with FPE7 samples only
    dds <- DESeqDataSetFromMatrix(
      countData = cts[, selected.samples], 
      colData = coldata[selected.samples, ], 
      design = ~ treatment)
    
    # Apply variance stabilizing transformation
    vsd <- vst(dds, blind = FALSE)
    
    # Generate PCA data
    pcaData <- plotPCA(
      object = vsd, 
      intgroup = c('FPE.num', 'participant_id', 'affected_status', 'treatment', 'TNFa.positive', 'treatment.time', 'replicate.num'), 
      returnData = TRUE)
    
    # Calculate percentage variance explained
    percentVar <- round(100 * attr(pcaData, 'percentVar'))
    
    # Create the plot with affected samples highlighted in black
    ggplot(data = pcaData, aes(x = PC1, y = PC2, fill = participant_id, shape = treatment, color = participant_id)) +
      geom_point(size = 3) +
      geom_point(data = pcaData[pcaData$affected_status == 'Affected', ], 
                 size = 3, 
                 color = 'black', 
                 show.legend = FALSE) +
      xlab(paste0('PC1: ', percentVar[1], '% variance')) +
      ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
      coord_fixed() +
      scale_fill_brewer(palette = 'Dark2') +
      scale_color_brewer(palette = 'Dark2') +
      scale_shape_manual(values = 21:24) + 
      ggtitle(label = 'Fibroblast Priming #7') +
      theme_minimal() +
      theme(
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      guides(
        fill = guide_legend(title = "Participant ID", override.aes = list(shape = 21)),
        color = guide_legend(title = "Participant ID", override.aes = list(shape = 21)),
        shape = guide_legend(title = "Treatment")
      )
    
  }, error = function(e) {
    # Return error plot if generation fails
    ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
      geom_text(label = paste("Error in FPE7 plot:", e$message), 
                size = 4, color = "red") +
      theme_void() +
      labs(title = "Plot Generation Error")
  })
}
