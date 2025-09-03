# ===================================================================================================
# EXPERIMENT 5 (FPE5) PCA PLOT MODULE
# ===================================================================================================
# This module contains the plot generation function for Fibroblast Priming Experiment #5
# Shows samples colored by sub-treatment and shaped by co-treatment, with special highlighting

generate_fpe5_plot <- function(cts, coldata) {
  tryCatch({
    # Filter samples for FPE5 only
    selected.samples <- coldata$FPE.num == 'FPE5'
    
    # Check if we have FPE5 samples
    if(sum(selected.samples) == 0) {
      stop("No FPE5 samples found in the data")
    }
    
    # Use precomputed PCA if available
    pca_file <- '../data/fpe5_pca.rds'
    if (file.exists(pca_file)) {
      pre <- readRDS(pca_file)
      pcaData <- pre$pcaData
      percentVar <- round(100 * pre$percentVar)
      
      # Apply FPE5-specific sub.treatment processing to match original script
      # Convert treatment to character for string processing
      treatment_char <- as.character(pcaData$treatment)
      pcaData$sub.treatment <- substr(x = treatment_char, start = ifelse(test = startsWith(x = treatment_char, prefix = 'TNFa+'), yes = 6, no = 0), stop = 100)
      pcaData$sub.treatment <- substr(x = pcaData$sub.treatment, start = ifelse(test = startsWith(x = pcaData$sub.treatment, prefix = 'TGFb+'), yes = 6, no = 0), stop = 100)
      pcaData$sub.treatment[pcaData$sub.treatment == 'TNFa'] <- 'none'
      pcaData$sub.treatment[pcaData$sub.treatment == 'TGFb'] <- 'none'
      pcaData$sub.treatment[pcaData$sub.treatment == 'TNFa+TGFb'] <- 'none'
      
      # Add co.treatment and none columns
      pcaData$TNFa.positive <- grepl(pattern = 'TNFa', x = treatment_char)
      pcaData$TGFb.positive <- grepl(pattern = 'TGFb', x = treatment_char)
      pcaData$co.treatment <- ifelse(
        test = pcaData$TNFa.positive & pcaData$TGFb.positive, 
        yes = 'TNFa+TGFb', 
        no = ifelse(test = pcaData$TNFa.positive, yes = 'TNFa', no = ifelse(test = pcaData$TGFb.positive, yes = 'TGFb', no = 'none')))
      pcaData$none <- grepl(pattern = 'none', x = treatment_char)
    } else {
      # Create DESeq2 dataset with FPE5 samples only
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
      
      # Apply FPE5-specific sub.treatment processing to match original script
      # Convert treatment to character for string processing
      treatment_char <- as.character(pcaData$treatment)
      pcaData$sub.treatment <- substr(x = treatment_char, start = ifelse(test = startsWith(x = treatment_char, prefix = 'TNFa+'), yes = 6, no = 0), stop = 100)
      pcaData$sub.treatment <- substr(x = pcaData$sub.treatment, start = ifelse(test = startsWith(x = pcaData$sub.treatment, prefix = 'TGFb+'), yes = 6, no = 0), stop = 100)
      pcaData$sub.treatment[pcaData$sub.treatment == 'TNFa'] <- 'none'
      pcaData$sub.treatment[pcaData$sub.treatment == 'TGFb'] <- 'none'
      pcaData$sub.treatment[pcaData$sub.treatment == 'TNFa+TGFb'] <- 'none'
      
      # Add co.treatment and none columns
      pcaData$TNFa.positive <- grepl(pattern = 'TNFa', x = treatment_char)
      pcaData$TGFb.positive <- grepl(pattern = 'TGFb', x = treatment_char)
      pcaData$co.treatment <- ifelse(
        test = pcaData$TNFa.positive & pcaData$TGFb.positive, 
        yes = 'TNFa+TGFb', 
        no = ifelse(test = pcaData$TNFa.positive, yes = 'TNFa', no = ifelse(test = pcaData$TGFb.positive, yes = 'TGFb', no = 'none')))
      pcaData$none <- grepl(pattern = 'none', x = treatment_char)
    }
    
    
    # Create the plot with special highlighting for 'none' samples
    ggplot(data = pcaData, aes(x = PC1, y = PC2, color = sub.treatment, shape = co.treatment)) +
      geom_point(size = 3) +
      geom_point(data = pcaData[pcaData$none, ], size = 3, color = 'black') +
      xlab(paste0('PC1: ', percentVar[1], '% variance')) +
      ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
      coord_fixed() +
      scale_color_brewer(palette = 'Dark2') +
      ggtitle(label = 'Fibroblast Priming #5') +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "gray95", color = NA),
        panel.grid.major = element_line(color = "white", size = 0.5),
        panel.grid.minor = element_line(color = "white", size = 0.25),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      guides(
        color = guide_legend(title = "Sub-treatment"),
        shape = guide_legend(title = "Co-treatment")
      )
    
  }, error = function(e) {
    # Return error plot if generation fails
    ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
      geom_text(label = paste("Error in FPE5 plot:", e$message), 
                size = 4, color = "red") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "gray95", color = NA)) +
      labs(title = "Plot Generation Error")
  })
}
