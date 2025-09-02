# ===================================================================================================
# EXPERIMENT 7 (FPE7) PCA PLOT MODULE
# ===================================================================================================
# This module contains the plot generation function for Fibroblast Priming Experiment #7
# Shows PCA plot and interactive gene expression plot with user-selectable genes
# Supports both gene symbols (e.g., 'DMD') and Ensembl IDs (e.g., 'ENSG00000198947')
# Uses genes reference database for bidirectional gene symbol ↔ Ensembl ID mapping

generate_fpe7_plot <- function(cts, coldata, fpe7_dds = NULL, gene_name = "DMD", genes = NULL) {
  tryCatch({
    # Filter samples for FPE7 only
    selected.samples <- coldata$FPE.num == 'FPE7'
    
    # Check if we have FPE7 samples
    if(sum(selected.samples) == 0) {
      stop("No FPE7 samples found in the data")
    }
    
    # Use pre-computed DESeq object if available, otherwise run analysis
    if(!is.null(fpe7_dds)) {
      # Use pre-computed DESeq results for much faster loading
      cat("  Using pre-computed DESeq analysis for FPE7\n")
      dds <- fpe7_dds
      
      # Apply variance stabilizing transformation for PCA
      vsd <- vst(dds, blind = FALSE)
    } else {
      # Fall back to running analysis on-demand (slower)
      cat("  Running DESeq analysis on-demand (slower - consider pre-computing)\n")
      
      # Create DESeq2 dataset with FPE7 samples only
      dds <- DESeqDataSetFromMatrix(
        countData = cts[, selected.samples], 
        colData = coldata[selected.samples, ], 
        design = ~ treatment)
      
      # Run DESeq2 analysis (needed for gene count plots)
      dds <- DESeq(dds)
      
      # Apply variance stabilizing transformation for PCA
      vsd <- vst(dds, blind = FALSE)
    }
    
    #########################################################################################
    ### PCA PLOT ############################################################################
    #########################################################################################
    
    # Generate PCA data
    pcaData <- plotPCA(
      object = vsd, 
      intgroup = c('FPE.num', 'participant_id', 'affected_status', 'treatment', 'TNFa.positive', 'treatment.time', 'replicate.num'), 
      returnData = TRUE)
    
    # Factor treatment levels for consistent ordering
    pcaData$treatment <- factor(x = pcaData$treatment, levels = c('none', 'iBET151', 'TNFa', 'TNFa+iBET151'))
    percentVar <- round(100 * attr(pcaData, 'percentVar'))
    
    # Create PCA plot
    pca_plot <- ggplot(data = pcaData, aes(x = PC1, y = PC2, fill = participant_id, shape = treatment, color = participant_id)) +
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
      ggtitle(label = 'PCA - Fibroblast Priming #7') +
      theme_minimal() +
      theme(
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "right"
      ) +
      guides(
        fill = guide_legend(title = "Participant", override.aes = list(shape = 21)),
        color = guide_legend(title = "Participant", override.aes = list(shape = 21)),
        shape = guide_legend(title = "Treatment")
      )
    
    #########################################################################################
    ### GENE COUNT PLOT ####################################################################
    #########################################################################################
    
    # Handle both gene names and Ensembl gene IDs
    user_input <- gene_name
    cat(paste("  Searching for gene:", user_input, "\n"))
    
    # Check if genes reference is available
    if(is.null(genes)) {
      cat("  ⚠️ Genes reference not available, using most highly expressed gene\n")
      fpe7_counts <- cts[, selected.samples]
      gene_means <- rowMeans(fpe7_counts)
      selected_gene <- names(sort(gene_means, decreasing = TRUE))[1]
      gene <- which(rownames(cts) == selected_gene)
      display_name <- paste(user_input, "(genes reference unavailable - showing top gene)")
    } else {
      # Check if input looks like an Ensembl gene ID (starts with ENSG)
      if(grepl("^ENSG", user_input, ignore.case = TRUE)) {
        # Input is an Ensembl gene ID - look up the corresponding gene name
        ensembl.gene.id <- user_input
        gene.name <- genes$external_gene_name[genes$ensembl_gene_id == ensembl.gene.id]
        
        if(length(gene.name) > 0 && !is.na(gene.name[1]) && gene.name[1] != "") {
          gene.name <- gene.name[1]  # Use first match
          cat(paste("  ✓ Ensembl ID", ensembl.gene.id, "maps to gene:", gene.name, "\n"))
        } else {
          cat(paste("  ⚠️ Ensembl ID", ensembl.gene.id, "not found in genes reference\n"))
          gene.name <- NULL
        }
      } else {
        # Input is a gene name - use directly
        gene.name <- user_input
      }
      
      # Find the gene using the genes mapping (if we have a valid gene name)
      if(!is.null(gene.name)) {
        gene <- which(startsWith(
          x = rownames(cts), 
          prefix = genes$ensembl_gene_id[genes$external_gene_name == gene.name]))
        
        if(length(gene) > 0) {
          # Use the first match if multiple found
          gene <- gene[1]  # Use first index
          selected_gene <- rownames(cts)[gene]
          cat(paste("  ✓ Found gene in count matrix:", selected_gene, "\n"))
          display_name <- gene.name
        } else {
          gene <- NULL
        }
      } else {
        gene <- NULL
      }
      
      # Fallback if gene not found
      if(is.null(gene) || length(gene) == 0) {
        cat(paste("  ⚠️ Gene", user_input, "not found, using most highly expressed gene\n"))
        fpe7_counts <- cts[, selected.samples]
        gene_means <- rowMeans(fpe7_counts)
        selected_gene <- names(sort(gene_means, decreasing = TRUE))[1]
        gene <- which(rownames(cts) == selected_gene)
        display_name <- paste(user_input, "(not found - showing top gene)")
      }
    }
    
    # Generate count data for the selected gene (use gene index)
    count_data <- plotCounts(dds, gene = gene, intgroup = 'treatment', returnData = TRUE)
    count_data$experiment_rna_short_read_id <- rownames(count_data)
    count_data <- merge(count_data, coldata[selected.samples, ])
    count_data$treatment <- factor(x = count_data$treatment, levels = c('none', 'iBET151', 'TNFa', 'TNFa+iBET151'))
    
    # Create gene count plot
    count_plot <- ggplot(count_data, aes(x = treatment, y = count, color = participant_id)) + 
      geom_boxplot(outlier.shape = NA, color = 'black', alpha = 0.6) +
      geom_point(position = position_jitter(w = 0.15, h = 0), size = 2.5) + 
      ggtitle(label = paste("Gene Expression -", display_name)) +
      scale_color_brewer(palette = 'Dark2') +
      scale_y_log10(limits = c(1, NA)) +
      theme_minimal() +
      theme(
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      ) +
      labs(
        x = "Treatment",
        y = "Normalized Count (log10)",
        color = "Participant"
      )
    
    #########################################################################################
    ### COMBINE PLOTS #######################################################################
    #########################################################################################
    
    # Combine both plots vertically
    combined_plot <- grid.arrange(pca_plot, count_plot, ncol = 1, heights = c(1, 1))
    
    return(combined_plot)
    
  }, error = function(e) {
    # Return error plot if generation fails
    ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
      geom_text(label = paste("Error in FPE7 plot:", e$message), 
                size = 4, color = "red") +
      theme_void() +
      labs(title = "Plot Generation Error")
  })
}
