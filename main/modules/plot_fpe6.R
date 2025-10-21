# ===================================================================================================
# EXPERIMENT 6 (FPE6) PCA PLOT MODULE
# ===================================================================================================
# This module contains the plot generation function for Fibroblast Priming Experiment #6
# Shows samples colored by sub-treatment and shaped by TNFa status, with special highlighting

generate_fpe6_plot <- function(cts, coldata, coldata_fpe6, fpe6_dds = NULL, gene_name = "CXCL8", genes = NULL) {
  tryCatch({
    # Filter samples for FPE6 only
    selected.samples <- coldata$FPE.num == 'FPE6'
    
    # Check if we have FPE6 samples
    if(sum(selected.samples) == 0) {
      stop("No FPE6 samples found in the data")
    }
    
    # Use the special FPE6 processed coldata
    coldata_subset <- coldata_fpe6[selected.samples, ]
    
    # Use pre-computed DESeq object if available, otherwise run analysis
    if(!is.null(fpe6_dds)) {
      # Use pre-computed DESeq results for much faster loading
      cat("  Using pre-computed DESeq analysis for FPE6\n")
      cat(paste("  Pre-computed object has", ncol(fpe6_dds), "samples and", nrow(fpe6_dds), "genes\n"))
      
      # Use pre-computed DESeq results directly (we trust it's complete since we just created it)
      cat("  ✓ Using pre-computed DESeq analysis for FPE6\n")
      dds <- fpe6_dds
    } else {
      # Fall back to running analysis on-demand (slower)
      cat("  Running DESeq analysis on-demand (slower - consider pre-computing)\n")
      
      # Create DESeq2 dataset with FPE6 samples only
      dds <- DESeqDataSetFromMatrix(
        countData = cts[, selected.samples], 
        colData = coldata_subset, 
        design = ~ treatment)
      
      # Run DESeq2 analysis (needed for gene count plots)
      dds <- DESeq(dds)
    }
    
    # Use precomputed PCA if available, otherwise generate it
    pca_file <- '../data/fpe6_pca.rds'
    if (file.exists(pca_file)) {
      pre <- readRDS(pca_file)
      pcaData <- pre$pcaData
      percentVar <- round(100 * pre$percentVar)
    } else {
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
    
    # Reorder factor levels for consistent plotting
    pcaData$sub.treatment <- factor(pcaData$sub.treatment, 
                                     levels = c('TNFa', 'rapa', 'TGFb+rapa', 'TGFb', 'TGFb+SB', 
                                               'SB', 'rapa+SB', 'KC7F2', 'TGFb+KC7F2', 
                                               'EGFRi', 'EGFRi+EGF', 'EGF'))
    pcaData$TNFa.positive <- factor(pcaData$TNFa.positive, levels = c(TRUE, FALSE))
    
    # Create PCA plot
    pca_plot <- ggplot(data = pcaData, aes(x = PC1, y = PC2, color = sub.treatment, shape = TNFa.positive)) +
      geom_point(size = 3) +
      geom_point(data = pcaData[pcaData$none, ], size = 3, color = 'black') +
      xlab(paste0('PC1: ', percentVar[1], '% variance')) +
      ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
      coord_fixed() +
      scale_shape_discrete(labels = c("Yes", "No")) +
      scale_color_brewer(palette = 'Set3', labels = c("TNFa", "Rapamycin", "TGFb+Rapa", "TGFb", "TGFb+SB", "SB", "Rapa+SB", "KC7F2", "TGFb+KC7F2", "EGFRi", "EGFRi+EGF", "EGF"), na.translate = FALSE) +
      ggtitle(label = 'Fibroblast Priming Experiment 6') +
      theme_gray(base_size = 20) +
      guides(
        shape = guide_legend(title = "TNF?", order = 1),
        color = guide_legend(title = "Sub-treatment", order = 2)
      )
    
    #########################################################################################
    ### GENE COUNT PLOT ####################################################################
    #########################################################################################
    
    # Handle both gene names and Ensembl gene IDs
    user_input <- gene_name
    if(is.null(user_input) || user_input == "" || is.na(user_input)) {
      user_input <- "DMD"  # Default gene for FPE6
      cat(paste("  Gene input was empty, using default:", user_input, "\n"))
    }
    cat(paste("  Searching for gene:", user_input, "\n"))
    
    # Check if genes reference is available
    if(is.null(genes)) {
      cat("  ⚠️ Genes reference not available, using most highly expressed gene\n")
      fpe6_counts <- cts[, selected.samples]
      gene_means <- rowMeans(fpe6_counts)
      selected_gene <- names(sort(gene_means, decreasing = TRUE))[1]
      gene <- which(rownames(cts) == selected_gene)
      display_name <- paste(user_input, "(genes reference unavailable - showing top gene)")
      cat(paste("  Fallback gene:", selected_gene, "at index:", gene, "\n"))
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
        # Find matching Ensembl IDs for the gene name
        matching_ensembl_ids <- genes$ensembl_gene_id[genes$external_gene_name == gene.name]
        
        if(length(matching_ensembl_ids) > 0) {
          # Look for any of the matching Ensembl IDs in the count matrix
          gene <- NULL
          for(ensembl_id in matching_ensembl_ids) {
            if(!is.na(ensembl_id) && ensembl_id != "") {
              gene_match <- which(startsWith(x = rownames(cts), prefix = ensembl_id))
              if(length(gene_match) > 0) {
                gene <- gene_match[1]  # Use first match
                selected_gene <- rownames(cts)[gene]
                cat(paste("  ✓ Found gene in count matrix:", selected_gene, "\n"))
                display_name <- gene.name
                break
              }
            }
          }
          
          if(is.null(gene)) {
            gene <- NULL
          }
        } else {
          gene <- NULL
        }
      } else {
        gene <- NULL
      }
      
      # Fallback if gene not found
      if(is.null(gene) || length(gene) == 0) {
        cat(paste("  ⚠️ Gene", user_input, "not found, using most highly expressed gene\n"))
        fpe6_counts <- cts[, selected.samples]
        gene_means <- rowMeans(fpe6_counts)
        selected_gene <- names(sort(gene_means, decreasing = TRUE))[1]
        gene <- which(rownames(cts) == selected_gene)
        display_name <- paste(user_input, "(not found - showing top gene)")
        cat(paste("  Fallback gene:", selected_gene, "at index:", gene, "\n"))
      }
    }
    
    # Generate count data for the selected gene (use gene index)
    cat(paste("  Gene index:", gene, "\n"))
    cat(paste("  Gene name from cts:", rownames(cts)[gene], "\n"))
    
    if(is.null(gene) || length(gene) == 0 || is.na(gene)) {
      stop("Gene index is invalid")
    }
    
    # Find the corresponding gene in the DESeq object
    gene_name_in_dds <- rownames(cts)[gene]
    gene_index_in_dds <- which(rownames(dds) == gene_name_in_dds)
    
    cat(paste("  Gene index in DESeq object:", gene_index_in_dds, "\n"))
    
    if(length(gene_index_in_dds) == 0) {
      stop(paste("Gene", gene_name_in_dds, "not found in DESeq object"))
    }
    
    count_data <- plotCounts(dds, gene = gene_index_in_dds, intgroup = 'treatment', returnData = TRUE)
    count_data$experiment_rna_short_read_id <- rownames(count_data)
    count_data <- merge(count_data, coldata_subset)
    
    # Apply FPE6-specific sub.treatment processing
    treatment_char <- as.character(count_data$treatment)
    count_data$sub.treatment <- substr(x = treatment_char, start = ifelse(test = startsWith(x = treatment_char, prefix = 'TNFa+'), yes = 6, no = 0), stop = 100)
    count_data$sub.treatment[count_data$sub.treatment == 'none'] <- 'TNFa'
    count_data$TNFa.positive <- grepl(pattern = 'TNFa', x = treatment_char)
    
    # Reorder sub.treatment factor levels
    count_data$sub.treatment <- factor(count_data$sub.treatment, 
                                       levels = c('TNFa', 'rapa', 'TGFb+rapa', 'TGFb', 'TGFb+SB', 
                                                 'SB', 'rapa+SB', 'KC7F2', 'TGFb+KC7F2', 
                                                 'EGFRi', 'EGFRi+EGF', 'EGF'))
    
    # Create gene count plot
    count_plot <- ggplot(count_data, aes(x = sub.treatment, y = count, color = TNFa.positive)) + 
      geom_boxplot(outlier.shape = NA, color = 'black', alpha = 0.6) +
      geom_point(position = position_jitter(w = 0.15, h = 0), size = 2.5) + 
      ggtitle(label = paste("RNA Abundance -", display_name, "(", selected_gene, ")")) +
      scale_color_manual(values = c('#f1a340', '#998ec3'), labels = c("No", "Yes")) +
      scale_x_discrete(labels = c("TNFa", "Rapamycin", "TGFb+Rapa", "TGFb", "TGFb+SB", "SB", "Rapa+SB", "KC7F2", "TGFb+KC7F2", "EGFRi", "EGFRi+EGF", "EGF")) +
      scale_y_log10(limits = c(1, NA)) +
      theme_gray(base_size = 20) +
      labs(
        x = "",
        y = "Normalized Count (log10)",
        color = "TNF?"
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
      geom_text(label = paste("Error in FPE6 plot:", e$message), 
                size = 4, color = "red") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "gray95", color = NA)) +
      labs(title = "Plot Generation Error")
  })
}
