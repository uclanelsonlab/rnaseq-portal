#
# CCRD RNA-Seq Portal - Shiny Web Application
# Interactive PCA Plot Viewer
#

library(shiny)
library(DESeq2)
library(foreach)
library(doMC)
library(ggplot2)
library(gridExtra)

# Note: Parallel processing removed for Shiny app stability

# Load data once at startup
####################################################################################################
### READ COUNT MATRIX ##############################################################################
####################################################################################################
tryCatch({
  cat("Loading count matrix...\n")
  path <- '../data/gene-counts_FPE/'
  files <- list.files(path = path)
  if(length(files) == 0) {
    stop("No count files found in data directory")
  }
  
  # Load first file to get gene names
  .file <- paste0(path, files[1])
  .cts <- read.table(file = .file, header = TRUE)
  genes <- .cts[,1]
  
  # Initialize count matrix
  cts <- matrix(nrow = length(genes), ncol = length(files))
  rownames(cts) <- genes
  colnames(cts) <- sub(pattern = '.gene_id.exon.ct.short.txt', replacement = '', x = files)
  
  # Load all files sequentially (more stable than parallel for Shiny)
  for(i in 1:length(files)) {
    file <- paste0(path, files[i])
    temp_cts <- read.table(file = file, header = TRUE)
    rownames(temp_cts) <- temp_cts[,1]
    cts[,i] <- temp_cts[genes, 2]
  }
  
  cat(paste("✓ Count matrix loaded:", nrow(cts), "genes x", ncol(cts), "samples\n"))
}, error = function(e) {
  cat(paste("ERROR loading count matrix:", e$message, "\n"))
  stop("Failed to load count matrix")
})

####################################################################################################
### READ COLDATA ###################################################################################
####################################################################################################
tryCatch({
  cat("Loading metadata...\n")
  coldata <- read.csv(file = '../data/metadata.csv')
  rownames(coldata) <- coldata$experiment_rna_short_read_id
  coldata <- coldata[colnames(cts),]
  coldata$FPE.num <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 1)
  coldata$treatment <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 3)
  coldata$TNFa.positive <- grepl(pattern = 'TNFa', x = coldata$treatment)
  coldata$treatment.time <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 4)
  coldata$replicate.num <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 5)
  cat(paste("✓ Metadata loaded:", nrow(coldata), "samples\n"))
}, error = function(e) {
  cat(paste("ERROR loading metadata:", e$message, "\n"))
  stop("Failed to load metadata")
})

tryCatch({
  cat("Processing additional metadata columns...\n")
  # Additional data processing for individual experiments
  cat("  - Processing source column...\n")
  coldata$source <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 2)
  coldata$source <- sub(pattern = 'MGD1679.', replacement = '', x = coldata$source)

  # For FPE5 and FPE6 specific processing
  cat("  - Processing treatment flags...\n")
  coldata$none <- grepl(pattern = 'none', x = coldata$treatment)
  coldata$TGFb.positive <- grepl(pattern = 'TGFb', x = coldata$treatment)
  
  cat("  - Processing co-treatment...\n")
  coldata$co.treatment <- ifelse(
    test = coldata$TNFa.positive & coldata$TGFb.positive, 
    yes = 'TNFa+TGFb', 
    no = ifelse(test = coldata$TNFa.positive, yes = 'TNFa', no = ifelse(test = coldata$TGFb.positive, yes = 'TGFb', no = 'none')))
  
  cat("  - Processing sub-treatment (step 1)...\n")
  coldata$sub.treatment <- substr(x = coldata$treatment, start = ifelse(test = startsWith(x = coldata$treatment, prefix = 'TNFa+'), yes = 6, no = 0), stop = 100)
  
  # For FPE6 specific processing (retain TGFb-derived sub-treatments; only strip TNFa+ prefix)
  cat("  - Creating FPE6 specific coldata...\n")
  coldata_fpe6 <- coldata
  coldata_fpe6$sub.treatment <- substr(x = coldata_fpe6$treatment, start = ifelse(test = startsWith(x = coldata_fpe6$treatment, prefix = 'TNFa+'), yes = 6, no = 0), stop = 100)
  coldata_fpe6$sub.treatment[coldata_fpe6$sub.treatment == 'none'] <- 'TNFa'
  cat("✓ Additional metadata processing complete\n")
}, error = function(e) {
  cat(paste("ERROR processing metadata:", e$message, "\n"))
  stop("Failed to process metadata")
})

# Load pre-computed FPE7 DESeq object for faster plot generation
tryCatch({
  cat("Loading pre-computed FPE7 DESeq analysis...\n")
  fpe7_deseq_file <- '../data/fpe7_deseq_precomputed.rds'
  
  if(file.exists(fpe7_deseq_file)) {
    fpe7_dds <- readRDS(fpe7_deseq_file)
    
    # Load analysis info
    info_file <- '../data/fpe7_deseq_info.rds'
    if(file.exists(info_file)) {
      fpe7_info <- readRDS(info_file)
      cat(paste("✓ Pre-computed FPE7 DESeq loaded (", fpe7_info$sample_count, "samples,", 
                fpe7_info$gene_count, "genes)\n"))
      cat(paste("  Analysis created:", fpe7_info$creation_date, "\n"))
    } else {
      cat("✓ Pre-computed FPE7 DESeq loaded\n")
    }
  } else {
    cat("⚠️  Pre-computed FPE7 DESeq not found - will run analysis on-demand (slower)\n")
    cat("   Run 'Rscript ../scripts/precompute_fpe7_deseq.R' to generate it\n")
    fpe7_dds <- NULL
  }
}, error = function(e) {
  cat(paste("⚠️  Error loading pre-computed FPE7 DESeq:", e$message, "\n"))
  cat("   Will fall back to on-demand analysis (slower)\n")
  fpe7_dds <- NULL
})

# Load pre-computed FPE6 DESeq object for faster plot generation
tryCatch({
  cat("Loading pre-computed FPE6 DESeq analysis...\n")
  fpe6_deseq_file <- '../data/fpe6_deseq_precomputed.rds'
  
  if(file.exists(fpe6_deseq_file)) {
    fpe6_dds <- readRDS(fpe6_deseq_file)
    
    # Load analysis info
    info_file <- '../data/fpe6_deseq_info.rds'
    if(file.exists(info_file)) {
      fpe6_info <- readRDS(info_file)
      cat(paste("✓ Pre-computed FPE6 DESeq loaded (", fpe6_info$sample_count, "samples,", 
                fpe6_info$gene_count, "genes)\n"))
      cat(paste("  Analysis created:", fpe6_info$creation_date, "\n"))
    } else {
      cat("✓ Pre-computed FPE6 DESeq loaded\n")
    }
  } else {
    cat("⚠️  Pre-computed FPE6 DESeq not found - will run analysis on-demand (slower)\n")
    cat("   Run 'Rscript ../scripts/precompute_fpe6_deseq.R' to generate it\n")
    fpe6_dds <- NULL
  }
}, error = function(e) {
  cat(paste("⚠️  Error loading pre-computed FPE6 DESeq:", e$message, "\n"))
  cat("   Will fall back to on-demand analysis (slower)\n")
  fpe6_dds <- NULL
})

# Load genes reference data for FPE7 gene mapping
tryCatch({
  cat("Loading genes reference data...\n")
  genes_file <- '../data/2025-03-19_genes.rds'
  if(file.exists(genes_file)) {
    genes <- readRDS(genes_file)
    cat(paste("✓ Genes reference loaded with", nrow(genes), "entries\n"))
  } else {
    cat("⚠️  Genes reference file not found - FPE7 gene search may be limited\n")
    genes <- NULL
  }
}, error = function(e) {
  cat(paste("⚠️  Error loading genes reference:", e$message, "\n"))
  genes <- NULL
})

# Confirm data loading is complete
cat("✓ All data loaded successfully! Shiny app is ready.\n")
cat(paste("Available FPE numbers:", paste(unique(coldata$FPE.num), collapse=", "), "\n"))

####################################################################################################
### LOAD PLOT MODULES ###############################################################################
####################################################################################################

cat("Loading plot modules...\n")

# Source individual plot modules
source("modules/plot_all.R")
cat("✓ All experiments plot module loaded\n")

source("modules/plot_fpe4.R") 
cat("✓ FPE4 plot module loaded\n")

source("modules/plot_fpe5.R")
cat("✓ FPE5 plot module loaded\n")

source("modules/plot_fpe6.R")
cat("✓ FPE6 plot module loaded\n")

source("modules/plot_fpe7.R")
cat("✓ FPE7 plot module loaded\n")

cat("✓ All plot modules loaded successfully!\n")

####################################################################################################
### SHINY UI ########################################################################################
####################################################################################################
ui <- fluidPage(
  tags$div(
    tags$h2("Nelson Lab"),
    tags$h4("RNA-Seq Portal", style = "margin-top: -10px; color: #666;")
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      tags$div(
        style = "display: flex; align-items: center; gap: 10px;",
        tags$label("Experiment:", style = "white-space: nowrap; margin: 0;"),
        tags$div(
          style = "flex: 1;",
          selectInput("experiment", 
                      label = NULL,
                      choices = list(
                        "All Experiments" = "all",
                        "Fibroblast Priming Experiment 4" = "fpe4", 
                        "Fibroblast Priming Experiment 5" = "fpe5",
                        "Fibroblast Priming Experiment 6" = "fpe6",
                        "Fibroblast Priming Experiment 7" = "fpe7"
                      ),
                      selected = "all")
        )
      ),
      
      # Conditional gene input for FPE6 and FPE7
      conditionalPanel(
        condition = "input.experiment == 'fpe6' || input.experiment == 'fpe7'",
        hr(),
        tags$div(
          style = "display: flex; align-items: center; gap: 10px;",
          tags$label("Gene:", style = "white-space: nowrap; margin: 0;"),
          tags$div(
            style = "flex: 0 1 50%; max-width: 200px;",
            textInput("gene_name", 
                      label = NULL,
                      value = "DMD",
                      placeholder = "Enter gene symbol (e.g., DMD, TTN) or Ensembl ID (e.g., ENSG00000198947)")
          )
        ),
        p(style = "font-size: 11px; color: #666;", 
          "Enter a gene symbol OR Ensembl ID"),
        tags$script("
          $(document).ready(function() {
            $('#experiment').on('change', function() {
              var experiment = $(this).val();
              var defaultValue = experiment === 'fpe6' ? 'DMD' : 'DMD';
              $('#gene_name').val(defaultValue);
            });
          });
        ")
        # p(style = "font-size: 11px; color: #666;", 
        #   "If gene not found, the most highly expressed gene will be shown.")
      ),
      
      hr()
    ),
    
    mainPanel(
      width = 9,
      uiOutput("plotUI")
    )
  )
)

####################################################################################################
### SHINY SERVER ####################################################################################
####################################################################################################
server <- function(input, output) {
  
  # Dynamic plot UI with different heights
  output$plotUI <- renderUI({
    plot_height <- switch(input$experiment,
                         "all" = "600px",
                         "fpe4" = "600px", 
                         "fpe5" = "600px",
                         "fpe6" = "1000px",  # FPE6 now has two plots stacked
                         "fpe7" = "1000px",  # FPE7 has two plots stacked
                         "600px")  # default
    
    plotOutput("pcaPlot", height = plot_height)
  })
  
  output$pcaPlot <- renderPlot({
    tryCatch({
      # Add validation to ensure data is loaded
      if(!exists("cts") || !exists("coldata")) {
        stop("Data not loaded properly")
      }
      
      # Get experiment name for logging
      experiment_name <- switch(input$experiment,
                               "all" = "All Experiments",
                               "fpe4" = "Experiment 4 (FPE4)",
                               "fpe5" = "Experiment 5 (FPE5)",
                               "fpe6" = "Experiment 6 (FPE6)",
                               "fpe7" = "Experiment 7 (FPE7)",
                               "Unknown Experiment")
      
      # Log start of plot generation
      cat(paste("\n🔄 [", Sys.time(), "] Starting to generate plot for:", experiment_name, "\n"))
      flush.console()
      
      # Get gene name for FPE6 and FPE7 (use default if empty)
      gene_name <- if(input$experiment == "fpe6" || input$experiment == "fpe7") {
        if(is.null(input$gene_name) || input$gene_name == "" || is.na(input$gene_name)) {
          "DMD"  # Default gene
        } else {
          trimws(input$gene_name)  # Clean whitespace
        }
      } else {
        NULL
      }
      
      # Generate the plot
      plot_result <- switch(input$experiment,
             "all" = generate_all_plot(cts, coldata),
             "fpe4" = generate_fpe4_plot(cts, coldata),
             "fpe5" = generate_fpe5_plot(cts, coldata),
             "fpe6" = generate_fpe6_plot(cts, coldata, coldata_fpe6, fpe6_dds, gene_name, genes),
             "fpe7" = generate_fpe7_plot(cts, coldata, fpe7_dds, gene_name, genes))
      
      # Log completion of plot generation
      cat(paste("✅ [", Sys.time(), "] Finished generating plot for:", experiment_name, "\n"))
      flush.console()
      
      return(plot_result)
      
    }, error = function(e) {
      # Log error
      cat(paste("❌ [", Sys.time(), "] Error generating plot for:", input$experiment, "- Error:", e$message, "\n"))
      flush.console()
      
      # Return an error plot instead of failing silently
      ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
        geom_text(label = paste("Error:", e$message), size = 6) +
        theme_minimal() +
        labs(title = "Plot Generation Error", 
             subtitle = "Check R console for details") +
        theme(panel.background = element_rect(fill = "gray95", color = NA),
              axis.text = element_blank(),
              axis.ticks = element_blank())
    })
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
