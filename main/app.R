#
# CCRD RNA-Seq Portal - Shiny Web Application
# Interactive PCA Plot Viewer
#

library(shiny)
library(DESeq2)
library(foreach)
library(doMC)
library(ggplot2)

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
  
  cat("  - Processing sub-treatment (step 2)...\n")
  coldata$sub.treatment <- substr(x = coldata$sub.treatment, start = ifelse(test = startsWith(x = coldata$sub.treatment, prefix = 'TGFb+'), yes = 6, no = 0), stop = 100)
  
  cat("  - Cleaning sub-treatment labels...\n")
  coldata$sub.treatment[coldata$sub.treatment == 'TNFa'] <- 'none'
  coldata$sub.treatment[coldata$sub.treatment == 'TGFb'] <- 'none'
  coldata$sub.treatment[coldata$sub.treatment == 'TNFa+TGFb'] <- 'none'

  # For FPE6 specific processing
  cat("  - Creating FPE6 specific coldata...\n")
  coldata_fpe6 <- coldata
  coldata_fpe6$sub.treatment[coldata_fpe6$sub.treatment == 'none'] <- 'TNFa'
  cat("✓ Additional metadata processing complete\n")
}, error = function(e) {
  cat(paste("ERROR processing metadata:", e$message, "\n"))
  stop("Failed to process metadata")
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
  titlePanel("CCRD RNA-Seq Portal - PCA Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Plot Selection"),
      selectInput("experiment", 
                  "Choose Experiment:",
                  choices = list(
                    "All Experiments" = "all",
                    "Experiment 4" = "fpe4", 
                    "Experiment 5" = "fpe5",
                    "Experiment 6" = "fpe6",
                    "Experiment 7" = "fpe7"
                  ),
                  selected = "all"),
      hr(),
      p("Select different experiments to view their respective PCA plots."),
      p(strong("All Experiments:"), "Overview of all data"),
      p(strong("Experiment 4-7:"), "Individual fibroblast priming experiments")
    ),
    
    mainPanel(
      width = 9,
      plotOutput("pcaPlot", height = "600px")
    )
  )
)

####################################################################################################
### SHINY SERVER ####################################################################################
####################################################################################################
server <- function(input, output) {
  
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
      
      # Generate the plot
      plot_result <- switch(input$experiment,
             "all" = generate_all_plot(cts, coldata),
             "fpe4" = generate_fpe4_plot(cts, coldata),
             "fpe5" = generate_fpe5_plot(cts, coldata),
             "fpe6" = generate_fpe6_plot(cts, coldata, coldata_fpe6),
             "fpe7" = generate_fpe7_plot(cts, coldata))
      
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
        theme(axis.text = element_blank(),
              axis.ticks = element_blank())
    })
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
