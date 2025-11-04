#!/usr/bin/env Rscript
# ===================================================================================================
# PRE-COMPUTE FPE6 DESEQ ANALYSIS
# ===================================================================================================
# This script pre-computes the DESeq2 analysis for FPE6 to speed up the Shiny app
# Run this script to generate fpe6_deseq_precomputed.rds and fpe6_deseq_info.rds

library(DESeq2)
library(foreach)
library(doParallel)

# Set up parallel processing
cl <- makeCluster(8)
registerDoParallel(cl)

cat("Starting FPE6 DESeq2 pre-computation...\n")
cat("=====================================\n")

####################################################################################################
### READ COUNT MATRIX ##############################################################################
####################################################################################################
cat("Loading count matrix...\n")
path <- 'data/gene-counts_FPE/'
files <- list.files(path = path)
.file <- paste0(path, files[1])
.cts <- read.table(file = .file, header = TRUE)
genes <- .cts[,1]

# Initialize count matrix
cts <- matrix(nrow = length(genes), ncol = length(files))
rownames(cts) <- genes
colnames(cts) <- sub(pattern = '.gene_id.exon.ct.short.txt', replacement = '', x = files)

# Load all files sequentially (more stable than parallel for this operation)
for(i in 1:length(files)) {
  file <- paste0(path, files[i])
  temp_cts <- read.table(file = file, header = TRUE)
  rownames(temp_cts) <- temp_cts[,1]
  cts[,i] <- temp_cts[genes, 2]
}

cat(paste("✓ Count matrix loaded:", nrow(cts), "genes x", ncol(cts), "samples\n"))

####################################################################################################
### READ COLDATA ###################################################################################
####################################################################################################
cat("Loading metadata...\n")
coldata <- read.csv(file = 'data/metadata.csv')
rownames(coldata) <- coldata$experiment_rna_short_read_id
coldata <- coldata[colnames(cts),]
coldata$FPE.num <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 1)
coldata$treatment <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 3)
coldata$TNFa.positive <- grepl(pattern = 'TNFa', x = coldata$treatment)
coldata$treatment.time <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 4)
coldata$replicate.num <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 5)

# Additional data processing for FPE6
cat("Processing FPE6-specific metadata...\n")
coldata$source <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 2)
coldata$source <- sub(pattern = 'MGD1679.', replacement = '', x = coldata$source)

# For FPE6 specific processing (retain TGFb-derived sub-treatments; only strip TNFa+ prefix)
coldata$sub.treatment <- substr(x = coldata$treatment, start = ifelse(test = startsWith(x = coldata$treatment, prefix = 'TNFa+'), yes = 6, no = 0), stop = 100)
coldata$sub.treatment[coldata$sub.treatment == 'none'] <- 'TNFa'
coldata$none <- grepl(pattern = 'none', x = coldata$treatment)

cat(paste("✓ Metadata loaded:", nrow(coldata), "samples\n"))

####################################################################################################
### FILTER FPE6 SAMPLES ############################################################################
####################################################################################################
cat("Filtering FPE6 samples...\n")
selected.samples <- coldata$FPE.num == 'FPE6'

if(sum(selected.samples) == 0) {
  stop("No FPE6 samples found in the data")
}

cat(paste("✓ Found", sum(selected.samples), "FPE6 samples\n"))

####################################################################################################
### DESeq2 ANALYSIS ################################################################################
####################################################################################################
cat("Creating DESeq2 dataset...\n")
dds <- DESeqDataSetFromMatrix(
  countData = cts[, selected.samples], 
  colData = coldata[selected.samples, ], 
  design = ~ treatment)

cat("Running DESeq2 analysis (this may take several minutes)...\n")
cat("  - Estimating size factors...\n")
cat("  - Estimating dispersions...\n")
cat("  - Gene-wise dispersion estimates...\n")
cat("  - Final dispersion estimates...\n")
cat("  - Fitting model and testing...\n")

start_time <- Sys.time()
dds <- DESeq(dds)
end_time <- Sys.time()

# Verify the analysis is complete
if(is.null(metadata(dds)$betaPriorVar)) {
  cat("⚠️ DESeq analysis appears incomplete, trying again...\n")
  dds <- DESeq(dds)
}

cat(paste("✓ DESeq2 analysis completed in", round(as.numeric(end_time - start_time, units = "mins"), 2), "minutes\n"))
cat(paste("✓ Analysis includes", ncol(dds), "samples and", nrow(dds), "genes\n"))

####################################################################################################
### SAVE RESULTS ####################################################################################
####################################################################################################
cat("Saving pre-computed results...\n")

# Save the DESeq object
saveRDS(dds, file = 'data/fpe6_deseq_precomputed.rds')
cat("✓ Saved DESeq object to data/fpe6_deseq_precomputed.rds\n")

# Save analysis info
analysis_info <- list(
  sample_count = sum(selected.samples),
  gene_count = nrow(cts),
  creation_date = as.character(Sys.time()),
  experiment = "FPE6",
  design_formula = "~ treatment"
)

saveRDS(analysis_info, file = 'data/fpe6_deseq_info.rds')
cat("✓ Saved analysis info to data/fpe6_deseq_info.rds\n")

# Clean up parallel processing
stopCluster(cl)

cat("\n=====================================\n")
cat("FPE6 DESeq2 pre-computation complete!\n")
cat("=====================================\n")
cat(paste("Analysis info:\n"))
cat(paste("  - Samples:", analysis_info$sample_count, "\n"))
cat(paste("  - Genes:", analysis_info$gene_count, "\n"))
cat(paste("  - Created:", analysis_info$creation_date, "\n"))
cat(paste("  - Design:", analysis_info$design_formula, "\n"))
cat("\nThe Shiny app will now load FPE6 plots much faster!\n")
