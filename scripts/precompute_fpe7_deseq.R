# ===================================================================================================
# PRE-COMPUTE DESEQ ANALYSIS FOR FPE7
# ===================================================================================================
# This script runs the computationally expensive DESeq analysis once and saves the result
# Run this script whenever you update your data to regenerate the pre-computed object

library(DESeq2)
library(foreach)
library(doMC)

cat("=== PRE-COMPUTING DESEQ ANALYSIS FOR FPE7 ===\n")

####################################################################################################
### LOAD DATA ######################################################################################
####################################################################################################

cat("Loading count matrix...\n")
path <- '../data/gene-counts_FPE/'
files <- list.files(path = path)
.file <- paste0(path, files[1])
.cts <- read.table(file = .file, header = TRUE)
genes <- .cts[,1]

# Initialize count matrix
cts <- matrix(nrow = length(genes), ncol = length(files))
rownames(cts) <- genes
colnames(cts) <- sub(pattern = '.gene_id.exon.ct.short.txt', replacement = '', x = files)

# Load all files sequentially
for(i in 1:length(files)) {
  file <- paste0(path, files[i])
  temp_cts <- read.table(file = file, header = TRUE)
  rownames(temp_cts) <- temp_cts[,1]
  cts[,i] <- temp_cts[genes, 2]
}

cat("Loading metadata...\n")
coldata <- read.csv(file = '../data/metadata.csv')
rownames(coldata) <- coldata$experiment_rna_short_read_id
coldata <- coldata[colnames(cts),]
coldata$FPE.num <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 1)
coldata$treatment <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 3)
coldata$TNFa.positive <- grepl(pattern = 'TNFa', x = coldata$treatment)
coldata$treatment.time <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 4)
coldata$replicate.num <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 5)

####################################################################################################
### RUN DESEQ ANALYSIS FOR FPE7 ####################################################################
####################################################################################################

cat("Filtering samples for FPE7...\n")
selected.samples <- coldata$FPE.num == 'FPE7'
fpe7_samples <- sum(selected.samples)
cat(paste("Found", fpe7_samples, "FPE7 samples\n"))

if(fpe7_samples == 0) {
  stop("No FPE7 samples found in the data!")
}

cat("Creating DESeq2 dataset...\n")
dds <- DESeqDataSetFromMatrix(
  countData = cts[, selected.samples], 
  colData = coldata[selected.samples, ], 
  design = ~ treatment)

cat("Running DESeq analysis (this may take several minutes)...\n")
start_time <- Sys.time()
dds_analyzed <- DESeq(dds)
end_time <- Sys.time()
cat(paste("DESeq analysis completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n"))

####################################################################################################
### SAVE PRE-COMPUTED RESULTS ######################################################################
####################################################################################################

output_dir <- "../data/"
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_file <- paste0(output_dir, "fpe7_deseq_precomputed.rds")
cat(paste("Saving pre-computed DESeq results to:", output_file, "\n"))

saveRDS(dds_analyzed, file = output_file)

# Also save some metadata about this run
analysis_info <- list(
  creation_date = Sys.time(),
  r_version = R.version.string,
  deseq2_version = packageVersion("DESeq2"),
  sample_count = fpe7_samples,
  gene_count = nrow(dds_analyzed),
  analysis_time_minutes = as.numeric(difftime(end_time, start_time, units = "mins"))
)

info_file <- paste0(output_dir, "fpe7_deseq_info.rds")
saveRDS(analysis_info, file = info_file)

cat("\n=== PRE-COMPUTATION COMPLETE ===\n")
cat("Files created:\n")
cat(paste("- DESeq object:", output_file, "\n"))
cat(paste("- Analysis info:", info_file, "\n"))
cat(paste("- Analysis time:", round(analysis_info$analysis_time_minutes, 2), "minutes\n"))
cat(paste("- Samples analyzed:", analysis_info$sample_count, "\n"))
cat(paste("- Genes analyzed:", analysis_info$gene_count, "\n"))
cat("\nYou can now use the fast-loading Shiny app!\n")
