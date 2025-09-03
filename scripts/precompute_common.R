#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
})

log_msg <- function(...) {
  cat(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""), "\n"))
  flush.console()
}

# Resolve project root (parent of this scripts/ folder) robustly
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0) {
  # Fallback to working directory if not run via Rscript
  project_root <- normalizePath("..", winslash = "/", mustWork = FALSE)
} else {
  project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
}

data_dir <- file.path(project_root, "data")
counts_dir <- file.path(data_dir, "gene-counts_FPE")
meta_file <- file.path(data_dir, "metadata.csv")
genes_file <- file.path(data_dir, "2025-03-19_genes.rds")

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

safe_save_rds <- function(object, file) {
  saveRDS(object, file = file, compress = FALSE)
  log_msg("✓ Saved ", basename(file))
}

################################################################################
# Load counts (combine once)                                                   #
################################################################################
log_msg("Loading count matrix from ", counts_dir, " ...")
stopifnot(dir.exists(counts_dir))
files <- list.files(path = counts_dir)
if (length(files) == 0) stop("No count files found in ", counts_dir)

.file <- file.path(counts_dir, files[1])
.cts <- read.table(file = .file, header = TRUE)
genes <- .cts[, 1]

cts <- matrix(nrow = length(genes), ncol = length(files))
rownames(cts) <- genes
colnames(cts) <- sub(pattern = '.gene_id.exon.ct.short.txt', replacement = '', x = files)

for (i in seq_along(files)) {
  f <- file.path(counts_dir, files[i])
  temp_cts <- read.table(file = f, header = TRUE)
  rownames(temp_cts) <- temp_cts[, 1]
  cts[, i] <- temp_cts[genes, 2]
}
log_msg("✓ Count matrix loaded: ", nrow(cts), " genes x ", ncol(cts), " samples")

safe_save_rds(cts, file.path(data_dir, "cts.rds"))

################################################################################
# Load and process coldata                                                     #
################################################################################
log_msg("Loading metadata from ", meta_file, " ...")
stopifnot(file.exists(meta_file))
coldata <- read.csv(file = meta_file)
rownames(coldata) <- coldata$experiment_rna_short_read_id
coldata <- coldata[colnames(cts), ]

# Base columns
coldata$FPE.num <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 1)
coldata$source <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 2)
coldata$source <- sub(pattern = 'MGD1679.', replacement = '', x = coldata$source)
coldata$treatment <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 3)
coldata$TNFa.positive <- grepl(pattern = 'TNFa', x = coldata$treatment)
coldata$TGFb.positive <- grepl(pattern = 'TGFb', x = coldata$treatment)
coldata$treatment.time <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 4)
coldata$replicate.num <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 5)
coldata$none <- grepl(pattern = 'none', x = coldata$treatment)

# Co-treatment label
coldata$co.treatment <- ifelse(
  test = coldata$TNFa.positive & coldata$TGFb.positive,
  yes = 'TNFa+TGFb',
  no = ifelse(test = coldata$TNFa.positive, yes = 'TNFa',
              no = ifelse(test = coldata$TGFb.positive, yes = 'TGFb', no = 'none')))

# Sub-treatment: strip only the leading 'TNFa+'; keep TGFb-derived combinations
coldata$sub.treatment <- substr(
  x = coldata$treatment,
  start = ifelse(test = startsWith(x = coldata$treatment, prefix = 'TNFa+'), yes = 6, no = 0),
  stop = 100)

# FPE6 specific coldata: derive sub.treatment similarly, map 'none' to 'TNFa'
coldata_fpe6 <- coldata
coldata_fpe6$sub.treatment <- substr(
  x = coldata_fpe6$treatment,
  start = ifelse(test = startsWith(x = coldata_fpe6$treatment, prefix = 'TNFa+'), yes = 6, no = 0),
  stop = 100)
coldata_fpe6$sub.treatment[coldata_fpe6$sub.treatment == 'none'] <- 'TNFa'

safe_save_rds(coldata,      file.path(data_dir, "coldata_processed.rds"))
safe_save_rds(coldata_fpe6, file.path(data_dir, "coldata_fpe6.rds"))

################################################################################
# Precompute per-experiment VST and PCA                                        #
################################################################################

compute_vsd_pca <- function(expr_name, sample_idx, intgroup, design_all = FALSE) {
  log_msg("Precomputing ", expr_name, " VST/PCA ...")
  if (!any(sample_idx)) {
    log_msg("  - Skipping ", expr_name, " (no samples)")
    return(invisible(NULL))
  }
  dds <- DESeqDataSetFromMatrix(
    countData = cts[, sample_idx],
    colData  = coldata[sample_idx, , drop = FALSE],
    design   = if (design_all) ~ treatment + treatment.time else ~ treatment
  )
  vsd <- vst(dds, blind = FALSE)
  pcaData <- plotPCA(
    object = vsd,
    intgroup = intgroup,
    returnData = TRUE
  )
  percentVar <- attr(pcaData, 'percentVar')
  safe_save_rds(vsd,     file.path(data_dir, paste0(expr_name, "_vsd.rds")))
  safe_save_rds(list(pcaData = pcaData, percentVar = percentVar),
                file.path(data_dir, paste0(expr_name, "_pca.rds")))
}

# Define experiment indices
idx_all <- rep(TRUE, ncol(cts))
idx_fpe4 <- coldata$FPE.num == 'FPE4'
idx_fpe5 <- coldata$FPE.num == 'FPE5'
idx_fpe6 <- coldata$FPE.num == 'FPE6'

# All experiments
compute_vsd_pca(
  expr_name = "fpe_all",
  sample_idx = idx_all,
  intgroup = c('FPE.num', 'participant_id', 'treatment', 'TNFa.positive', 'treatment.time', 'replicate.num'),
  design_all = TRUE
)

# FPE4
compute_vsd_pca(
  expr_name = "fpe4",
  sample_idx = idx_fpe4,
  intgroup = c('FPE.num', 'participant_id', 'treatment', 'treatment.time', 'replicate.num')
)

# FPE5
compute_vsd_pca(
  expr_name = "fpe5",
  sample_idx = idx_fpe5,
  intgroup = c('FPE.num', 'participant_id', 'treatment', 'treatment.time', 'replicate.num')
)

# FPE6 (use full set including sub.treatment, TNFa.positive, none)
compute_vsd_pca(
  expr_name = "fpe6",
  sample_idx = idx_fpe6,
  intgroup = c('FPE.num', 'participant_id', 'treatment', 'treatment.time', 'replicate.num',
               'sub.treatment', 'TNFa.positive', 'none')
)

# Build gene index map if genes reference exists
if (file.exists(genes_file)) {
  log_msg("Building gene index map from ", basename(genes_file), " ...")
  genes <- readRDS(genes_file)
  # Ensembl -> index
  ensembl_ids <- unique(genes$ensembl_gene_id)
  ensembl_to_index <- vapply(ensembl_ids, function(eid) {
    w <- which(startsWith(rownames(cts), eid))
    if (length(w) > 0) w[1] else NA_integer_
  }, integer(1))
  # Symbol -> index (via first ensembl match for that symbol)
  symbols <- unique(genes$external_gene_name)
  symbol_to_index <- vapply(symbols, function(sym) {
    eids <- genes$ensembl_gene_id[genes$external_gene_name == sym]
    w <- integer(0)
    for (eid in eids) {
      w <- which(startsWith(rownames(cts), eid))
      if (length(w) > 0) break
    }
    if (length(w) > 0) w[1] else NA_integer_
  }, integer(1))
  gene_index_map <- list(
    ensembl_to_index = ensembl_to_index,
    symbol_to_index  = symbol_to_index
  )
  safe_save_rds(gene_index_map, file.path(data_dir, "gene_index_map.rds"))
} else {
  log_msg("⚠️  Genes reference file not found; skipping gene index map")
}

################################################################################
# Legend specifications (levels)                                               #
################################################################################
legend_specs <- list(
  treatment_levels = c('none', 'iBET151', 'TNFa', 'TNFa+iBET151'),
  sub_treatment_levels_fpe6 = c(
    'TNFa', 'rapa', 'TGFb+rapa', 'TGFb', 'TGFb+SB', 'SB', 'rapa+SB',
    'KC7F2', 'TGFb+KC7F2', 'EGFRi', 'EGFRi+EGF', 'EGF'
  )
)
safe_save_rds(legend_specs, file.path(data_dir, "legend_specs.rds"))

log_msg("✅ Precompute complete.")


