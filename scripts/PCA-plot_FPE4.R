rm(list = ls());

library(DESeq2);
library(foreach);
library(doMC);
library(ggplot2);

registerDoMC(cores = 8);

####################################################################################################
### READ COUNT MATRIX ##############################################################################
####################################################################################################
path <- 'data/gene-counts_FPE/';
files <- list.files(path = path);
.file <- paste0(path, files[1]);
.cts <- read.table(file = .file, header = TRUE);
genes <- .cts[,1];
cts <- foreach(i=1:length(files), .combine = cbind) %dopar% {
  file <- paste0(path, files[i]);
  cts <- read.table(file = file, header = TRUE);
  rownames(cts) <- cts[,1];
  cts <- cts[genes,2];
}
rownames(cts) <- genes;
colnames(cts) <- sub(pattern = '.gene_id.exon.ct.short.txt', replacement = '', x = files);

####################################################################################################
### READ COLDATA ###################################################################################
####################################################################################################
coldata <- read.csv(file = 'data/metadata.csv');
rownames(coldata) <- coldata$experiment_rna_short_read_id;
coldata <- coldata[ colnames(cts),];
coldata$FPE.num <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 1);
coldata$source <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 2);
coldata$source <- sub(pattern = 'MGD1679.', replacement = '', x = coldata$source);
coldata$treatment <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 3);
coldata$treatment.time <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 4);
coldata$replicate.num <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 5);

####################################################################################################
### DESeq2 #########################################################################################
####################################################################################################
selected.samples <- coldata$FPE.num == 'FPE4';
dds <- DESeqDataSetFromMatrix(
  countData = cts[, selected.samples], 
  colData = coldata[selected.samples, ], 
  design = ~ treatment);
vsd <- vst(dds, blind = FALSE);

####################################################################################################
### PCA plot #######################################################################################
####################################################################################################
pcaData <- plotPCA(
  object = vsd, 
  intgroup = c(
    'FPE.num',
    'participant_id',
    'treatment',
    'treatment.time',
    'replicate.num',
    'source'), 
  returnData = TRUE);
percentVar <- round(100 * attr(pcaData, 'percentVar'));
ggplot(
  data = pcaData, 
  aes(x = PC1, y = PC2, shape = treatment, color = source)) +
  geom_point(size = 3) +
  xlab(paste0('PC1: ', percentVar[1], '% variance')) +
  ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
  coord_fixed() +
  scale_color_brewer(palette = 'Dark2') +
  ggtitle(label = 'Fibroblast Priming #4');

