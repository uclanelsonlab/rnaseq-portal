rm(list = ls());

library(DESeq2);
library(foreach);
library(doMC);
library(ggplot2);

registerDoMC(cores = 8);

####################################################################################################
### READ COUNT MATRIX ##############################################################################
####################################################################################################
path <- '../data/gene-counts_FPE/';
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
coldata <- read.csv(file = '../data/metadata.csv');
rownames(coldata) <- coldata$experiment_rna_short_read_id;
coldata <- coldata[ colnames(cts),];
coldata$FPE.num <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 1);
coldata$treatment <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 3);
coldata$TNFa.positive <- grepl(pattern = 'TNFa', x = coldata$treatment);
coldata$treatment.time <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 4);
coldata$replicate.num <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 5);

####################################################################################################
### DESeq2 #########################################################################################
####################################################################################################
dds <- DESeqDataSetFromMatrix(
  countData = cts, 
  colData = coldata, 
  design = ~ treatment + treatment.time);
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
    'TNFa.positive', 
    'treatment.time', 
    'replicate.num'), 
  returnData = TRUE);
percentVar <- round(100 * attr(pcaData, 'percentVar'));
ggplot(pcaData, aes(x = PC1, y = PC2, color = FPE.num, shape = TNFa.positive)) +
  geom_point(size = 3) +
  xlab(paste0('PC1: ', percentVar[1], '% variance')) +
  ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
  coord_fixed() +
  scale_color_brewer(palette = 'Set1')

