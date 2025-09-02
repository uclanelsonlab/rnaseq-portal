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
coldata$treatment <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 3);
coldata$none <- grepl(pattern = 'none', x = coldata$treatment);
coldata$TNFa.positive <- grepl(pattern = 'TNFa', x = coldata$treatment);
coldata$sub.treatment <- substr(x = coldata$treatment, start = ifelse(test = startsWith(x = coldata$treatment, prefix = 'TNFa+'), yes = 6, no = 0), stop = 100);
coldata$sub.treatment[coldata$sub.treatment == 'none'] <- 'TNFa';
coldata$treatment.time <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 4);
coldata$replicate.num <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 5);

####################################################################################################
### DESeq2 #########################################################################################
####################################################################################################
selected.samples <- coldata$FPE.num == 'FPE6';
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
    'sub.treatment',
    'TNFa.positive',
    'none'), 
  returnData = TRUE);
percentVar <- round(100 * attr(pcaData, 'percentVar'));
ggplot(
  data = pcaData, 
  aes(x = PC1, y = PC2, color = sub.treatment, shape = TNFa.positive)) +
  geom_point(size = 3) +
  geom_point(data = pcaData[pcaData$none, ], size = 3, color = 'black') +
  xlab(paste0('PC1: ', percentVar[1], '% variance')) +
  ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
  coord_fixed() +
  scale_color_brewer(palette = 'Set3') +
  ggtitle(label = 'Fibroblast Priming #6');

