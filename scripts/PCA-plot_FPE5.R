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
coldata$TGFb.positive <- grepl(pattern = 'TGFb', x = coldata$treatment);
coldata$co.treatment <- ifelse(
  test = coldata$TNFa.positive & coldata$TGFb.positive, 
  yes = 'TNFa+TGFb', 
  no = ifelse(test = coldata$TNFa.positive, yes = 'TNFa', no = ifelse(test = coldata$TGFb.positive, yes = 'TGFb', no = 'none')));
coldata$sub.treatment <- substr(x = coldata$treatment, start = ifelse(test = startsWith(x = coldata$treatment, prefix = 'TNFa+'), yes = 6, no = 0), stop = 100);
coldata$sub.treatment <- substr(x = coldata$sub.treatment, start = ifelse(test = startsWith(x = coldata$sub.treatment, prefix = 'TGFb+'), yes = 6, no = 0), stop = 100);
coldata$sub.treatment[coldata$sub.treatment == 'TNFa'] <- 'none';
coldata$sub.treatment[coldata$sub.treatment == 'TGFb'] <- 'none';
coldata$sub.treatment[coldata$sub.treatment == 'TNFa+TGFb'] <- 'none';
coldata$treatment.time <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 4);
coldata$replicate.num <- sapply(strsplit(coldata$experiment_rna_short_read_id, "-"), `[`, 5);

####################################################################################################
### DESeq2 #########################################################################################
####################################################################################################
selected.samples <- coldata$FPE.num == 'FPE5';
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
    'co.treatment',
    'none'), 
  returnData = TRUE);
percentVar <- round(100 * attr(pcaData, 'percentVar'));
ggplot(
  data = pcaData, 
  aes(x = PC1, y = PC2, color = sub.treatment, shape = co.treatment)) +
  geom_point(size = 3) +
  xlab(paste0('PC1: ', percentVar[1], '% variance')) +
  ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
  coord_fixed() +
  scale_color_brewer(palette = 'Dark2') +
  ggtitle(label = 'Fibroblast Priming #5') +
  guides(
    shape = guide_legend(title = "Co-treatment", order = 1),
    color = guide_legend(title = "Sub-treatment", order = 2)
  );

