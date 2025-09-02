library(BoutrosLab.utilities);

genes <- readRDS(file = 'data/2025-03-19_genes.rds');

####################################################################################################
### DESeq2 #########################################################################################
####################################################################################################
selected.samples <- coldata$FPE.num == 'FPE7';
dds <- DESeqDataSetFromMatrix(
  countData = cts[, selected.samples], 
  colData = coldata[selected.samples, ], 
  design = ~ treatment);

####################################################################################################
### DEA ############################################################################################
####################################################################################################
dds <- DESeq(dds);

####################################################################################################
### Plot counts ####################################################################################
####################################################################################################
plot.counts <- function(gene, gene.name){
  d <- plotCounts(dds, gene = gene, intgroup = 'treatment', returnData = TRUE);
  d$experiment_rna_short_read_id <- rownames(d);
  d <- merge(d, coldata);
  d$treatment <- factor(x = d$treatment, levels = c('none', 'iBET151', 'TNFa', 'TNFa+iBET151'));
  ggplot(d, aes(x = treatment, y = count, color = participant_id)) + 
    geom_boxplot(outliers = FALSE, color = 'black') +
    geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) + 
    ggtitle(label = gene.name) +
    scale_color_brewer(palette = 'Dark2') +
    scale_y_log10(limits = c(1, NA));
}

####################################################################################################
### Plot counts by gene name #######################################################################
####################################################################################################
gene.name <- 'DMD';
gene <- which(startsWith(
  x = rownames(cts), 
  prefix = genes$ensembl_gene_id[genes$external_gene_name == gene.name]));
plot.counts(gene, gene.name);

####################################################################################################
### Plot counts by ensembl gene id #################################################################
####################################################################################################
ensembl.gene.id <- 'ENSG00000109971';
gene.name <- genes$external_gene_name[genes$ensembl_gene_id == ensembl.gene.id];
gene <- which(startsWith(
  x = rownames(cts), 
  prefix = genes$ensembl_gene_id[genes$external_gene_name == gene.name]));
plot.counts(gene, gene.name);