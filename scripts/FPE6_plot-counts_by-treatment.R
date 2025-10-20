library(BoutrosLab.utilities);

genes <- readRDS(file = 'data/2025-03-19_genes.rds');

####################################################################################################
### DESeq2 #########################################################################################
####################################################################################################
selected.samples <- coldata$FPE.num == 'FPE6';
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
  
  d$FPE.num <- sapply(strsplit(d$experiment_rna_short_read_id, "-"), `[`, 1);
  d$treatment <- sapply(strsplit(d$experiment_rna_short_read_id, "-"), `[`, 3);
  d$TNFa.positive <- grepl(pattern = 'TNFa', x = d$treatment);
  d$sub.treatment <- substr(
    x = d$treatment, 
    start = ifelse(test = startsWith(x = d$treatment, prefix = 'TNFa+'), yes = 6, no = 0), 
    stop = 100);
  d$sub.treatment[d$sub.treatment == 'none'] <- 'TNFa';
  d$treatment.time <- sapply(strsplit(d$experiment_rna_short_read_id, "-"), `[`, 4);
  d$replicate.num <- sapply(strsplit(d$experiment_rna_short_read_id, "-"), `[`, 5);
  d$sub.treatment <- factor(
    x = d$sub.treatment, 
    levels = c(
      'TNFa', 'rapa', 'TGFb+rapa', 'TGFb', 'TGFb+SB', 'SB', 'rapa+SB', 'KC7F2', 'TGFb+KC7F2',
      'EGFRi', 'EGFRi+EGF', 'EGF'));
  
  ggplot(d, aes(x = sub.treatment, y = count, color = TNFa.positive)) + 
    geom_boxplot(outliers = FALSE, show.legend = FALSE) +
    geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) + 
    ggtitle(label = gene.name) +
    scale_color_manual(values = c('#f1a340', '#998ec3')) +
    scale_y_log10(limits = c(1, NA));
}

####################################################################################################
### Plot counts by gene name #######################################################################
####################################################################################################
gene.name <- 'CXCL8';
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