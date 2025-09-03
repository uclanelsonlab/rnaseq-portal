# 🚀 Precomputed Data Generation Guide

This document explains how to generate precomputed data files to dramatically speed up the CCRD RNA-Seq Portal Shiny app.

## 📋 Overview

The Shiny app can use precomputed data files to skip expensive computations (DESeq2 analysis, VST transformation, PCA calculation) and load plots almost instantly. Without precomputation, each plot takes 10-30 seconds to generate. With precomputation, plots load in under 2 seconds.

## 🗂️ Precomputed Files Generated

### Core Data Files
- `data/cts.rds` - Combined count matrix (62,757 genes × 231 samples)
- `data/coldata_processed.rds` - Processed metadata with all derived columns
- `data/coldata_fpe6.rds` - FPE6-specific metadata processing

### Per-Experiment Files
- `data/fpe_all_vsd.rds` / `data/fpe_all_pca.rds` - All experiments VST and PCA
- `data/fpe4_vsd.rds` / `data/fpe4_pca.rds` - Experiment 4 VST and PCA
- `data/fpe5_vsd.rds` / `data/fpe5_pca.rds` - Experiment 5 VST and PCA
- `data/fpe6_vsd.rds` / `data/fpe6_pca.rds` - Experiment 6 VST and PCA
- `data/fpe7_vsd.rds` / `data/fpe7_pca.rds` - Experiment 7 VST and PCA

### FPE7 Optimization Files
- `data/fpe7_deseq_precomputed.rds` - Pre-analyzed DESeq2 object for FPE7
- `data/fpe7_deseq_info.rds` - Analysis metadata (creation date, sample count, etc.)
- `data/fpe7_norm_counts.rds` - Normalized counts matrix for gene expression plots

### Reference Files
- `data/gene_index_map.rds` - Gene symbol ↔ Ensembl ID mapping for fast lookups
- `data/legend_specs.rds` - Factor levels and legend specifications

## 🔧 Prerequisites

### Required R Packages
```r
install.packages(c("DESeq2", "ggplot2"))
# DESeq2 requires Bioconductor:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
```

### Required Data Files
Ensure these files exist in your `data/` directory:
- `data/gene-counts_FPE/` - Directory with 231 count files (`.gene_id.exon.ct.short.txt`)
- `data/metadata.csv` - Sample metadata file
- `data/2025-03-19_genes.rds` - Gene reference database (86,402 entries)

## 🚀 Generation Scripts

### 1. Generate Core Precomputed Data

**Script:** `scripts/precompute_common.R`

**What it does:**
- Combines 231 individual count files into a single matrix
- Processes metadata with all derived columns
- Generates VST and PCA data for all experiments
- Creates gene index mapping for fast lookups
- Builds legend specifications

**Run from project root:**
```bash
Rscript scripts/precompute_common.R
```

**Expected runtime:** ~3-5 minutes

**Output example:**
```
2025-09-02 22:10:06 | Loading count matrix from data/gene-counts_FPE ...
2025-09-02 22:10:13 | ✓ Count matrix loaded: 62757 genes x 231 samples
2025-09-02 22:10:13 | ✓ Saved cts.rds
2025-09-02 22:10:13 | ✓ Saved coldata_processed.rds
2025-09-02 22:10:13 | ✓ Saved coldata_fpe6.rds
2025-09-02 22:10:45 | ✓ Saved fpe_all_vsd.rds
2025-09-02 22:10:45 | ✓ Saved fpe_all_pca.rds
2025-09-02 22:10:46 | ✓ Saved fpe4_vsd.rds
2025-09-02 22:10:46 | ✓ Saved fpe4_pca.rds
2025-09-02 22:10:51 | ✓ Saved fpe5_vsd.rds
2025-09-02 22:10:51 | ✓ Saved fpe5_pca.rds
2025-09-02 22:11:07 | ✓ Saved fpe6_vsd.rds
2025-09-02 22:11:07 | ✓ Saved fpe6_pca.rds
2025-09-02 22:12:38 | ✓ Saved gene_index_map.rds
2025-09-02 22:12:38 | ✓ Saved legend_specs.rds
2025-09-02 22:12:38 | ✅ Precompute complete.
```

### 2. Generate FPE7 DESeq Analysis (Optional but Recommended)

**Script:** `scripts/precompute_fpe7_deseq.R`

**What it does:**
- Runs the computationally expensive DESeq2 analysis for FPE7 samples
- Saves the analyzed DESeq object for instant loading
- Generates normalized counts for gene expression plots

**Run from project root:**
```bash
Rscript scripts/precompute_fpe7_deseq.R
```

**Expected runtime:** ~4-6 minutes

**Output example:**
```
=== PRE-COMPUTING DESEQ ANALYSIS FOR FPE7 ===
Loading count matrix...
Loading metadata...
Filtering samples for FPE7...
Found 60 FPE7 samples
Creating DESeq2 dataset...
Running DESeq analysis (this may take several minutes)...
DESeq analysis completed in 4.23 minutes

=== PRE-COMPUTATION COMPLETE ===
Files created:
- DESeq object: ../data/fpe7_deseq_precomputed.rds
- Analysis info: ../data/fpe7_deseq_info.rds
- Analysis time: 4.23 minutes
- Samples analyzed: 60
- Genes analyzed: 62757

You can now use the fast-loading Shiny app!
```

## ⚡ Performance Impact

### Without Precomputation
- **App startup:** ~15-20 seconds (loading 231 count files)
- **Plot generation:** 10-30 seconds per plot
- **FPE7 plots:** 30-60 seconds (includes DESeq analysis)

### With Precomputation
- **App startup:** ~3-5 seconds (loading RDS files)
- **Plot generation:** <2 seconds per plot
- **FPE7 plots:** <2 seconds (uses precomputed DESeq)

### Storage Requirements
- **Raw data:** ~50MB (231 text files + metadata)
- **Precomputed data:** ~200MB (all RDS files)
- **Total speedup:** 10-15x faster plot generation

## 🔄 When to Regenerate

Regenerate precomputed data when:

### Always Regenerate
- **New data added:** New samples or experiments
- **Metadata changes:** Sample annotations updated
- **Gene reference updated:** New version of `2025-03-19_genes.rds`

### Selectively Regenerate
- **Code changes to modules:** Only if plot logic changes
- **UI changes:** No regeneration needed
- **Styling changes:** No regeneration needed

### Quick Regeneration Commands
```bash
# Regenerate everything
Rscript scripts/precompute_common.R
Rscript scripts/precompute_fpe7_deseq.R

# Regenerate only specific experiments (if needed)
# Edit precompute_common.R to comment out unneeded experiments
```

## 🛠️ Troubleshooting

### Common Issues

**Error: "No count files found"**
```
Solution: Ensure data/gene-counts_FPE/ directory exists with .txt files
Check: ls data/gene-counts_FPE/ | wc -l  # Should show 231
```

**Error: "Genes reference file not found"**
```
Solution: Ensure data/2025-03-19_genes.rds exists
Impact: Gene index mapping will be skipped (FPE7 gene search limited)
```

**Error: "some variables in design formula are characters"**
```
Status: This is just a warning, not an error
Impact: DESeq2 automatically converts to factors
```

**Memory issues during generation**
```
Solution: Increase R memory limit or run on a machine with more RAM
Alternative: Generate experiments individually by editing the script
```

### Validation

**Check precomputed files exist:**
```bash
ls -la data/*.rds
# Should show all the files listed in "Precomputed Files Generated"
```

**Test app with precomputed data:**
```bash
cd main/
R -e "shiny::runApp('app.R', port=3838)"
# Should start in ~3-5 seconds and plots should load quickly
```

**Verify FPE7 optimization:**
```
Look for this message in app startup:
"✓ Pre-computed FPE7 DESeq loaded (60 samples, 62757 genes)"
```

## 📁 File Structure After Generation

```
rnaseq-portal/
├── data/
│   ├── cts.rds                          # Combined count matrix
│   ├── coldata_processed.rds             # Processed metadata
│   ├── coldata_fpe6.rds                  # FPE6-specific metadata
│   ├── fpe_all_vsd.rds / fpe_all_pca.rds # All experiments
│   ├── fpe4_vsd.rds / fpe4_pca.rds       # Experiment 4
│   ├── fpe5_vsd.rds / fpe5_pca.rds       # Experiment 5
│   ├── fpe6_vsd.rds / fpe6_pca.rds       # Experiment 6
│   ├── fpe7_vsd.rds / fpe7_pca.rds       # Experiment 7
│   ├── fpe7_deseq_precomputed.rds        # FPE7 DESeq object
│   ├── fpe7_deseq_info.rds               # FPE7 analysis metadata
│   ├── fpe7_norm_counts.rds              # FPE7 normalized counts
│   ├── gene_index_map.rds                # Gene mapping
│   └── legend_specs.rds                  # Legend specifications
├── scripts/
│   ├── precompute_common.R               # Main precomputation script
│   └── precompute_fpe7_deseq.R          # FPE7 DESeq precomputation
└── main/
    └── app.R                            # Shiny app (uses precomputed data)
```

## 🎯 Best Practices

1. **Run precomputation on a powerful machine** - The initial generation is CPU-intensive
2. **Store precomputed files in version control** - Consider using Git LFS for large RDS files
3. **Document when you regenerate** - Keep track of data versions and regeneration dates
4. **Test after regeneration** - Always verify the app works correctly with new precomputed data
5. **Monitor file sizes** - Precomputed files should be consistent in size between runs

---

**Need help?** Check the main `README.md` for general app documentation or open an issue if you encounter problems with precomputation.
