# CCRD RNA-Seq Portal 🧬

An interactive Shiny web application for exploring PCA plots from RNA-seq analysis across different fibroblast priming experiments.

## 🚀 **Quick Start**

```r
# Navigate to the main directory
cd /path/to/CCRD-RNAseq-portal/main/

# Launch the Shiny app
R -e "shiny::runApp('app.R')"
```

## 📊 **Features**

- **Interactive Plot Selection**: Switch between different experiment views with a dropdown
- **All Experiments Overview**: Shows all samples colored by experiment number  
- **Individual Experiments**: Dedicated visualizations for Experiments 4, 5, 6, and 7
- **Interactive Gene Selection**: Enter custom gene names for FPE7 gene expression plots
- **Modular Code Structure**: Easy customization of individual experiment plots
- **Real-time Updates**: Plots change instantly when selections are made
- **Intelligent Gene Search**: Automatic gene matching with fallback to most expressed genes

## 📋 **Required R Packages**

### CRAN packages:
```r
install.packages(c("shiny", "ggplot2", "foreach", "doMC", "gridExtra"))
```

### Bioconductor packages:
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
```

## 📁 **Project Structure**

```
CCRD-RNAseq-portal/
├── README.md                   # This comprehensive guide
├── data/
│   ├── gene-counts_FPE/       # Count matrix files (231 files)
│   └── metadata.csv           # Sample metadata
├── main/
│   ├── app.R                  # Main Shiny application
│   └── modules/               # Modular plot functions
│       ├── plot_all.R         # All experiments overview
│       ├── plot_fpe4.R        # Experiment 4 plot
│       ├── plot_fpe5.R        # Experiment 5 plot
│       ├── plot_fpe6.R        # Experiment 6 plot
│       └── plot_fpe7.R        # Experiment 7 plot
└── scripts/                   # Original individual plotting scripts
    ├── PCA-plot_all.R
    ├── PCA-plot_FPE4.R
    ├── PCA-plot_FPE5.R
    ├── PCA-plot_FPE6.R
    └── PCA-plot_FPE7.R
```

## 🎯 **Plot Types**

| Experiment | Visualization | Color | Shape | Special Features |
|------------|---------------|-------|-------|------------------|
| **All** | Overview of all samples | FPE number | TNFα treatment | Complete dataset view |
| **Experiment 4** | FPE4 samples only | Source | Treatment | Focus on treatment effects |
| **Experiment 5** | FPE5 samples only | Sub-treatment | Co-treatment | Black highlighting for 'none' |
| **Experiment 6** | FPE6 samples only | Sub-treatment | TNFα status | Black highlighting for 'none' |  
| **Experiment 7** | FPE7 PCA + Gene Count | Participant ID | Treatment | Combined: PCA plot + Interactive gene expression |

## 🧬 **Interactive Gene Selection (Experiment 7)**

Experiment 7 includes a unique **interactive gene selection** feature:

### 🎯 **How to Use**
1. Select "**Experiment 7**" from the dropdown
2. Enter a **gene symbol** (e.g., `DMD`, `ACTA1`) **OR** an **Ensembl ID** (e.g., `ENSG00000198947`) in the **Gene Name** text box
3. Press Enter or click elsewhere to update the plot
4. The gene expression plot updates instantly with your selected gene

### 🔍 **Gene Search Method**
The app supports **both gene symbols and Ensembl IDs** using a comprehensive genes reference database (86,402 entries):

#### **For Gene Symbols (e.g., `DMD`, `TTR`)**
1. **Direct Lookup**: Uses `genes$external_gene_name` to find the corresponding Ensembl ID
2. **Count Matrix Search**: Searches count matrix row names using `startsWith()` with the Ensembl ID prefix
3. **Exact Matching**: Finds genes where row names start with the correct Ensembl ID

#### **For Ensembl IDs (e.g., `ENSG00000198947`)**
1. **ID Detection**: Automatically detects Ensembl IDs (starting with "ENSG")
2. **Gene Name Lookup**: Uses `genes$ensembl_gene_id` to find the corresponding gene symbol
3. **Count Matrix Search**: Proceeds with normal search using the mapped gene symbol

#### **Fallback Strategy**
If gene not found in reference, shows the most highly expressed gene with a warning

### ✨ **Examples**

#### **Gene Symbols** (any gene from the reference database):
- **`DMD`**: Duchenne muscular dystrophy gene (default)
- **`TTR`**: Transthyretin 
- **`ACTA1`**: Skeletal muscle alpha-actin
- **`MYH7`**: Cardiac/skeletal muscle myosin
- **`GAPDH`**: Housekeeping gene control

#### **Ensembl IDs** (automatically detected):
- **`ENSG00000198947`**: DMD (Duchenne muscular dystrophy)
- **`ENSG00000132664`**: TTR (Transthyretin)
- **`ENSG00000143632`**: ACTA1 (Skeletal muscle alpha-actin)
- **`ENSG00000092054`**: MYH7 (Cardiac/skeletal muscle myosin)

## ✏️ **Customizing Individual Experiments**

### 🔧 **Easy Customization Workflow**

1. **Open the specific module**: `modules/plot_fpe4.R`
2. **Make your changes** (see examples below)
3. **Save the file**
4. **Restart the app**: `shiny::runApp('app.R')`

### 🎨 **Common Customizations**

#### Change Colors
```r
# In any module, find:
scale_color_brewer(palette = 'Dark2')

# Replace with custom colors:
scale_color_manual(values = c("red", "blue", "green", "purple"))
```

#### Update Plot Title
```r
# Find:
ggtitle(label = 'Fibroblast Priming #4')

# Replace with:
ggtitle(label = 'My Custom Experiment 4 Title')
```

#### Add Plot Elements
```r
# Add sample labels:
+ geom_text(aes(label = sample_id), size = 2, hjust = 1.1)

# Add trend lines:
+ geom_smooth(method = "lm", se = FALSE, alpha = 0.5)

# Modify point size:
geom_point(size = 4)  # Change from size = 3
```

#### Change Color Palettes
```r
# Available Brewer palettes:
scale_color_brewer(palette = 'Set1')      # Bright colors
scale_color_brewer(palette = 'Dark2')     # Darker colors  
scale_color_brewer(palette = 'Set3')      # Pastel colors
scale_color_brewer(palette = 'Spectral')  # Rainbow colors
```

## 🏗️ **Modular Architecture Benefits**

✅ **Independent Editing**: Modify one experiment without affecting others  
✅ **Error Isolation**: Problems in one plot won't break the entire app  
✅ **Team Collaboration**: Multiple people can work on different experiments  
✅ **Clean Code**: Each function is self-contained and well-documented  
✅ **Easy Testing**: Test individual modules independently  

## 🔧 **Development**

### Adding a New Experiment

1. **Create new module**: `modules/plot_fpe8.R`
```r
generate_fpe8_plot <- function(cts, coldata) {
  tryCatch({
    # Filter for your experiment
    selected.samples <- coldata$FPE.num == 'FPE8'
    
    # Check data availability
    if(sum(selected.samples) == 0) {
      stop("No FPE8 samples found")
    }
    
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(
      countData = cts[, selected.samples], 
      colData = coldata[selected.samples, ], 
      design = ~ treatment)
    vsd <- vst(dds, blind = FALSE)
    
    # Generate PCA
    pcaData <- plotPCA(object = vsd, intgroup = c('treatment'), returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, 'percentVar'))
    
    # Create plot
    ggplot(pcaData, aes(x = PC1, y = PC2, color = treatment)) +
      geom_point(size = 3) +
      xlab(paste0('PC1: ', percentVar[1], '% variance')) +
      ylab(paste0('PC2: ', percentVar[2], '% variance')) + 
      coord_fixed() +
      ggtitle('Fibroblast Priming #8') +
      theme_minimal()
      
  }, error = function(e) {
    # Error handling
    ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
      geom_text(label = paste("Error in FPE8 plot:", e$message), color = "red") +
      theme_void()
  })
}
```

2. **Update main app.R**:
   - Add `source("modules/plot_fpe8.R")` in the module loading section
   - Add `"Experiment 8" = "fpe8"` to the UI dropdown choices
   - Add `"fpe8" = generate_fpe8_plot(cts, coldata)` to the server switch statement

### Testing Individual Modules
```r
# Test a specific module
source("modules/plot_fpe4.R")

# Test with your data (after loading in R console)
plot <- generate_fpe4_plot(cts, coldata)
print(plot)
```

## 🔍 **Troubleshooting**

### App Won't Start
- ✅ Check all R packages are installed (`DESeq2`, `shiny`, `ggplot2`, etc.)
- ✅ Verify data files exist: `../data/metadata.csv` and `../data/gene-counts_FPE/`
- ✅ Ensure you're in the `main/` directory when running the app
- ✅ Check R console for specific error messages

### Plots Not Displaying
- ✅ Look for error messages in R console
- ✅ Verify the experiment has samples (check FPE numbers in metadata)
- ✅ Test individual modules with `source("modules/plot_name.R")`

### Module Errors
- ✅ Check function syntax (missing commas, unmatched parentheses)
- ✅ Ensure required data columns exist in your metadata
- ✅ Verify all `aes()` mappings reference valid column names

### Data Loading Issues
- ✅ Confirm file paths are correct (`../data/` relative to `main/` directory)
- ✅ Check file permissions and readability
- ✅ Verify count files are in expected format

## 📊 **Data Requirements**

### Count Files (`data/gene-counts_FPE/`)
- Format: Tab-separated text files
- Columns: `gene_id` and count values
- Naming: Must end with `.gene_id.exon.ct.short.txt`

### Metadata (`data/metadata.csv`)
- Required columns: `experiment_rna_short_read_id`, `participant_id`, `affected_status`
- Sample IDs should match count file names (after removing suffix)
- FPE numbers extracted from sample IDs (format: `FPE4-...`, `FPE5-...`, etc.)

## 🎊 **Advanced Customization**

### Custom Color Scales
```r
# Define custom colors for specific values
custom_colors <- c("control" = "#1f77b4", "treated" = "#ff7f0e", "combo" = "#2ca02c")
scale_color_manual(values = custom_colors)
```

### Interactive Elements
```r
# Add hover information (requires plotly)
library(plotly)
ggplotly(your_plot, tooltip = c("x", "y", "colour"))
```

### Export Options
```r
# Save plots programmatically
ggsave("experiment4_plot.png", plot, width = 8, height = 6, dpi = 300)
```

## 🏃‍♂️ **Performance Notes**

- **Data Loading**: Count matrix (62,757 genes × 231 samples) loads once at startup
- **Plot Generation**: PCA calculations performed on-demand when switching experiments  
- **FPE7 Optimization**: Pre-computed DESeq analysis loaded at startup for instant plot generation
- **Memory Usage**: ~250MB RAM for typical dataset size (includes pre-computed results)
- **Load Time**: ~30-60 seconds initial startup, <2 seconds plot switching (FPE7 now instant!)

### 🚀 **FPE7 Speed Optimization**

Experiment 7 uses pre-computed DESeq2 results for lightning-fast plot generation:

- **Pre-computed**: DESeq analysis runs once offline (~4 minutes) 
- **Runtime**: Plot generation now takes <2 seconds instead of 3-4 minutes
- **Auto-loading**: Pre-computed results loaded automatically at app startup

### 🔄 **Regenerating Pre-computed Results**

If you update your data, regenerate the pre-computed DESeq analysis:

```bash
# Navigate to scripts directory
cd scripts/

# Run pre-computation (takes ~4 minutes)
Rscript precompute_fpe7_deseq.R

# Restart your Shiny app to use new results
```

**When to regenerate:**
- After adding/removing FPE7 samples
- After updating count data or metadata
- When changing DESeq2 analysis parameters

## 📞 **Support**

For questions or issues:
1. Check the troubleshooting section above
2. Look for error messages in the R console
3. Test individual modules independently
4. Verify data file formats and paths

---

**🧬 Happy RNA-seq analysis! This modular structure makes your portal highly customizable and maintainable.** ✨
