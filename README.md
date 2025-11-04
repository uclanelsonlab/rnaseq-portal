# CCRD RNA-Seq Portal 🧬

An interactive Shiny web application for exploring PCA plots from RNA-seq analysis across different fibroblast priming experiments.

## 🚀 **Quick Start**

### Option 1: Docker (Recommended)

The easiest way to run the app is using Docker:

```bash
# 1. Download and extract the data directory
aws s3 cp s3://ucla-rare-diseases/UCLA-UDN/gcarvalho_test/rnaseq/portal/rnaseq-portal-data.zip .
unzip rnaseq-portal-data.zip

# 2. Build and run with Docker Compose
docker-compose up -d

# 3. Access the app at http://localhost:3838
```

Or using Docker directly:

```bash
# Build the image
docker build -t rnaseq-portal .

# Run the container
docker run -d -p 3838:3838 -v $(pwd)/data:/app/data:ro --name rnaseq-portal rnaseq-portal

# Access the app at http://localhost:3838
```

### Option 2: Local R Installation

```bash
# 1. Download and extract the data directory
aws s3 cp s3://ucla-rare-diseases/UCLA-UDN/gcarvalho_test/rnaseq/portal/rnaseq-portal-data.zip .
unzip rnaseq-portal-data.zip

# 2. Navigate to the main directory
cd /path/to/CCRD-RNAseq-portal/main/

# 3. Launch the Shiny app
R -e "shiny::runApp('app.R')"
```

## 📊 **Features**

- **Interactive Plot Selection**: Switch between different experiment views with a dropdown
- **All Experiments Overview**: Shows all samples colored by experiment number  
- **Individual Experiments**: Dedicated visualizations for Experiments 4, 5, 6, and 7
- **Interactive Gene Selection**: Enter custom gene names for FPE6 and FPE7 gene expression plots
- **Modular Code Structure**: Easy customization of individual experiment plots
- **Real-time Updates**: Plots change instantly when selections are made
- **Intelligent Gene Search**: Automatic gene matching with fallback to most expressed genes
- **Grouped Box Plots**: FPE6 gene counts displayed with separate box plots for TNFα+ and TNFα- samples

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
├── README.md                              # This comprehensive guide
├── DOCKER.md                              # Docker deployment guide
├── Dockerfile                             # Docker image definition
├── docker-compose.yml                     # Docker Compose configuration
├── .dockerignore                          # Docker build exclusions
├── data/                                  # Data directory (download from S3)
│   ├── gene-counts_FPE/                  # Count matrix files (231 files)
│   ├── metadata.csv                      # Sample metadata
│   ├── fpe6_deseq_precomputed.rds        # Pre-computed FPE6 DESeq analysis
│   └── fpe7_deseq_precomputed.rds        # Pre-computed FPE7 DESeq analysis
├── main/
│   ├── app.R                             # Main Shiny application
│   └── modules/                          # Modular plot functions
│       ├── plot_all.R                    # All experiments overview
│       ├── plot_fpe4.R                   # Experiment 4 plot
│       ├── plot_fpe5.R                   # Experiment 5 plot
│       ├── plot_fpe6.R                   # Experiment 6 plot (with gene counts)
│       └── plot_fpe7.R                   # Experiment 7 plot (with gene counts)
└── scripts/                              # Plotting and pre-computation scripts
    ├── PCA-plot_all.R
    ├── PCA-plot_FPE4.R
    ├── PCA-plot_FPE5.R
    ├── PCA-plot_FPE6.R
    ├── PCA-plot_FPE7.R
    ├── FPE6_plot-counts_by-treatment.R   # FPE6 gene count visualization
    ├── FPE7_plot-counts_by-treatment.R   # FPE7 gene count visualization
    ├── precompute_fpe6_deseq.R           # Pre-compute FPE6 DESeq analysis
    └── precompute_fpe7_deseq.R           # Pre-compute FPE7 DESeq analysis
```

## 🎯 **Plot Types**

| Experiment | Visualization | Color | Shape | Special Features |
|------------|---------------|-------|-------|------------------|
| **All** | Overview of all samples | FPE number | TNFα treatment | Complete dataset view |
| **Experiment 4** | FPE4 samples only | Source | Treatment | Focus on treatment effects |
| **Experiment 5** | FPE5 samples only | Sub-treatment | Co-treatment | Black highlighting for 'none' |
| **Experiment 6** | FPE6 PCA + Gene Count | Sub-treatment | TNFα status | Combined: PCA plot + Interactive gene expression with grouped box plots |  
| **Experiment 7** | FPE7 PCA + Gene Count | Participant ID | Treatment | Combined: PCA plot + Interactive gene expression |

## 🧬 **Interactive Gene Selection**

### **Experiment 6 - Grouped Box Plots**

Experiment 6 features **interactive gene selection with grouped box plots**:

#### 🎯 **How to Use**
1. Select "**Experiment 6**" from the dropdown
2. Enter a **gene symbol** (e.g., `CXCL8`, `DMD`) **OR** an **Ensembl ID** (e.g., `ENSG00000169429`) in the **Gene Name** text box
3. Press Enter or click elsewhere to update the plot
4. The gene expression plot updates instantly with your selected gene

#### 📊 **Visualization Features**
- **Dual Box Plots**: Shows **two box plots per sub-treatment** - one for TNFα+ samples and one for TNFα- samples
- **Color Coding**: Orange (#f1a340) for TNFα+ and Purple (#998ec3) for TNFα-
- **Individual Points**: Overlaid jittered points show individual sample values
- **Log Scale**: Y-axis uses log10 scale for better visualization of expression ranges
- **12 Sub-treatments**: TNFa, Rapamycin, TGFb+Rapa, TGFb, TGFb+SB, SB, Rapa+SB, KC7F2, TGFb+KC7F2, EGFRi, EGFRi+EGF, EGF

#### 🔍 **Gene Search Method**
The app supports **both gene symbols and Ensembl IDs** using a comprehensive genes reference database:

- **Gene Symbols**: Direct lookup using `genes$external_gene_name` to find corresponding Ensembl ID
- **Ensembl IDs**: Automatically detected (starting with "ENSG"), mapped to gene symbol
- **Fallback Strategy**: If gene not found, shows the most highly expressed gene with a warning

#### ✨ **Examples**

**Gene Symbols** (inflammation and fibrosis markers):
- **`CXCL8`**: Interleukin-8 (default) - inflammatory chemokine
- **`IL6`**: Interleukin-6 - pro-inflammatory cytokine
- **`COL1A1`**: Collagen Type I Alpha 1 - fibrosis marker
- **`ACTA2`**: Alpha-smooth muscle actin - myofibroblast marker
- **`TGFB1`**: Transforming growth factor beta 1

**Ensembl IDs**:
- **`ENSG00000169429`**: CXCL8 (Interleukin-8)
- **`ENSG00000136244`**: IL6 (Interleukin-6)
- **`ENSG00000108821`**: COL1A1 (Collagen Type I Alpha 1)

### **Experiment 7 - Time Series**

Experiment 7 includes **interactive gene selection with time series analysis**:

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

## 🐳 **Docker Deployment**

For detailed Docker documentation, troubleshooting, and production deployment guide, see **[DOCKER.md](DOCKER.md)**.

### Building the Image

```bash
# Build with docker-compose
docker-compose build

# Or build directly
docker build -t rnaseq-portal .
```

### Running the Container

```bash
# Start with docker-compose (recommended)
docker-compose up -d

# Or run directly
docker run -d -p 3838:3838 -v $(pwd)/data:/app/data:ro --name rnaseq-portal rnaseq-portal
```

### Managing the Container

```bash
# View logs
docker-compose logs -f
# Or: docker logs -f rnaseq-portal

# Stop the container
docker-compose down
# Or: docker stop rnaseq-portal

# Restart the container
docker-compose restart
# Or: docker restart rnaseq-portal

# Remove the container
docker-compose down -v
# Or: docker rm -f rnaseq-portal
```

### Accessing the App

Once running, access the app at: **http://localhost:3838**

### Docker Image Details

- **Base Image**: `rocker/r-ver:4.3`
- **Installed Packages**: shiny, ggplot2, foreach, doMC, gridExtra, DESeq2
- **Port**: 3838
- **Data Volume**: Mounted at `/app/data` (read-only)
- **Health Check**: Enabled with 30s interval

### Updating Data

To update the data without rebuilding:

```bash
# Stop the container
docker-compose down

# Update your data directory
aws s3 cp s3://ucla-rare-diseases/UCLA-UDN/gcarvalho_test/rnaseq/portal/rnaseq-portal-data.zip .
unzip -o rnaseq-portal-data.zip

# Restart the container
docker-compose up -d
```

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
      geom_text(label = paste("Error in FPE8 plot:", e$message), color = "red")
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

### Docker Issues

#### Container Won't Start
- ✅ Check Docker is running: `docker ps`
- ✅ Verify data directory exists and is in the right location
- ✅ Check logs: `docker-compose logs` or `docker logs rnaseq-portal`
- ✅ Ensure port 3838 is not already in use: `lsof -i :3838` (macOS/Linux)

#### Can't Access the App
- ✅ Verify container is running: `docker ps`
- ✅ Check health status: `docker inspect rnaseq-portal | grep -A 5 Health`
- ✅ Confirm port mapping: should see `0.0.0.0:3838->3838/tcp`
- ✅ Try accessing: http://localhost:3838

#### Data Not Loading
- ✅ Verify data volume is mounted: `docker inspect rnaseq-portal | grep -A 10 Mounts`
- ✅ Check data directory contains expected files
- ✅ Ensure data was downloaded and extracted correctly
- ✅ Verify file permissions (data should be readable)

### Local R Installation Issues

#### App Won't Start
- ✅ Check all R packages are installed (`DESeq2`, `shiny`, `ggplot2`, etc.)
- ✅ Verify data files exist: `../data/metadata.csv` and `../data/gene-counts_FPE/`
- ✅ Ensure you're in the `main/` directory when running the app
- ✅ Check R console for specific error messages

#### Plots Not Displaying
- ✅ Look for error messages in R console
- ✅ Verify the experiment has samples (check FPE numbers in metadata)
- ✅ Test individual modules with `source("modules/plot_name.R")`

#### Module Errors
- ✅ Check function syntax (missing commas, unmatched parentheses)
- ✅ Ensure required data columns exist in your metadata
- ✅ Verify all `aes()` mappings reference valid column names

#### Data Loading Issues
- ✅ Confirm file paths are correct (`../data/` relative to `main/` directory)
- ✅ Check file permissions and readability
- ✅ Verify count files are in expected format

## 📊 **Data Requirements**

### Downloading the Data

The complete `data/` directory is available from S3:

```bash
aws s3 cp s3://ucla-rare-diseases/UCLA-UDN/gcarvalho_test/rnaseq/portal/rnaseq-portal-data.zip .
unzip rnaseq-portal-data.zip
```

This archive includes all necessary files: count matrices, metadata, pre-computed DESeq analyses, and gene reference databases.

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
- **FPE6 & FPE7 Optimization**: Pre-computed DESeq analysis loaded at startup for instant plot generation
- **Memory Usage**: ~250MB RAM for typical dataset size (includes pre-computed results)
- **Load Time**: ~30-60 seconds initial startup, <2 seconds plot switching (FPE6 & FPE7 now instant!)

### 🚀 **Speed Optimization (FPE6 & FPE7)**

Experiments 6 and 7 use pre-computed DESeq2 results for lightning-fast plot generation:

- **Pre-computed**: DESeq analysis runs once offline (~4 minutes per experiment) 
- **Runtime**: Plot generation now takes <2 seconds instead of 3-4 minutes
- **Auto-loading**: Pre-computed results loaded automatically at app startup
- **Grouped Visualization**: FPE6 displays dual box plots (TNFα+ vs TNFα-) for each sub-treatment

### 🔄 **Regenerating Pre-computed Results**

If you update your data, regenerate the pre-computed DESeq analysis:

```bash
# Navigate to scripts directory
cd scripts/

# Run pre-computation for FPE6 (takes ~4 minutes)
Rscript precompute_fpe6_deseq.R

# Run pre-computation for FPE7 (takes ~4 minutes)
Rscript precompute_fpe7_deseq.R

# Restart your Shiny app to use new results
```

**When to regenerate:**
- After adding/removing FPE6 or FPE7 samples
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
