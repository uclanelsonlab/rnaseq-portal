# CCRD RNA-Seq Portal Docker Image
# Base image with R 4.3
FROM rocker/r-ver:4.3

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libicu-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN R -e "install.packages(c('shiny', 'ggplot2', 'foreach', 'doMC', 'gridExtra'), repos='https://cloud.r-project.org/')"

# Install BiocManager and DESeq2
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')" && \
    R -e "BiocManager::install('DESeq2', ask=FALSE, update=FALSE)"

# Create app directory
RUN mkdir -p /app

# Set working directory
WORKDIR /app

# Copy application files
COPY main/ /app/main/
COPY scripts/ /app/scripts/

# The data directory will be mounted as a volume at runtime
# This keeps the image size small and allows easy data updates

# Expose port for Shiny
EXPOSE 3838

# Set environment variables
ENV SHINY_LOG_LEVEL=INFO

# Run the Shiny app
CMD ["R", "-e", "shiny::runApp('/app/main/app.R', host='0.0.0.0', port=3838)"]

