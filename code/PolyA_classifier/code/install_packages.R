# List of CRAN packages
cran_packages <- c("dplyr", "caTools", "randomForest", "stringr", 
                   "ranger", "tidyr", "ggplot2", "readr", "patchwork", "gridExtra", "plotROC")

# List of Bioconductor packages
bioconductor_packages <- c("GenomicRanges", "GenomicFeatures", 
                           "BSgenome.Hsapiens.UCSC.hg38", "Biostrings", 
                           "rtracklayer", "memes")

# Function to check and install CRAN packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Install missing CRAN packages
invisible(sapply(cran_packages, install_if_missing))

# Install Bioconductor Manager if not already installed
if (!require("BiocManager", character.only = TRUE)) {
  install.packages("BiocManager")
}

# Function to check and install Bioconductor packages
install_bioc_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Install missing Bioconductor packages
invisible(sapply(bioconductor_packages, install_bioc_if_missing))

# Install devtools if not already installed (needed for GitHub packages)
if (!require("devtools", character.only = TRUE)) {
  install.packages("devtools", dependencies = TRUE)
}

# Install LATER package from GitHub
if (!require("LATER", character.only = TRUE)) {
  devtools::install_github("hilgers-lab/LATER")
}
