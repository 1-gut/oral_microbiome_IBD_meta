# Package loading & setup -------------------------------------------------

# CRAN + Bioconductor packages required for the project
packages <- c(
  "stringr", "microbiome", "ggrepel", "rstatix", "phyloseq", "Biostrings", 
  "ggplot2", "tidyverse", "cowplot", "purrr", "pairwiseAdonis", "vegan", 
  "dplyr", "Maaslin2", "reshape2", "picante", "knitr", "tibble","RColorBrewer", 
  "ggside", "viridis", "scales", "ggsignif", "ggExtra", "ggpicrust2", "renv"
)

# Bioconductor-only packages
bioc_packages <- c("phyloseq", "Biostrings")

# Ensure BiocManager is available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install missing packages
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    
    if (pkg %in% bioc_packages) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

# Load all packages quietly
invisible(
  lapply(packages, function(pkg) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  })
)