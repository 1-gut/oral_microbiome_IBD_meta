## ===============================
## analysis/06_figS7_beta_continent.R
## Diversity analysis

source("analysis/00_config.R")
source("R/packages.R")
source("R/functions_phyloseq.R")

# Load data
metadata <- read.table(FILE_PS_META, sep = "\t", header = TRUE, row.names = 1)
rel <- readRDS(FILE_PS_PROCESSED_REL)
rel_sgn <- readRDS(FILE_PS_PROCESSED_REL_SGN)

# Figure S7 - beta continent ----------------------------------------------

# Beta all
con_beta_all_plot <- plot_beta_diversity(rel, comparator = "continent")

# SGN
con_beta_sgn_plot <- plot_beta_diversity_faceted(
  physeq_object = rel_sgn,
  comparator = "continent",
  facet_variable = "study_group_name"
)

# Run the function with default p-value adjustment
con_sgn_adonis_results <- adonis_with_subsets(
  physeq_object = rel_sgn,
  comparator = "continent"
)

con_beta_all_plot <- con_beta_all_plot$plot
con_beta_sgn_plot <- con_beta_sgn_plot$plot

ggsave(file.path(DIR_FIG_SUPP, "figS7A_beta_continent.png"), plot = con_beta_all_plot, width = 8.2, height = 8, dpi = 300)
ggsave(file.path(DIR_FIG_SUPP, "figS7B_beta_continent_SGN.png"), plot = con_beta_sgn_plot, width = 14, height = 7.5, dpi = 300)
