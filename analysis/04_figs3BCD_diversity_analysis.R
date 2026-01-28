## ===============================
## analysis/04_fig3BCD_diversity_analysis.R
## Perform diversity analysis and plot

source("analysis/00_config.R")
source("R/packages.R")
source("R/functions_phyloseq.R")

# Load data
metadata <- read.table(FILE_PS_META, sep = "\t", header = TRUE, row.names = 1)
rel <- readRDS(FILE_PS_PROCESSED_REL)
rel_sgn <- readRDS(FILE_PS_PROCESSED_REL_SGN)

# Figure 3B - alpha boxplot  -----------------------------------------------------
# Alpha
alpha <- microbiome::alpha(rel, index = "all")
merge <- merge(metadata, alpha, by = "row.names")

df_long <- merge %>%
  pivot_longer(cols = c(study_group_name, IBD_group_name),
               names_to = "source", values_to = "group") %>%
  filter(!(source == "study_group_name" & group == "IBD"))

df_long <- df_long %>%
  mutate(group = factor(group, levels = c("HC","IBD","CD","UC")))

df_long <- df_long %>%
  mutate(group = recode(group,
                        "HC" = "HC"))

# Define colour palette
colours <- c(
  "HC" = "#648fff",
  "IBD" = "#785ef0",
  "CD" = "#dc267f",
  "UC" = "#ffb000"
)

# Plot the data
boxplot <- ggplot(df_long, aes(x = group, y = observed)) +
  geom_boxplot(aes(fill = group), alpha = 0.9) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +  # Increase y limit for annotation space
  labs(x = NULL, y = "Richness") +
  geom_signif(
    comparisons = list(
      c("HC", "IBD"),
      c("HC", "CD"),
      c("CD", "UC"),
      c("HC", "UC")
    ),
    annotations = c("***", "***", "***", "ns"),
    y_position = c(65, 68, 71, 74),  # Raise the final "ns" label
    tip_length = 0.01,
    textsize = 6,
    vjust = 0.5
  ) +
  theme_cowplot(20) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Save plot
ggsave(file.path(DIR_FIG_MAIN, "Fig3B_alpha.pdf"), plot = boxplot, width = 3, height = 5, dpi = 300)

# Figures 3C + D - beta all and SGN ------------------------------------------------
# Run beta function IBD
beta_all_plot <- plot_beta_diversity(rel, comparator = "IBD_group_name")

# Store plot
beta_all_plot <- beta_all_plot$plot

# Save
ggsave(file.path(DIR_FIG_MAIN, "Fig3C_beta_IBD.png"), plot = beta_all_plot, width = 7, height = 7, dpi = 300)

# Run beta function SGN
beta_sgn_plot <- plot_beta_diversity(sgn_rel, comparator = "study_group_name")

# Store plot
beta_sgn_plot <- beta_sgn_plot$plot

#Save
ggsave(file.path(DIR_FIG_MAIN, "fig3D_sgn_ad.png"), plot = beta_sgn_plot, width = 7, height = 7, dpi = 300)