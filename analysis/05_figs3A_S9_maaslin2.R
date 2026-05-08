## ===============================
## analysis/05_fig4B_maaslin2.R
## Picrust2 analysis

source("analysis/00_config.R")
source("R/packages.R")
source("R/functions_phyloseq.R")

# Load data
metadata <- read.table(FILE_PS_META, sep = "\t", header = TRUE, row.names = 1)
rel <- readRDS(FILE_PS_PROCESSED_REL)
rel_sgn <- readRDS(FILE_PS_PROCESSED_REL_SGN)

# Figures 4B - Maaslin2 -----------------------------------------------------

# IBD vs HC
ps.ob <- rel

metadata_df <- as.data.frame(sample_data(ps.ob))
metadata_df <- as.matrix(metadata_df)
metadata_df <- as.data.frame(metadata_df)

# Extract OTU table
otu <- otu_table(ps.ob)
feat <- as.data.frame(otu)
feat <- t(feat)

## Tax table
tax <- tax_table(ps.ob)
tax <- as.data.frame(tax)

# Ensure row names of df1 match column names of df2
matching_indices <- match(colnames(feat), rownames(tax))

# Replace column names in df2 with values from the "Genus" column in df1
colnames(feat)[!is.na(matching_indices)] <- tax$Genus[matching_indices[!is.na(matching_indices)]]

# Run maaslin2
fit_data_random = Maaslin2(input_data     = feat, 
                           input_metadata = metadata_df,
                           max_significance = 0.05,
                           normalization  = "NONE",
                           output         = "Maaslin2_IBD_rel",
                           fixed_effects  = c("IBD_group_name","population"),
                           random_effects = c("continent", "V_region"),
                           reference = "IBD_group_name,HC")

# Read in signficant taxa
df.ibd <- read.table("Maaslin2_IBD_rel/significant_results.tsv", sep = "\t", header = TRUE)

# Keep significant features
df.ibd <- df.ibd[df.ibd$qval <= 0.05, ]
df.ibd$log_qval <- -log10(df.ibd$qval)  # Compute -log10(q_val) for size scaling

# Add a significance threshold line
significance_threshold <- -log10(0.05)

# Create a new column to classify positive and negative coefficients
df.ibd$effect_direction <- ifelse(df.ibd$coef > 0, "Positive", "Negative")

# Define colors
positive_color <- "#00AFBB"
negative_color <- "#FC4E07"

# Volcano Plot
volcano_plot <- ggplot(df.ibd, aes(x = coef, y = log_qval, label = feature)) +
  geom_point(aes(fill = effect_direction), 
             color = "black",
             shape = 21,
             size = 5,
             stroke = 1.2) +
  geom_hline(yintercept = significance_threshold, linetype = "dashed", color = "red") +
  geom_text_repel(data = subset(df.ibd, qval < 0.05), 
                  size = 4, 
                  max.overlaps = 10, 
                  box.padding = 0.4) +
  scale_fill_manual(values = c("Positive" = positive_color, "Negative" = negative_color)) +
  labs(x = "Model Coefficient", 
       y = "-log10(q-value)", 
       fill = "") +
  theme_cowplot(20) +
  theme(
    legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_line(color = "grey85", size = 0.5),
    panel.grid.minor = element_blank()
  )

# OFFSETTING
# Adjust the vertical offset
offset_y <- 0.5 

# Add Annotations with Adjustments
volcano_plot <- volcano_plot +
  # IBD - Draw arrow
  annotate("text", x = max(df.ibd$coef) - 0.1, y = max(df.ibd$log_qval) + 0.2 + offset_y, 
           label = "IBD", size = 10, fontface = "bold", color = positive_color, hjust = 1) +
  geom_segment(aes(x = max(df.ibd$coef), y = max(df.ibd$log_qval) + offset_y, 
                   xend = max(df.ibd$coef) + 0.4, yend = max(df.ibd$log_qval) + offset_y), 
               arrow = arrow(length = unit(0.3, "inches"), type = "closed"), 
               color = positive_color, size = 3, lineend = "round") +
  # HC - Draw arrow
  annotate("text", x = min(df.ibd$coef) + 0.1, y = max(df.ibd$log_qval) + 0.2 + offset_y, 
           label = "HC", size = 10, fontface = "bold", color = negative_color, hjust = 0) +
  geom_segment(aes(x = min(df.ibd$coef), y = max(df.ibd$log_qval) + offset_y, 
                   xend = min(df.ibd$coef) - 0.4, yend = max(df.ibd$log_qval) + offset_y), 
               arrow = arrow(length = unit(0.3, "inches"), type = "closed"), 
               color = negative_color, size = 3, lineend = "round")

# Save plot
ggsave(file.path(DIR_FIG_MAIN, "fig4B_IBD.png"), plot = volcano_plot, width = 10, height = 10, dpi = 300)

# Figures 5C + D - Maaslin2 ---------------------------------------------------
# CD/UC vs HC
ps.ob <- rel_sgn

# extract metadata from phyloseq object
metadata_df <- as.data.frame(sample_data(ps.ob))
metadata_df <- as.matrix(metadata_df)
metadata_df <- as.data.frame(metadata_df)

# Extract OTU table
otu <- otu_table(ps.ob)
feat <- as.data.frame(otu)
feat <- t(feat)

## Tax table
tax <- tax_table(ps.ob)
tax <- as.data.frame(tax)

# Ensure row names of df1 match column names of df2
matching_indices <- match(colnames(feat), rownames(tax))

# Replace column names in df2 with values from the "Genus" column in df1
colnames(feat)[!is.na(matching_indices)] <- tax$Genus[matching_indices[!is.na(matching_indices)]]

# Run maaslin2
fit_data_random = Maaslin2(input_data     = feat, 
                           input_metadata = metadata_df, 
                           max_significance = 0.05,
                           normalization  = "NONE",
                           output         = "Maaslin2_SGN_rel", 
                           fixed_effects  = c("study_group_name","population"),
                           random_effects = c("continent", "V_region"),
                           reference = "study_group_name,HC")

# Read in significant taxa CD
df.sgn <- read.table("Maaslin2_SGN_rel/significant_results.tsv", sep = "\t", header = TRUE)

# Filter significant feature associated with CD
df.sgn <- df.sgn %>%
  dplyr::filter(value == "CD")

# Keep significant features
df.sgn <- df.sgn[df.sgn$qval <= 0.05, ]
df.sgn$log_qval <- -log10(df.sgn$qval)  # Compute -log10(q_val) for size scaling

# Add a significance threshold line
significance_threshold <- -log10(0.05)

# Create a new column to classify positive and negative coefficients
df.sgn$effect_direction <- ifelse(df.sgn$coef > 0, "Positive", "Negative")

# Define colors
positive_color <- "#00AFBB"
negative_color <- "#FC4E07"

# Volcano Plot
volcano_plot <- ggplot(df.sgn, aes(x = coef, y = log_qval, label = feature)) +
  geom_point(aes(fill = effect_direction), 
             color = "black",
             shape = 21,
             size = 5,
             stroke = 1.2) +
  geom_hline(yintercept = significance_threshold, linetype = "dashed", color = "red") +
  geom_text_repel(data = subset(df.sgn, qval < 0.05), 
                  size = 4, 
                  max.overlaps = 10, 
                  box.padding = 0.4) +
  scale_fill_manual(values = c("Positive" = positive_color, "Negative" = negative_color)) +
  labs(x = "Model Coefficient", 
       y = "-log10(q-value)", 
       fill = "") +
  theme_cowplot(20) +
  theme(
    legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_line(color = "grey85", size = 0.5),
    panel.grid.minor = element_blank()
  )

# OFFSETTING
# Adjust the vertical offset
offset_y <- 0.5  # Adjust this value to control how much you move it upwards

# Add Annotations
volcano_plot <- volcano_plot +
  # CD - Draw arrow
  annotate("text", x = max(df.sgn$coef) - 0.1, y = max(df.sgn$log_qval) + 0.2 + offset_y, 
           label = "CD", size = 10, fontface = "bold", color = positive_color, hjust = 1) +
  geom_segment(aes(x = max(df.sgn$coef), y = max(df.sgn$log_qval) + offset_y, 
                   xend = max(df.sgn$coef) + 0.4, yend = max(df.sgn$log_qval) + offset_y), 
               arrow = arrow(length = unit(0.3, "inches"), type = "closed"), 
               color = positive_color, size = 3, lineend = "round") +  
  # HC - Draw arrow
  annotate("text", x = min(df.sgn$coef) + 0.1, y = max(df.sgn$log_qval) + 0.2 + offset_y, 
           label = "HC", size = 10, fontface = "bold", color = negative_color, hjust = 0) +
  geom_segment(aes(x = min(df.sgn$coef), y = max(df.sgn$log_qval) + offset_y, 
                   xend = min(df.sgn$coef) - 0.4, yend = max(df.sgn$log_qval) + offset_y), 
               arrow = arrow(length = unit(0.3, "inches"), type = "closed"), 
               color = negative_color, size = 3, lineend = "round")

# Save the plot
ggsave(file.path(DIR_FIG_MAIN, "fig5C_CD.png"), plot = volcano_plot, width = 10, height = 10, dpi = 300)


# Read in significant taxa UC
df.sgn <- read.table("Maaslin2_SGN_rel/significant_results.tsv", sep = "\t", header = TRUE)

# Filter significant feature associated with UC
df.sgn <- df.sgn %>%
  dplyr::filter(value == "UC")

# Keep significant features
df.sgn <- df.sgn[df.sgn$qval <= 0.05, ]
df.sgn$log_qval <- -log10(df.sgn$qval)  # Compute -log10(q_val) for size scaling

# Add a significance threshold line
significance_threshold <- -log10(0.05)

# Create a new column to classify positive and negative coefficients
df.sgn$effect_direction <- ifelse(df.sgn$coef > 0, "Positive", "Negative")

# Define colors
positive_color <- "#00AFBB"
negative_color <- "#FC4E07"

# Volcano Plot
volcano_plot <- ggplot(df.sgn, aes(x = coef, y = log_qval, label = feature)) +
  geom_point(aes(fill = effect_direction), 
             color = "black",
             shape = 21,
             size = 5,
             stroke = 1.2) +
  geom_hline(yintercept = significance_threshold, linetype = "dashed", color = "red") +
  geom_text_repel(data = subset(df.sgn, qval < 0.05), 
                  size = 4, 
                  max.overlaps = 10, 
                  box.padding = 0.4) +
  scale_fill_manual(values = c("Positive" = positive_color, "Negative" = negative_color)) +
  labs(x = "Model Coefficient", 
       y = "-log10(q-value)", 
       fill = "") +
  theme_cowplot(20) +
  theme(
    legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_line(color = "grey85", size = 0.5),
    panel.grid.minor = element_blank()
  )

# OFFSETTING
# Adjust the vertical offset
offset_y <- 0.5  # Adjust this value to control how much you move it upwards

# Add Annotations
volcano_plot <- volcano_plot +
  # UC - Draw arrow
  annotate("text", x = max(df.sgn$coef) - 0.1, y = max(df.sgn$log_qval) + 0.2 + offset_y, 
           label = "UC", size = 10, fontface = "bold", color = positive_color, hjust = 1) +
  geom_segment(aes(x = max(df.sgn$coef), y = max(df.sgn$log_qval) + offset_y, 
                   xend = max(df.sgn$coef) + 0.4, yend = max(df.sgn$log_qval) + offset_y), 
               arrow = arrow(length = unit(0.3, "inches"), type = "closed"), 
               color = positive_color, size = 3, lineend = "round") + 
  # HC - Draw arrow
  annotate("text", x = min(df.sgn$coef) + 0.1, y = max(df.sgn$log_qval) + 0.2 + offset_y, 
           label = "HC", size = 10, fontface = "bold", color = negative_color, hjust = 0) +
  geom_segment(aes(x = min(df.sgn$coef), y = max(df.sgn$log_qval) + offset_y, 
                   xend = min(df.sgn$coef) - 0.4, yend = max(df.sgn$log_qval) + offset_y), 
               arrow = arrow(length = unit(0.3, "inches"), type = "closed"), 
               color = negative_color, size = 3, lineend = "round") 

# Save the plot
ggsave(file.path(DIR_FIG_MAIN, "fig5D_UC.png"), plot = volcano_plot, width = 10, height = 10, dpi = 300)
