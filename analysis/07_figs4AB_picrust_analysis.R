## ===============================
## analysis/07_figs4AB_diversity_analysis.R
## Perform diversity analysis and plot

source("analysis/00_config.R")
source("R/packages.R")
source("R/functions_phyloseq.R")

# Load and clean data
# metadata
metadata <- read.table(FILE_PS_META, sep = "\t", header = TRUE, row.names = 1)

# Filter for IBD sub-type
metadata_df <- metadata %>%
  filter(study_group_name != "IBD")

# Picrust metacyc output
metacyc <- fread(FILE_PICRUST_METACYC, sep = "\t", data.table = FALSE)

# remove full description column
metacyc_filt <- metacyc[, -1]

# Filter data to keep IBD sub-type samples
metacyc_filt <- metacyc_filt[, colnames(metacyc_filt) %in% rownames(metadata_df)]

# Run maaslin2 ------------------------------------------------------------

fit_out <- Maaslin2(input_data     = metacyc_filt, 
                    input_metadata = metadata_df,
                    max_significance = 0.05,
                    normalization = "TSS",
                    transform = "LOG",
                    output         = "SGN_rel_metacyc",
                    fixed_effects  = c("study_group_name","Population"),
                    random_effects = c("continent", "V_region"),
                    reference = "study_group_name,HC",
                    cores = 8)

# Figure 4A - Picrust metacyc CD vs HC -----------------------------------

# read significant results in
df <- read.table("SGN_rel_metacyc/significant_results.tsv", 
                 sep = "\t", header = TRUE)

# change number of pathway to match exisiting df
df$feature[df$feature == "X3.HYDROXYPHENYLACETATE.DEGRADATION.PWY"] <- 
  "3.HYDROXYPHENYLACETATE.DEGRADATION.PWY"

# change "." for hyphen
df$feature <- gsub("\\.", "-", df$feature)

# Filter significant feature associated with CD
df <- df %>%
  dplyr::filter(value == "CD")

# Keep significant features <= 0.05
df <- df[df$qval <= 0.05, ]

# Compute -log10(q_val) for size scaling
df$log_qval <- -log10(df$qval)  

# Add a significance threshold line
significance_threshold <- -log10(0.05)

# Create a new column to classify positive and negative coefficients
df$effect_direction <- ifelse(df$coef > 0, "Positive", "Negative")

# pathway mapping file
mapping = metacyc

mapping <- mapping %>%
  dplyr::select(pathway, description)

mapping <- mapping %>%
  rename(feature = pathway)

# adding full pathway description
df <- df %>%
  left_join(mapping, by = "feature")

# filtering for pathways above and below 0.8 threshold
df <- df %>%
  filter(coef >= 0.8 | coef <= -0.8)

df$description <- factor(
  df$description,
  levels = df$description[order(df$coef, decreasing = FALSE)]
)

# Define colors
colours <- c(
  "HC" = "#648fff",
  "IBD" = "#785ef0",
  "CD" = "#dc267f",
  "UC" = "#ffb000"
)

positive_color <- "#dc267f"
negative_color <- "#648fff"

CD <- "#dc267f"
HC  <- "#648fff"

# Plotting data
p_beta <- ggplot(df, aes(x = coef, y = description)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_segment(
    aes(x = 0, xend = coef, y = description, yend = description),
    linewidth = 0.5,
    colour = "grey70"
  ) +
  geom_point(
    aes(fill = effect_direction, size = log_qval),
    shape = 21,
    colour = "black",
    stroke = 0.8
  ) +
  scale_fill_manual(
    values = c(
      Positive = positive_color,
      Negative = negative_color
    ),
    guide = "none"
  ) +
  scale_size_continuous(range = c(2, 6), name = expression(-log[10](FDR))) +
  labs(
    x = expression(beta ~ "coefficient"),
    y = NULL
  ) +
  theme_cowplot(16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

legend_plot <- ggplot(
  data.frame(group = c("CD", "HC"), x = 1, y = c(2, 1)),
  aes(x = x, y = y, fill = group)
) +
  geom_point(shape = 21, size = 4, colour = "black", stroke = 0.8) +
  scale_fill_manual(
    name = NULL,
    values = c(CD = CD, HC = HC)
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0.1, "cm"),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size = 11)
  )

# save plot
ggsave(file.path(DIR_FIG_MAIN, "poster_metacyc_CD.png"), plot = p_beta,  width = 12, height = 6, dpi = 300)

# Figure 4B - Picrust metacyc UC vs HC -----------------------------------

# read significant results in
df <- read.table("SGN_rel_metacyc/significant_results.tsv", 
                 sep = "\t", header = TRUE)

# change number of pathway to match exisiting df
df$feature[df$feature == "X3.HYDROXYPHENYLACETATE.DEGRADATION.PWY"] <- 
  "3.HYDROXYPHENYLACETATE.DEGRADATION.PWY"

# change "." for hyphen
df$feature <- gsub("\\.", "-", df$feature)

# Filter significant feature associated with UC
df <- df %>%
  dplyr::filter(value == "UC")

# Keep significant features <= 0.05
df <- df[df$qval <= 0.05, ]

# Compute -log10(q_val) for size scaling
df$log_qval <- -log10(df$qval)  

# Add a significance threshold line
significance_threshold <- -log10(0.05)

# Create a new column to classify positive and negative coefficients
df$effect_direction <- ifelse(df$coef > 0, "Positive", "Negative")

# pathway mapping file
mapping = metacyc

mapping <- mapping %>%
  dplyr::select(pathway, description)

mapping <- mapping %>%
  rename(feature = pathway)

# adding full pathway description
df <- df %>%
  left_join(mapping, by = "feature")

df$description <- factor(
  df$description,
  levels = df$description[order(df$coef, decreasing = FALSE)]
)

# Define colors
colours <- c(
  "HC" = "#648fff",
  "IBD" = "#785ef0",
  "CD" = "#dc267f",
  "UC" = "#ffb000"
)

positive_color <- "#ffb000"
negative_color <- "#648fff"

UC <- "#ffb000"
HC  <- "#648fff"

# Plotting data
p_beta <- ggplot(df, aes(x = coef, y = description)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_segment(
    aes(x = 0, xend = coef, y = description, yend = description),
    linewidth = 0.5,
    colour = "grey70"
  ) +
  geom_point(
    aes(fill = effect_direction, size = log_qval),
    shape = 21,
    colour = "black",
    stroke = 0.8
  ) +
  scale_fill_manual(
    values = c(
      Positive = positive_color,
      Negative = negative_color
    ),
    guide = "none"
  ) +
  scale_size_continuous(range = c(2, 6), name = expression(-log[10](FDR))) +
  labs(
    x = expression(beta ~ "coefficient"),
    y = NULL
  ) +
  theme_cowplot(14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

legend_plot <- ggplot(
  data.frame(group = c("UC", "HC"), x = 1, y = c(2, 1)),
  aes(x = x, y = y, fill = group)
) +
  geom_point(shape = 21, size = 4, colour = "black", stroke = 0.8) +
  scale_fill_manual(
    name = NULL,
    values = c(UC = UC, HC = HC)
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0.1, "cm"),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size = 11)
  )

# save plot
ggsave(file.path(DIR_FIG_MAIN, "fig4B_metacyc_UC.png"), plot = p_beta,  width = 9, height = 8, dpi = 300)
