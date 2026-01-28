## ===============================
## analysis/02_fig1B_study_sample_details.R
## Create fig1B showing study details

source("analysis/00_config.R")
source("R/packages.R")
source("R/functions_phyloseq.R")

# Load data
metadata <- read.table(FILE_PS_META, sep = "\t", header = TRUE, row.names = 1)

# Figure 1B - Study sample details  ----------------------------------------
facet_order <- c("study_group_name", "sample_type", "V_region", "continent", "Population")

# Reshape data to long format for faceting
metadata_long <- metadata %>%
  pivot_longer(cols = all_of(facet_order), names_to = "Facet_Variable", values_to = "Facet_Value") %>%
  mutate(
    Facet_Variable = factor(Facet_Variable, levels = facet_order),  # Order facets
    Facet_Value = fct_infreq(Facet_Value)  # Order bars within each facet by frequency (largest to smallest)
  )

set.seed(4)

# 1) Count per category per study + study numbers
counts <- metadata_long %>%
  group_by(Facet_Variable, Facet_Value, study_alias) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(study_number = as.integer(factor(study_alias)))

# 2) Create facet-specific numeric coding for Facet_Value
counts_j <- counts %>%
  group_by(Facet_Variable) %>%
  mutate(
    Facet_Value_factor  = factor(Facet_Value, levels = unique(Facet_Value)),
    Facet_Value_numeric = as.numeric(Facet_Value_factor),
    Facet_Value_jittered = Facet_Value_numeric + runif(n(), -0.2, 0.2)
  ) %>%
  ungroup()

# Legend lookup
lookup <- counts_j %>%
  distinct(study_number, study_alias) %>%
  arrange(study_number)

# 3) Plot
demo_dots <- ggplot(counts_j, aes(x = Facet_Value_jittered, y = n)) +
  geom_point(aes(color = factor(study_number)), alpha = 0, size = 4, show.legend = TRUE) +
  geom_point(shape = 21, size = 8, fill = "steelblue", color = "black") +
  geom_text(aes(label = study_number), color = "white", size = 3.8) +
  facet_grid(~ Facet_Variable, scales = "free_x", space = "free") +
  theme_cowplot(20) +
  labs(x = NULL, y = "Number of samples") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.background = element_rect(color = "white", fill = "white", size = 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, 400), clip = "off") +
  scale_color_discrete(
    name = "Studies",
    breaks = as.character(lookup$study_number),
    labels = paste0(lookup$study_number, " = ", lookup$study_alias)
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4))) +
  scale_x_continuous(
    breaks = sort(unique(counts_j$Facet_Value_numeric)),
    labels = levels(counts_j$Facet_Value_factor),
    expand = expansion(mult = c(0.1, 0.1))
  )

# Save to file
ggsave(file.path(DIR_FIG_MAIN, "fig1B_demo.png"), plot = demo_dots, width = 14, height = 5, dpi = 300)

# Legend
legend_df <- lookup %>%
  distinct(study_number, study_alias) %>%
  arrange(as.numeric(study_number)) %>%
  mutate(y = rev(seq_along(study_number)))  # vertical order

# Legend plot
legend_plot <- ggplot(legend_df, aes(x = 1, y = y)) +
  geom_point(color = "steelblue", size = 10) +
  geom_text(aes(label = study_number), color = "white", size = 5) +
  geom_text(aes(x = 1.1, label = study_alias), hjust = 0, size = 7) +
  scale_y_continuous(NULL, breaks = NULL) +
  scale_x_continuous(NULL, breaks = NULL, limits = c(0.8, 2)) +
  theme_void(base_size = 14) +
  theme(plot.margin = margin(10, 10, 10, 10))

# Save plot to file
ggsave(file.path(DIR_FIG_MAIN, "fig1B_demo_legend.png"), plot = legend_plot, width = 5, height = 4, dpi = 300)

