## ===============================
## analysis/03_fig3A_hedges_smd.R
## Perform hedges smd on alpha diversity metrics

source("analysis/00_config.R")
source("R/packages.R")
source("R/functions_phyloseq.R")

# Load data
metadata <- read.table(FILE_PS_META, sep = "\t", header = TRUE, row.names = 1)
rel <- readRDS(FILE_PS_PROCESSED_REL)
rel_sgn <- readRDS(FILE_PS_PROCESSED_REL_SGN)

# Figure 3A - Hedges SMD  ------------------------------------------------------
# Shannon diversity
# Calculate alpha diversity metrics
alpha <- microbiome::alpha(rel, index = "all")
merge <- merge(metadata, alpha, by = "row.names")

# IBD vs HC
summary_df <- merge %>%
  group_by(study_alias, IBD_group_name) %>%
  summarise(
    n = n(),
    mean = mean(diversity_shannon),
    sd = sd(diversity_shannon),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = IBD_group_name, values_from = c(n, mean, sd),
              names_glue = "{.value}.{IBD_group_name}")

# Calculate Hedges SMD
meta_result <- metacont(n.IBD, mean.IBD, sd.IBD,
                        n.HC, mean.HC, sd.HC,
                        data = summary_df,
                        sm = "SMD",
                        method.smd = "Hedges",
                        studlab = study_alias)

# Create table of results
meta_df <- data.frame(
  study_alias = meta_result[["studlab"]],
  Hedges_g = meta_result[["TE"]],
  se = meta_result[["seTE"]],
  CI_lower = meta_result[["lower"]],
  CI_upper = meta_result[["upper"]],
  p_value = meta_result[["pval"]],
  n_e = meta_result[["n.e"]],
  n_c = meta_result[["n.c"]],
  mean_e = meta_result[["mean.e"]],
  mean_c = meta_result[["mean.c"]],
  sd_e = meta_result[["sd.e"]],
  sd_c = meta_result[["sd.c"]]
)

# Overall mode
overall_row <- data.frame(
  study_alias = "Overall",
  Hedges_g = meta_result$TE.random,
  se = meta_result$seTE.random,
  CI_lower = meta_result$lower.random,
  CI_upper = meta_result$upper.random,
  p_value = meta_result$pval.random,
  n_e = NA,
  n_c = NA,
  mean_e = NA,
  mean_c = NA,
  sd_e = NA,
  sd_c = NA
)

meta_df_with_overall_IBD_HC <- rbind(meta_df, overall_row)

meta_df_with_overall_IBD_HC$comparison <- "IBD vs HC"

# CD vs HC

merge_cd <- merge %>%
  filter(study_group_name %in% c("CD", "HC")) %>%
  group_by(study_alias) %>%
  filter(all(c("CD", "HC") %in% study_group_name)) %>%
  ungroup()

summary_df <- merge_cd %>%
  group_by(study_alias, study_group_name) %>%
  summarise(
    n = n(),
    mean = mean(diversity_shannon),
    sd = sd(diversity_shannon),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = study_group_name, values_from = c(n, mean, sd),
              names_glue = "{.value}.{study_group_name}")

meta_result <- metacont(n.CD, mean.CD, sd.CD,
                        n.HC, mean.HC, sd.HC,
                        data = summary_df,
                        sm = "SMD",
                        method.smd = "Hedges",
                        studlab = study_alias)

meta_df <- data.frame(
  study_alias = meta_result[["studlab"]],
  Hedges_g = meta_result[["TE"]],
  se = meta_result[["seTE"]],
  CI_lower = meta_result[["lower"]],
  CI_upper = meta_result[["upper"]],
  p_value = meta_result[["pval"]],
  n_e = meta_result[["n.e"]],
  n_c = meta_result[["n.c"]],
  mean_e = meta_result[["mean.e"]],
  mean_c = meta_result[["mean.c"]],
  sd_e = meta_result[["sd.e"]],
  sd_c = meta_result[["sd.c"]]
)

overall_row <- data.frame(
  study_alias = "Overall",
  Hedges_g = meta_result$TE.random,
  se = meta_result$seTE.random,
  CI_lower = meta_result$lower.random,
  CI_upper = meta_result$upper.random,
  p_value = meta_result$pval.random,
  n_e = NA,
  n_c = NA,
  mean_e = NA,
  mean_c = NA,
  sd_e = NA,
  sd_c = NA
)

meta_df_with_overall_CD_HC <- rbind(meta_df, overall_row)

meta_df_with_overall_CD_HC$comparison <- "CD vs HC"

# UC vs HC

merge_uc <- merge %>%
  filter(study_group_name %in% c("UC", "HC")) %>%
  group_by(study_alias) %>%
  filter(all(c("UC", "HC") %in% study_group_name)) %>%
  ungroup()

summary_df <- merge_uc %>%
  group_by(study_alias, study_group_name) %>%
  summarise(
    n = n(),
    mean = mean(diversity_shannon),
    sd = sd(diversity_shannon),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = study_group_name, values_from = c(n, mean, sd),
              names_glue = "{.value}.{study_group_name}")

meta_result <- metacont(n.UC, mean.UC, sd.UC,
                        n.HC, mean.HC, sd.HC,
                        data = summary_df,
                        sm = "SMD",
                        method.smd = "Hedges",
                        studlab = study_alias)

meta_df <- data.frame(
  study_alias = meta_result[["studlab"]],
  Hedges_g = meta_result[["TE"]],
  se = meta_result[["seTE"]],
  CI_lower = meta_result[["lower"]],
  CI_upper = meta_result[["upper"]],
  p_value = meta_result[["pval"]],
  n_e = meta_result[["n.e"]],
  n_c = meta_result[["n.c"]],
  mean_e = meta_result[["mean.e"]],
  mean_c = meta_result[["mean.c"]],
  sd_e = meta_result[["sd.e"]],
  sd_c = meta_result[["sd.c"]]
)

overall_row <- data.frame(
  study_alias = "Overall",
  Hedges_g = meta_result$TE.random,
  se = meta_result$seTE.random,
  CI_lower = meta_result$lower.random,
  CI_upper = meta_result$upper.random,
  p_value = meta_result$pval.random,
  n_e = NA,
  n_c = NA,
  mean_e = NA,
  mean_c = NA,
  sd_e = NA,
  sd_c = NA
)

meta_df_with_overall_UC_HC <- rbind(meta_df, overall_row)

meta_df_with_overall_UC_HC$comparison <- "UC vs HC"

# CD vs UC

merge_cd_uc <- merge %>%
  filter(study_group_name %in% c("CD", "UC")) %>%
  group_by(study_alias) %>%
  filter(all(c("CD", "UC") %in% study_group_name)) %>%
  ungroup()

summary_df <- merge_cd_uc %>%
  group_by(study_alias, study_group_name) %>%
  summarise(
    n = n(),
    mean = mean(diversity_shannon),
    sd = sd(diversity_shannon),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = study_group_name, values_from = c(n, mean, sd),
              names_glue = "{.value}.{study_group_name}")

meta_result <- metacont(n.CD, mean.CD, sd.CD,
                        n.UC, mean.UC, sd.UC,
                        data = summary_df,
                        sm = "SMD",
                        method.smd = "Hedges",
                        studlab = study_alias)

meta_df <- data.frame(
  study_alias = meta_result[["studlab"]],
  Hedges_g = meta_result[["TE"]],
  se = meta_result[["seTE"]],
  CI_lower = meta_result[["lower"]],
  CI_upper = meta_result[["upper"]],
  p_value = meta_result[["pval"]],
  n_e = meta_result[["n.e"]],
  n_c = meta_result[["n.c"]],
  mean_e = meta_result[["mean.e"]],
  mean_c = meta_result[["mean.c"]],
  sd_e = meta_result[["sd.e"]],
  sd_c = meta_result[["sd.c"]]
)

overall_row <- data.frame(
  study_alias = "Overall",
  Hedges_g = meta_result$TE.random,
  se = meta_result$seTE.random,
  CI_lower = meta_result$lower.random,
  CI_upper = meta_result$upper.random,
  p_value = meta_result$pval.random,
  n_e = NA,
  n_c = NA,
  mean_e = NA,
  mean_c = NA,
  sd_e = NA,
  sd_c = NA
)

meta_df_with_overall_CD_UC <- rbind(meta_df, overall_row)

meta_df_with_overall_CD_UC$comparison <- "CD vs UC"

## All dfs combine ###
meta_output_all <- rbind(meta_df_with_overall_IBD_HC, meta_df_with_overall_CD_HC, meta_df_with_overall_UC_HC, meta_df_with_overall_CD_UC)

# Plotting
df = meta_output_all

# 1. Extract "Overall" rows
overall_df <- df %>%
  dplyr::filter(study_alias == "Overall") %>%
  dplyr::select(
    comparison,
    Hedges_g_overall = Hedges_g,
    CI_lower = CI_lower,
    CI_upper = CI_upper,
    p_value_overall = p_value
  )

# 2. Filter to only study-level rows
hedges_all <- df %>%
  filter(study_alias != "Overall")

# 3. Join overall data back to study rows (for plotting)
hedges_all <- left_join(hedges_all, overall_df, by = "comparison")

# 4. Set desired y-axis (comparison) order
desired_order <- c("CD vs UC","UC vs HC","CD vs HC","IBD vs HC")

hedges_all <- hedges_all %>%
  mutate(
    comparison = factor(comparison, levels = desired_order)
  )

# 5. Add jitter and consistent study number
set.seed(4)
hedges_all_jittered <- hedges_all %>%
  mutate(
    comparison_factor = comparison,
    comparison_numeric = as.numeric(comparison_factor),
    comparison_jittered = comparison_numeric + runif(n(), -0.2, 0.2),
    study_number = as.numeric(factor(study_alias))  # consistent per study
  )

# 6. Prepare data for overall diamonds and CI lines
overall_plot_df <- hedges_all_jittered %>%
  distinct(comparison, comparison_numeric, Hedges_g_overall, CI_lower.y, CI_upper.y, p_value_overall) %>%
  rename(y = comparison_numeric) %>%
  mutate(
    sig_fill = ifelse(p_value_overall < 0.05,  "#00BFC4", "#F8766D")
  )

# 7. Plot
plt <- ggplot(hedges_all_jittered, aes(x = Hedges_g, y = comparison_jittered)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  # Study-level points
  geom_point(size = 10, color = "steelblue") +
  geom_text(aes(label = study_number), color = "white", size = 5) +
  # Overall 95% CI bars
  geom_segment(
    data = overall_plot_df,
    aes(x = CI_lower.y, xend = CI_upper.y, y = y, yend = y),
    inherit.aes = FALSE,
    color = "black", size = 1.2
  ) +
  # Overall Hedges' g point coloured by significance
  geom_point(
    data = overall_plot_df,
    aes(x = Hedges_g_overall, y = y, fill = I(sig_fill)),
    shape = 23, size = 5, color = "black", stroke = 1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  # Formatting axes and theme
  scale_y_continuous(
    breaks = sort(unique(hedges_all_jittered$comparison_numeric)),
    labels = levels(hedges_all_jittered$comparison_factor),
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  coord_cartesian(clip = "off") +
  theme_cowplot(20) +
  labs(x = "Hedges' SMD", y = NULL, title = "Shannon index") +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

# Save plot
ggsave(file.path(DIR_FIG_MAIN, "Fig3A_hedges_shannon.png"), plot = plt, width = 8, height = 4, dpi = 300)

# Richness
# IBD vs HC 
summary_df <- merge %>%
  group_by(study_alias, IBD_group_name) %>%
  summarise(
    n = n(),
    mean = mean(observed),
    sd = sd(observed),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = IBD_group_name, values_from = c(n, mean, sd),
              names_glue = "{.value}.{IBD_group_name}")

meta_result <- metacont(n.IBD, mean.IBD, sd.IBD,
                        n.HC, mean.HC, sd.HC,
                        data = summary_df,
                        sm = "SMD",
                        method.smd = "Hedges",
                        studlab = study_alias)

meta_df <- data.frame(
  study_alias = meta_result[["studlab"]],
  Hedges_g = meta_result[["TE"]],
  se = meta_result[["seTE"]],
  CI_lower = meta_result[["lower"]],
  CI_upper = meta_result[["upper"]],
  p_value = meta_result[["pval"]],
  n_e = meta_result[["n.e"]],
  n_c = meta_result[["n.c"]],
  mean_e = meta_result[["mean.e"]],
  mean_c = meta_result[["mean.c"]],
  sd_e = meta_result[["sd.e"]],
  sd_c = meta_result[["sd.c"]]
)

overall_row <- data.frame(
  study_alias = "Overall",
  Hedges_g = meta_result$TE.random,
  se = meta_result$seTE.random,
  CI_lower = meta_result$lower.random,
  CI_upper = meta_result$upper.random,
  p_value = meta_result$pval.random,
  n_e = NA,
  n_c = NA,
  mean_e = NA,
  mean_c = NA,
  sd_e = NA,
  sd_c = NA
)

meta_df_with_overall_IBD_HC <- rbind(meta_df, overall_row)

meta_df_with_overall_IBD_HC$comparison <- "IBD vs HC"

# CD vs HC

merge_cd <- merge %>%
  filter(study_group_name %in% c("CD", "HC")) %>%
  group_by(study_alias) %>%
  filter(all(c("CD", "HC") %in% study_group_name)) %>%
  ungroup()

summary_df <- merge_cd %>%
  group_by(study_alias, study_group_name) %>%
  summarise(
    n = n(),
    mean = mean(observed),
    sd = sd(observed),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = study_group_name, values_from = c(n, mean, sd),
              names_glue = "{.value}.{study_group_name}")

meta_result <- metacont(n.CD, mean.CD, sd.CD,
                        n.HC, mean.HC, sd.HC,
                        data = summary_df,
                        sm = "SMD",
                        method.smd = "Hedges",
                        studlab = study_alias)

meta_df <- data.frame(
  study_alias = meta_result[["studlab"]],
  Hedges_g = meta_result[["TE"]],
  se = meta_result[["seTE"]],
  CI_lower = meta_result[["lower"]],
  CI_upper = meta_result[["upper"]],
  p_value = meta_result[["pval"]],
  n_e = meta_result[["n.e"]],
  n_c = meta_result[["n.c"]],
  mean_e = meta_result[["mean.e"]],
  mean_c = meta_result[["mean.c"]],
  sd_e = meta_result[["sd.e"]],
  sd_c = meta_result[["sd.c"]]
)

overall_row <- data.frame(
  study_alias = "Overall",
  Hedges_g = meta_result$TE.random,
  se = meta_result$seTE.random,
  CI_lower = meta_result$lower.random,
  CI_upper = meta_result$upper.random,
  p_value = meta_result$pval.random,
  n_e = NA,
  n_c = NA,
  mean_e = NA,
  mean_c = NA,
  sd_e = NA,
  sd_c = NA
)

meta_df_with_overall_CD_HC <- rbind(meta_df, overall_row)

meta_df_with_overall_CD_HC$comparison <- "CD vs HC"

# UC vs HC

merge_uc <- merge %>%
  filter(study_group_name %in% c("UC", "HC")) %>%
  group_by(study_alias) %>%
  filter(all(c("UC", "HC") %in% study_group_name)) %>%
  ungroup()

summary_df <- merge_uc %>%
  group_by(study_alias, study_group_name) %>%
  summarise(
    n = n(),
    mean = mean(observed),
    sd = sd(observed),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = study_group_name, values_from = c(n, mean, sd),
              names_glue = "{.value}.{study_group_name}")

meta_result <- metacont(n.UC, mean.UC, sd.UC,
                        n.HC, mean.HC, sd.HC,
                        data = summary_df,
                        sm = "SMD",
                        method.smd = "Hedges",
                        studlab = study_alias)

meta_df <- data.frame(
  study_alias = meta_result[["studlab"]],
  Hedges_g = meta_result[["TE"]],
  se = meta_result[["seTE"]],
  CI_lower = meta_result[["lower"]],
  CI_upper = meta_result[["upper"]],
  p_value = meta_result[["pval"]],
  n_e = meta_result[["n.e"]],
  n_c = meta_result[["n.c"]],
  mean_e = meta_result[["mean.e"]],
  mean_c = meta_result[["mean.c"]],
  sd_e = meta_result[["sd.e"]],
  sd_c = meta_result[["sd.c"]]
)

overall_row <- data.frame(
  study_alias = "Overall",
  Hedges_g = meta_result$TE.random,
  se = meta_result$seTE.random,
  CI_lower = meta_result$lower.random,
  CI_upper = meta_result$upper.random,
  p_value = meta_result$pval.random,
  n_e = NA,
  n_c = NA,
  mean_e = NA,
  mean_c = NA,
  sd_e = NA,
  sd_c = NA
)

meta_df_with_overall_UC_HC <- rbind(meta_df, overall_row)

meta_df_with_overall_UC_HC$comparison <- "UC vs HC"

# CD vs UC

merge_cd_uc <- merge %>%
  filter(study_group_name %in% c("CD", "UC")) %>%
  group_by(study_alias) %>%
  filter(all(c("CD", "UC") %in% study_group_name)) %>%
  ungroup()

summary_df <- merge_cd_uc %>%
  group_by(study_alias, study_group_name) %>%
  summarise(
    n = n(),
    mean = mean(observed),
    sd = sd(observed),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = study_group_name, values_from = c(n, mean, sd),
              names_glue = "{.value}.{study_group_name}")

meta_result <- metacont(n.CD, mean.CD, sd.CD,
                        n.UC, mean.UC, sd.UC,
                        data = summary_df,
                        sm = "SMD",
                        method.smd = "Hedges",
                        studlab = study_alias)

meta_df <- data.frame(
  study_alias = meta_result[["studlab"]],
  Hedges_g = meta_result[["TE"]],
  se = meta_result[["seTE"]],
  CI_lower = meta_result[["lower"]],
  CI_upper = meta_result[["upper"]],
  p_value = meta_result[["pval"]],
  n_e = meta_result[["n.e"]],
  n_c = meta_result[["n.c"]],
  mean_e = meta_result[["mean.e"]],
  mean_c = meta_result[["mean.c"]],
  sd_e = meta_result[["sd.e"]],
  sd_c = meta_result[["sd.c"]]
)

overall_row <- data.frame(
  study_alias = "Overall",
  Hedges_g = meta_result$TE.random,
  se = meta_result$seTE.random,
  CI_lower = meta_result$lower.random,
  CI_upper = meta_result$upper.random,
  p_value = meta_result$pval.random,
  n_e = NA,
  n_c = NA,
  mean_e = NA,
  mean_c = NA,
  sd_e = NA,
  sd_c = NA
)

meta_df_with_overall_CD_UC <- rbind(meta_df, overall_row)

meta_df_with_overall_CD_UC$comparison <- "CD vs UC"

## All dfs combine ###
meta_output_all <- rbind(meta_df_with_overall_IBD_HC, meta_df_with_overall_CD_HC, meta_df_with_overall_UC_HC, meta_df_with_overall_CD_UC)

# Plotting
df = meta_output_all

# 1. Extract "Overall" rows
overall_df <- df %>%
  dplyr::filter(study_alias == "Overall") %>%
  dplyr::select(
    comparison,
    Hedges_g_overall = Hedges_g,
    CI_lower = CI_lower,
    CI_upper = CI_upper,
    p_value_overall = p_value
  )

# 2. Filter to only study-level rows
hedges_all <- df %>%
  filter(study_alias != "Overall")

# 3. Join overall data back to study rows (for plotting)
hedges_all <- left_join(hedges_all, overall_df, by = "comparison")

# 4. Set desired y-axis (comparison) order
desired_order <- c("CD vs UC","UC vs HC","CD vs HC","IBD vs HC")

hedges_all <- hedges_all %>%
  mutate(
    comparison = factor(comparison, levels = desired_order)
  )

# 5. Add jitter and consistent study number
set.seed(4)
hedges_all_jittered <- hedges_all %>%
  mutate(
    comparison_factor = comparison,
    comparison_numeric = as.numeric(comparison_factor),
    comparison_jittered = comparison_numeric + runif(n(), -0.2, 0.2),
    study_number = as.numeric(factor(study_alias))  # consistent per study
  )

# 6. Prepare data for overall diamonds and CI lines
overall_plot_df <- hedges_all_jittered %>%
  distinct(comparison, comparison_numeric, Hedges_g_overall, CI_lower.y, CI_upper.y, p_value_overall) %>%
  rename(y = comparison_numeric) %>%
  mutate(
    sig_fill = ifelse(p_value_overall < 0.05,  "#00BFC4", "#F8766D")
  )

# 7. Plot
plt <- ggplot(hedges_all_jittered, aes(x = Hedges_g, y = comparison_jittered)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  # Study-level points
  geom_point(size = 10, color = "steelblue") +
  geom_text(aes(label = study_number), color = "white", size = 5) +
  # Overall 95% CI bars
  geom_segment(
    data = overall_plot_df,
    aes(x = CI_lower.y, xend = CI_upper.y, y = y, yend = y),
    inherit.aes = FALSE,
    color = "black", size = 1.2
  ) +
  # Overall Hedges' g point (diamond) colored by significance
  geom_point(
    data = overall_plot_df,
    aes(x = Hedges_g_overall, y = y, fill = I(sig_fill)),
    shape = 23, size = 5, color = "black", stroke = 1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  # Formatting axes and theme
  scale_y_continuous(
    breaks = sort(unique(hedges_all_jittered$comparison_numeric)),
    labels = levels(hedges_all_jittered$comparison_factor),
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  coord_cartesian(clip = "off") +
  theme_cowplot(20) +
  labs(x = "Hedges' SMD", y = NULL, title = "Richness") +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5)  # Center the title
  )

# Save plot
ggsave(file.path(DIR_FIG_MAIN, "Fig3A_hedges_richness.png"), plot = plt, width = 8, height = 4, dpi = 300)

# Legend
legend_df <- hedges_all_jittered %>%
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

# Save legend plot
ggsave(file.path(DIR_FIG_MAIN, "Fig3A_legend_plot.png"), plot = legend_plot, width = 5, height = 4, dpi = 300)