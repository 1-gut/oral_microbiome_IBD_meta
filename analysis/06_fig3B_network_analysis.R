## ===============================
## analysis/05_fig4B_maaslin2.R
## Picrust2 analysis

source("analysis/00_config.R")
source("R/packages.R")
source("R/functions_phyloseq.R")

# Import data and process -------------------------------------------------
# Read phyloseq object
raw_abs <- readRDS(FILE_PS_PROCESSED_ABS)

# Filtering data to remove taxa present in less than 5% of samples ---------
prev_prop <- 0.05      # 5% of samples
ra_prop   <- 0.0001    # 0.01% as proportion

otu <- as(otu_table(raw_abs), "matrix")
if (!taxa_are_rows(raw_abs)) otu <- t(otu)

lib <- colSums(otu)

# per-sample absolute count threshold equivalent to 0.01% RA
abs_thresh_per_sample <- ceiling(ra_prop * lib)  # length = n_samples

# For each ASV, in how many samples does it meet that sample's raw_abs threshold?
meets <- sweep(otu, 2, abs_thresh_per_sample, FUN = ">=")  # taxa x samples (logical)

prev <- rowSums(meets) / ncol(otu)

taxa_keep <- rownames(otu)[prev >= prev_prop]

ps_filtered <- prune_taxa(taxa_keep, raw_abs)

ps_filtered <- tax_glom(ps_filtered, taxrank = "Genus")

# Clean object -----------------------------------------------------------
# Store absolute abundances object
abs = ps_filtered

# Check taxonomy table columns
colnames(tax_table(abs))

# Pull out the Genus column as a character vector
genus_vec <- as.character(tax_table(abs)[, "Genus"])

# Rename taxa to Genus
taxa_names(abs) <- genus_vec

# Perform network analysis ------------------------------------------------
# IBD analysis and plotting -----------------------------------------------

# Subset phyloseq object into disease-specific groups
IBD <- subset_samples(abs, IBD_group_name == "IBD")
HC <- subset_samples(abs, IBD_group_name == "HC")

# Select phyloseq object for SPIEC-EASI analysis
ps.speic = IBD

# Run SPIEC-EASI network inference using Meinshausen-Bühlmann method
se.abs <- spiec.easi(ps.speic, method='mb', lambda.min.ratio=1e-2,
                     nlambda=20, pulsar.params=list(rep.num=20, ncores=4))

# Convert inferred adjacency matrix into igraph object
g <- adj2igraph(getRefit(se.abs))

# Extract edge association weights from SPIEC-EASI model
beta <- getOptBeta(se.abs)

# Initialise vector to store edge weights
weights <- numeric(ecount(g))

# Assign beta coefficients to corresponding network edges
for (i in seq_len(ecount(g))) {
  e <- ends(g, i)
  weights[i] <- beta[e[1], e[2]]
}

# Store edge weights within igraph object
E(g)$weight <- weights

# Remove isolated nodes with no network connections
g_filt <- delete_vertices(g, degree(g) == 0)

# Calculate node degree for each remaining vertex
deg <- degree(g_filt)

# Incorporate bacterial taxa ----------------------------------------------

# Extract taxonomy table from phyloseq object
tax <- tax_table(ps.speic)

# Convert taxonomy table to dataframe
tax <- as.data.frame(tax)

# Generate genus-level labels for network nodes
lab_map <- with(tax, ifelse(!is.na(Genus) & Genus != "",
                            gsub("^g__", "", Genus),
                            rownames(tax)))

# Match taxonomy labels to feature IDs
names(lab_map) <- rownames(tax)

# Assign bacterial genus labels to network vertices
V(g_filt)$label <- lab_map[V(g_filt)$name]

# Colour edges by association direction
# Positive associations = blue, negative associations = red
E(g_filt)$color <- ifelse(E(g_filt)$weight > 0, "blue", "red")

# Assign grey colour to edges with missing or zero weights
E(g_filt)$color[is.na(E(g_filt)$weight) | E(g_filt)$weight == 0] <- "grey"

# Use absolute edge weights for graph layout calculations
w_lay <- abs(E(g_filt)$weight)

# Replace NA or zero weights with smallest positive value
# Required for weighted Fruchterman-Reingold layout
min_pos <- min(w_lay[w_lay > 0], na.rm = TRUE)
w_lay[is.na(w_lay) | w_lay <= 0] <- min_pos

# Generate weighted network layout
lay <- layout_with_fr(g_filt, weights = w_lay)

# Open PNG graphics device for high-resolution network plot
png(
  file.path(DIR_FIG_MAIN, "fig3B_IBD_network.png"),
  width = 3000,
  height = 3000,
  res = 300
)

# Plot inferred microbial association network
plot(
  g_filt,
  layout = lay,
  vertex.label = ifelse(deg > 0, V(g_filt)$label, NA),
  vertex.size  = 10,
  vertex.label.cex = 0.7,
  vertex.label.color = "black",
  vertex.label.family = "sans",
  vertex.label.font = 2,
  edge.color  = E(g_filt)$color,
  edge.width  = scales::rescale(abs(E(g_filt)$weight), to = c(1, 8))
)

# Close graphics device and save plot
dev.off()


# HC analysis and plotting ------------------------------------------------

# Select phyloseq object for SPIEC-EASI analysis
ps.speic = HC

# Run SPIEC-EASI network inference using Meinshausen-Bühlmann method
se.abs <- spiec.easi(ps.speic, method='mb', lambda.min.ratio=1e-2,
                     nlambda=20, pulsar.params=list(rep.num=20, ncores=4))

# Convert inferred adjacency matrix into igraph object
g <- adj2igraph(getRefit(se.abs))

# Extract edge association weights from SPIEC-EASI model
beta <- getOptBeta(se.abs)

# Initialise vector to store edge weights
weights <- numeric(ecount(g))

# Assign beta coefficients to corresponding network edges
for (i in seq_len(ecount(g))) {
  e <- ends(g, i)
  weights[i] <- beta[e[1], e[2]]
}

# Store edge weights within igraph object
E(g)$weight <- weights

# Remove isolated nodes with no network connections
g_filt <- delete_vertices(g, degree(g) == 0)

# Calculate node degree for each remaining vertex
deg <- degree(g_filt)

# Incorporate bacterial taxa ----------------------------------------------

# Extract taxonomy table from phyloseq object
tax <- tax_table(ps.speic)

# Convert taxonomy table to dataframe
tax <- as.data.frame(tax)

# Generate genus-level labels for network nodes
lab_map <- with(tax, ifelse(!is.na(Genus) & Genus != "",
                            gsub("^g__", "", Genus),
                            rownames(tax)))

# Match taxonomy labels to feature IDs
names(lab_map) <- rownames(tax)

# Assign bacterial genus labels to network vertices
V(g_filt)$label <- lab_map[V(g_filt)$name]

# Colour edges by association direction
# Positive associations = blue, negative associations = red
E(g_filt)$color <- ifelse(E(g_filt)$weight > 0, "blue", "red")

# Assign grey colour to edges with missing or zero weights
E(g_filt)$color[is.na(E(g_filt)$weight) | E(g_filt)$weight == 0] <- "grey"

# Use absolute edge weights for graph layout calculations
w_lay <- abs(E(g_filt)$weight)

# Replace NA or zero weights with smallest positive value
# Required for weighted Fruchterman-Reingold layout
min_pos <- min(w_lay[w_lay > 0], na.rm = TRUE)
w_lay[is.na(w_lay) | w_lay <= 0] <- min_pos

# Generate weighted network layout
lay <- layout_with_fr(g_filt, weights = w_lay)

# Open PNG graphics device for high-resolution network plot
png(
  file.path(DIR_FIG_MAIN, "fig3B_HC_network.png"),
  width = 3000,
  height = 3000,
  res = 300
)

# Plot inferred microbial association network
plot(
  g_filt,
  layout = lay,
  vertex.label = ifelse(deg > 0, V(g_filt)$label, NA),
  vertex.size  = 10,
  vertex.label.cex = 0.7,
  vertex.label.color = "black",
  vertex.label.family = "sans",
  vertex.label.font = 2,
  edge.color  = E(g_filt)$color,
  edge.width  = scales::rescale(abs(E(g_filt)$weight), to = c(1, 8))
)

# Close graphics device and save plot
dev.off()
