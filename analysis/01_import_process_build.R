## ===============================
## analysis/01_import_process_build.R
## Create phyloseq object, process, and save processed objects

source("analysis/00_config.R")
source("R/packages.R")
source("R/functions_phyloseq.R")

# Import data and process -------------------------------------------------
# Read phyloseq object
data <- readRDS(FILE_PS_START)

# Read in metadata
metadata <- read.table(FILE_PS_META, sep = "\t", header = TRUE, row.names = 1)

# create phyloseq metadata object
metadata.p <- sample_data(metadata)

# Update the metadata in the phyloseq object
sample_data(data) <- metadata.p

# Read read count file
nonchim_count <- read.table(FILE_PS_READS, header = TRUE, sep = "\t")

# Identify samples with < 10,000 non-chimeric reads
low_samples_10000 <- nonchim_count$sample[nonchim_count$nonchim < 10000]

data_filt <- prune_samples(!(sample_names(data) %in% low_samples_10000), data)

# Specify phyloseq objects  ------------------------------------------------
# Store absolute abundances object
abs = data_filt

# Convert into total sum scaling (relative abundance)
rel = transform_sample_counts(abs, function(x) x / sum(x))

# Remove taxa occuring <0.1%
rel <- filter_taxa(rel, function(x) sum(x) > .001, TRUE)

# Convert (AGAIN) into total sum scaling
rel <- transform_sample_counts(rel, function(x) x / sum(x))

# collapse to genus level
rel <- tax_glom(rel, taxrank = "Genus")

# subset object to keep only IBD-subtypes + HC
sgn_rel <- subset_samples(rel, study_group_name=="CD" | study_group_name=="UC" | study_group_name=="HC")

# Save objects
saveRDS(abs, FILE_PS_PROCESSED_ABS)
saveRDS(rel, FILE_PS_PROCESSED_REL)
saveRDS(sgn_rel, FILE_PS_PROCESSED_REL_SGN)

