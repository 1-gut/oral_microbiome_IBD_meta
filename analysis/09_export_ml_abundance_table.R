## ===============================
## analysis/07_export_ml_abundance_table.R
## Diversity analysis

source("analysis/00_config.R")
source("R/packages.R")
source("R/functions_phyloseq.R")

# Load data
metadata <- read.table(FILE_PS_META, sep = "\t", header = TRUE, row.names = 1)
rel <- readRDS(FILE_PS_PROCESSED_REL)
rel_sgn <- readRDS(FILE_PS_PROCESSED_REL_SGN)

# Export IBD abundance table for ML -----------------------------------------
# Extract metadata
metadata_df <- as.data.frame(sample_data(rel))
metadata_df <- as.matrix(metadata_df)
metadata_df <- as.data.frame(metadata_df)

# Extract OTU table
otu <- otu_table(rel)
feat <- as.data.frame(otu)
feat <- t(feat)

## Tax table
tax <- tax_table(rel)
tax <- as.data.frame(tax)

tax$feature <- rownames(tax)

# Ensure row names of df1 match column names of df2
matching_indices <- match(colnames(feat), rownames(tax))

# Replace column names in df2 with values from the "Genus" column in df1
colnames(feat)[!is.na(matching_indices)] <- tax$Genus[matching_indices[!is.na(matching_indices)]]

all_dfs <- merge(metadata_df, feat, by = "row.names")
colnames(all_dfs)[1] <- "sampleID"

write.csv(all_dfs, FILE_ML_IBD, row.names = FALSE)


# Export CD abundance table for ML -----------------------------------------
# export CD study samples
cd_studies <- c("Scotland - New", "PRJEB70783", "PRJNA773617", "SRP385133")
cd <- subset_samples(rel, study_alias %in% cd_studies)
cd <- subset_samples(cd, study_group_name == "CD" | study_group_name == "HC")

# Extract metadata
metadata_df <- as.data.frame(sample_data(cd))
metadata_df <- as.matrix(metadata_df)
metadata_df <- as.data.frame(metadata_df)

# Extract OTU table
otu <- otu_table(cd)
feat <- as.data.frame(otu)
feat <- t(feat)

## Tax table
tax <- tax_table(cd)
tax <- as.data.frame(tax)

tax$feature <- rownames(tax)

# Ensure row names of df1 match column names of df2
matching_indices <- match(colnames(feat), rownames(tax))

# Replace column names in df2 with values from the "Genus" column in df1
colnames(feat)[!is.na(matching_indices)] <- tax$Genus[matching_indices[!is.na(matching_indices)]]

all_dfs <- merge(metadata_df, feat, by = "row.names")
colnames(all_dfs)[1] <- "sampleID"

write.csv(all_dfs, FILE_ML_CD, row.names = FALSE)


# Export CD abundance table for ML -----------------------------------------
# export UC study samples
uc_studies <- c("Scotland - New", "PRJEB70783", "PRJNA684508", "PRJNA749643", "SRP385133")
uc <- subset_samples(rel, study_alias %in% uc_studies)
uc <- subset_samples(uc, study_group_name == "UC" | study_group_name == "HC")

# Extract metadata
metadata_df <- as.data.frame(sample_data(uc))
metadata_df <- as.matrix(metadata_df)
metadata_df <- as.data.frame(metadata_df)

# Extract OTU table
otu <- otu_table(uc)
feat <- as.data.frame(otu)
feat <- t(feat)

## Tax table
tax <- tax_table(uc)
tax <- as.data.frame(tax)

tax$feature <- rownames(tax)

# Ensure row names of df1 match column names of df2
matching_indices <- match(colnames(feat), rownames(tax))

# Replace column names in df2 with values from the "Genus" column in df1
colnames(feat)[!is.na(matching_indices)] <- tax$Genus[matching_indices[!is.na(matching_indices)]]

all_dfs <- merge(metadata_df, feat, by = "row.names")
colnames(all_dfs)[1] <- "sampleID"

write.csv(all_dfs, FILE_ML_UC, row.names = FALSE)
