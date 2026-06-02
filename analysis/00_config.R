## ===============================
## analysis/00_config.R
## Project-wide config (paths, seeds, filenames)

# If you ever need to run from outside the repo root, you can set:
# Sys.setenv(PROJECT_ROOT = "/full/path/to/repo")
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = ".")

# --- Core directories ---
DIR_DATA_RAW          <- file.path(PROJECT_ROOT, "data", "raw")
DIR_DATA_INTERIM      <- file.path(PROJECT_ROOT, "data", "interim")
DIR_DATA_AMPLISEQ     <- file.path(PROJECT_ROOT, "ampliseq_out")
DIR_R_OBJECTS         <- file.path(DIR_DATA_INTERIM, "r_objects")

DIR_AMPLISEQ_PHYLOSEQ <- file.path(DIR_DATA_AMPLISEQ, "phyloseq")
DIR_AMPLISEQ_META     <- file.path(DIR_DATA_AMPLISEQ, "input")
DIR_DADA2_READS       <- file.path(DIR_DATA_AMPLISEQ, "dada2")
DIR_RDS_PHYLOSEQ      <- file.path(DIR_R_OBJECTS, "phyloseq")
DIR_PICRUST_METACYC   <- file.path(DIR_DATA_AMPLISEQ, "picrust")

DIR_OUTPUTS           <- file.path(PROJECT_ROOT, "outputs")
DIR_FIG_MAIN          <- file.path(DIR_OUTPUTS, "figures", "main")
DIR_FIG_SUPP          <- file.path(DIR_OUTPUTS, "figures", "supplement")
DIR_TABLES            <- file.path(DIR_OUTPUTS, "tables")
DIR_TABLES_ML         <- file.path(DIR_TABLES, "ml")
DIR_LOGS              <- file.path(DIR_OUTPUTS, "logs")

# Create dirs (safe if they already exist)
dir.create(DIR_RDS_PHYLOSEQ, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG_MAIN,     recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG_SUPP,     recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TABLES_ML,    recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_LOGS,         recursive = TRUE, showWarnings = FALSE)

# --- Canonical filenames ---
FILE_PS_START             <- file.path(DIR_AMPLISEQ_PHYLOSEQ, "dada2_phyloseq.rds")
FILE_PS_META              <- file.path(DIR_AMPLISEQ_META, "oral_metadata.tsv")
FILE_PS_READS             <- file.path(DIR_DADA2_READS, "DADA2_stats.tsv")

FILE_PS_PROCESSED_ABS     <- file.path(DIR_RDS_PHYLOSEQ, "ps_abs.rds")
FILE_PS_PROCESSED_REL     <- file.path(DIR_RDS_PHYLOSEQ, "ps_rel.rds")
FILE_PS_PROCESSED_REL_SGN <- file.path(DIR_RDS_PHYLOSEQ, "ps_rel_sgn.rds")

FILE_PICRUST_METACYC      <- file.path(DIR_PICRUST_METACYC, "METACYC_path_abun_unstrat_descrip.tsv")

FILE_ML_IBD               <- file.path(DIR_TABLES_ML, "ml_ibd_features.csv")
FILE_ML_CD                <- file.path(DIR_TABLES_ML, "ml_cd_features.csv")
FILE_ML_UC                <- file.path(DIR_TABLES_ML, "ml_uc_features.csv")

# --- Reproducibility defaults ---
set.seed(1)
options(stringsAsFactors = FALSE)

# Ensure figure scripts stay tidy
fig_main <- function(filename) file.path(DIR_FIG_MAIN, filename)
fig_supp <- function(filename) file.path(DIR_FIG_SUPP, filename)