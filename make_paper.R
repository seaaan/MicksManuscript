# This file runs all of the scripts that analyze the data and generate the statistics and figures
# for the paper. 

# This script is the only place containing external file paths, so it is the only file that may 
# need manual changes to reproduce the analysis. (For example, after downloading the data from
# GEO, you may need to adjust the file paths to point to the right places.) This script assumes
# that your working directory is the same folder in which it is contained. 

# Load all packages and utility functions
source("source_all.R")

###################################################################################################
# Analyze raw data
###################################################################################################
# Microarray data
microarray_batch_1_file <- "data/raw/GSE142058_GenomeStudio_output_batch_1.txt.gz"
microarray_batch_1_metadata <- "data/raw/GSE142058_Metadata_batch_1.csv.gz"
microarray_batch_2_file <- "data/raw/GSE142058_GenomeStudio_output_batch_2.txt.gz"
microarray_batch_2_metadata <- "data/raw/GSE142058_Metadata_batch_2.csv.gz"
source("R/clean_microarray_data.R")
source("R/modeling_microarray.R")

# ddPCR data
ddpcr_data_folder <- "data/raw/ddpcr_Granulysin_YWHAZ_raw/"
source("R/clean_ddpcr_data.R")

# Protein data
elisa_data_folder <- "data/raw/ELISA_Granulysin_raw/"
source("R/clean_elisa_granulysin_data.R")
msd_data_folder <- "data/raw/MSD_raw/"
source("R/clean_msd_data.R")

source("R/modeling_ddpcr_protein.R")

###################################################################################################
# Analyze other studies
###################################################################################################

# GSE80455
# GSE80455_non_normalized.txt from GEO (needed for detection p-values)
GSE80455_raw_file <- "data/raw/other_studies/GSE80455/GSE80455_non_normalized.txt"
# Code itself downloads the normalized data (GSE80455_series_matrix.txt.gz) and platform
# data (GPL10558.soft) via getGEO

# GSE29981
# metadata created from the GEO2R window on GEO
GSE29981_metadata_file <- "data/raw/other_studies/GSE29981/gse29981_pData.csv"
GSE29981_cel_directory <- "data/raw/other_studies/GSE29981/"

# GSE28044
# metadata created from the GEO2R window on GEO
GSE28044_metadata_file <- "data/raw/other_studies/GSE28044/gse28044_pData.csv"
GSE28044_cel_directory <- "data/raw/other_studies/GSE28044/"

# GSE86491
GSE86491_raw_file <- "data/raw/other_studies/GSE86491/GSE86491_totalRNA_tissue.counts.txt.gz"

# GSE6364
# metadata created from the GEO2R window on GEO
GSE6364_metadata_file <- "data/raw/other_studies/GSE6364/gse6364-meta.csv"
GSE6364_cel_directory <- "data/raw/other_studies/GSE6364/"

# GSE122248
# metadata created from the supplementary data of the paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6749057/
GSE122248_raw_file <- "data/raw/other_studies/GSE122248/GSE122248_counts.matrix.txt.gz"
GSE122248_pData_file <- "data/raw/other_studies/GSE122248/GSE122248_meta.csv"

# PMID26750085
ANALYZE_PMID26750085 <- FALSE # set to TRUE if you have the PMID26750085 data
PMID26750085_124_study_cel_directory <- "data/raw/other_studies/PMID26750085/124_study/"
PMID26750085_119_study_cel_directory <- "data/raw/other_studies/PMID26750085/119_study/"

PMID26750085_124_study_metadata_file <- "data/raw/other_studies/PMID26750085/124_study/FOL_LUT_ expt description.xlsx"
PMID26750085_119_study_metadata_file <- "data/raw/other_studies/PMID26750085/119_study/Experiment Description.xlsx"

source("R/modeling_other_peoples_studies.R")

###################################################################################################
# Make figures and statistics
###################################################################################################
# Generate figures
source("R/figures.R")

# Generate supplement
source("R/generate_supplement.R")
