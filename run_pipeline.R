#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script:       run_pipeline.R
# Purpose:      Runs the integrated differential pipeline on the provided dataset
# Author:       Sophia Li
# Affiliation:  Sadreyev lab
# Date:         2025-05-09
#
# Inputs:       
#   - <dataset>_config.R:                   Pipeline configuration file
#   - <dataset>_beta_values.csv:            Raw beta values
#   - <dataset>_detection_pval.csv:         Raw detection P-values
#   - <dataset>_star_gene_counts.csv:       Raw STAR-aligned gene-level counts
#   - <dataset>_clinical.xlsx:              Clinical Data
#   - <dataset>_sample_mapping.xlsx:        Sample ID mapping
#
# Outputs:   
#   - <dataset>_dmgs.rds/xlsx:              Differentially methylated genes (DMGs)
#   - <dataset>_degs.rds/xlsx:              Differentially expressed genes (DEGs)
#   - <dataset>_correlated_genes.rds/xlsx:  Correlated genes (between DMG & DEG)
#
# Version:      R version 4.4.2 (2024-10-31)
# ------------------------------------------------------------------------------

# ===| 1. Load Dependencies |===================================================

# --- Load all setup, adapter, and workflow functions ---
source(here::here("utilities/setup.R"))
load_sources()

# --- Load all necessary library dependencies ---
setup_libraries()

# ===| 2. Load Configurations |=================================================

# --- Parse and validate that a valid configuration was provided ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript run_pipeline.R <dataset_config.R>\n",
       " <dataset_config.R> must be a path to your dataset-specific configuration file.")
}

# --- Load the configuration file ---
config_path <- file.path(, args[[1]])
config <- load_config(config_path)

# --- Set up key directories associated with the given dataset ---
dataset_name <- config$general$dataset_name
paths <- setup_directories(dataset_name)

# ===| 3. Adapt Data |==========================================================

# --- Parse and validate the provided raw data ---
raw_data_paths <- c(
  file.path(paths$data_raw, paste0(dataset_name, "_beta_values.csv")),
  file.path(paths$data_raw, paste0(dataset_name, "_detection_pval.csv")),
  file.path(paths$data_raw, paste0(dataset_name, "_star_gene_counts.csv"))
)

missing_files <- raw_data_paths[!file.exists(raw_data_paths)]
if (length(missing_files)) {
  stop("Required input files missing:\n ", paste(missing_files, collapse = "\n "))
} else {
  message("Raw data successfully found for dataset: ", dataset_name)
}

# --- Process the data using the adapter ---
methylation_data <- run_methylation_adapter(raw_data_paths, dataset_name, paths)
assays(methylation_data)$detection_pval <- run_pval_adapter(raw_data_paths, dataset_name, paths)
expression_data <- run_expression_adapter(raw_data_paths, dataset_name, paths)

# --- Parse and validate the provided metadata ---
metadata_paths <- c(
  file.path(paths$data_raw, paste0(config$general$dataset_name, "_clinical.xlsx")),
  file.path(paths$data_raw, paste0(config$general$dataset_name, "_sample_mapping.xlsx"))
)

missing_files <- metadata_paths[!file.exists(metadata_paths)]
if (length(missing_files)) {
  stop("Required input files missing:\n ", paste(missing_files, collapse = "\n "))
} else {
  message("Metadata successfully found for dataset: ", dataset_name)
}

# --- Generate ID mappings using the adapter ---
data_samples <- list(
  beta_ids = gsub("\\.", "-", sub("^X", "", colnames(methylation_data))),
  counts_ids = colnames(expression_data)
)
mapping <- generate_mapping(config, metadata_paths, data_samples)


# --- Clean intermediates ---
rm(args, config_path, raw_data_paths, metadata_paths, missing_files, data_samples)

# ===| 4. Run Pipeline |========================================================

# --- A. Subset the data by cohort ---
methylation_cohorts <- subset_by_cohort(methylation_data, "methylation", mapping, config)
expression_cohorts <- subset_by_cohort(expression_data, "expression", mapping, config)

# --- C. Preprocess the data ---

# process the adapted methylation data
methylation_cohorts <- remove_low_quality_probes(methylation_cohorts, config)
methylation_cohorts <- remove_probes_by_pval(methylation_cohorts, config)
cpg_anno <- generate_probe_annotation(methylation_cohorts, config)

methylation_cohorts <- filter_probes(methylation_cohorts, cpg_anno)
methylation_cohorts <- normalize_beta_values(methylation_cohorts, config)
methylation_cohorts <- impute_beta_values(methylation_cohorts, cpg_anno, config)
save_processed_methylation(methylation_cohorts, paths)

# process the adapted expression data
promoters <- generate_promoters(config)
expression_cohorts <- filter_genes(expression_cohorts, promoters)
expression_cohorts <- remove_low_quality_genes(expression_cohorts, config)
expression_cohorts <- remove_lowly_expressed_genes(expression_cohorts, config)
expression_cohorts <- normalize_counts(expression_cohorts, config)
save_processed_expression(expression_cohorts, paths)

# --- D. Differential analysis ---

# generate differentially methylated genes (DMGs)
generate_slurm_files(methylation_cohorts, cpg_anno, promoters, config)
raw_dmrs <- fetch_bumphunter_results(paths)
significant_dmrs <- filter_dmrs(raw_dmrs, paths, config)
dmgs <- identify_dmgs(significant_dmrs, promoters, paths, config)

dmgs <- split_dmgs(dmgs, paths, config)

# generate differentially expressed genes (DEGs)
expression_differentials <- generate_differentials(expression_cohorts)
degs <- identify_degs(expression_differentials, paths, config)

degs <- split_degs(degs, paths, config)

# generate background genes for external pathway enrichment analysis
generate_background_genes(cpg_anno, promoters, paths, config)

# --- E. Integrated analysis ---
common_genes <- find_common_genes(dmgs, degs, paths, config)
common_methylation <- aggregate_methylation(common_genes, methylation_cohort, 
                                            promoters, cpg_anno)
common_expression <- process_common_expression(common_genes, expression_cohort)

correlation <- calculate_pearson_correlation(common_methylation, common_expression,
                                             mapping, paths, config)
correlated_genes <- identify_correlated_genes(correlation, paths, config)

# ===| 5. Clean Up |============================================================

# --- remove intermediate objects ---
rm(config, dataset_name, paths, methylation_data, expression_data, mapping,
   methylation_cohorts, expression_cohorts, cpg_anno, promoters,
   raw_dmrs, significant_dmrs, expression_differentials, common_genes, 
   common_methylation, common_expression, correlation)
gc()




