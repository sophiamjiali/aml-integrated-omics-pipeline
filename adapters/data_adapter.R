# ------------------------------------------------------------------------------
# Script:       default_config.R
# Purpose:      Adapter for methylation and expression data for the pipeline
# Author:       Sophia Li
# Affiliation:  Sadreyev lab
# Date:         2025-05-11
#
# Inputs:       
#   - beta_values_path:                Path to a CSV of β values
#   - gene_counts_path:                Path to a CSV of gene-level counts 
#   - dataset_name:                    Name of the current dataset
#   - paths:                           List of primary directories
#
# Outputs:   
#   - <dataset>_methylation_data.rds:  SummarizedExperiment of β values
#   - <dataset>_expression_data.rds:   SummarizedExperiment of gene-level counts
#
# Version:      R version 4.4.2 (2024-10-31)
# ------------------------------------------------------------------------------

# ===| A. Process Methylation Data |============================================

run_methylation_adapter <- function(data_paths, dataset_name, paths) {
  
  # --- Read the CSV file into a matrix of beta values---
  beta_values <- vroom(
    file = data_paths[1],
    delim = ",",
    col_names = TRUE,
    show_col_types = FALSE
  )
  beta_values <- column_to_rownames(beta_values, var = colnames(beta_values)[1])
  beta_values <- as.matrix(beta_values)
  
  # --- Sanity checks for data validity ---
  
  # 1. non‑empty  
  if (nrow(beta_values) == 0 || ncol(beta_values) == 0) {
    stop("Beta values are empty: ", beta_values_path)
  }
  
  # 2. numeric  
  if (!is.numeric(beta_values)) {
    stop("Beta values contain non‑numeric values")
  }
  
  # 3. range [0,1]  
  if (any(beta_values < 0 | beta_values > 1, na.rm = TRUE)) {
    stop("Beta values contain values outside [0,1]")
  }
  
  # 4. missing values  
  na_frac <- mean(is.na(beta_values))
  if (na_frac > 0.01) {
    stop(sprintf("Beta values contain too many missing values: %.1f%%", na_frac * 100))
  } else if (na_frac > 0) {
    message(sprintf("Warning: %.1f%% missing beta values; imputation will be performed", na_frac * 100))
  }
  
  # 5. unique names  
  if (any(duplicated(rownames(beta_values)))) {
    stop("Beta values contain duplicate probe IDs")
  }
  if (any(duplicated(colnames(beta_values)))) {
    stop("Beta values contain duplicate sample names")
  }
  
  # --- Process the data into a SummarizedExperiment ---
  methylation_data <- SummarizedExperiment(
    assays = SimpleList(beta = beta_values),
    rowData = DataFrame(probeID = rownames(beta_values)),
    colData = DataFrame(sample_id = colnames(beta_values)),
    metadata = list(
      reference_genome = "hg19",
      dataset_name = dataset_name,
      date_processed = Sys.Date())
  )
  
  # --- Save and return the processed data ---
  processed_path <- file.path(
    paths$data_processed, 
    paste0(dataset_name, "_adapted_methylation.rds")
  )
  saveRDS(methylation_data, processed_path)
  message("Wrote adapted methylation data to: ", processed_path)
  
  return(methylation_data)
}

# ===| B. Process Detection P-values |==========================================

run_pval_adapter <- function(data_paths, dataset_name, paths) {
  
  # --- Read the CSV file into a matrix of detection P-values ---
  detection_pval <- vroom(
    file = data_paths[2],
    delim = ",",
    col_names = TRUE,
    show_col_types = FALSE
  )
  detection_pval <- column_to_rownames(detection_pval, var = colnames(detection_pval)[1])
  detection_pval <- as.matrix(detection_pval)
  
  # --- Save and return the processed data ---
  processed_path <- file.path(
    paths$data_processed, 
    paste0(dataset_name, "_adapted_detection_pval.rds")
  )
  saveRDS(detection_pval, processed_path)
  message("Wrote adapted detection p-values to: ", processed_path)
  
  return(detection_pval)
}


# ===| C. Process Expression Data |=============================================

run_expression_adapter <- function(data_paths, dataset_name, paths) {
  
  # --- Read the CSV file into a matrix of gene-level counts ---
  gene_counts <- vroom(
    file = data_paths[3],
    delim = ",",
    col_names = TRUE,
    show_col_types = FALSE
  )
  gene_counts <- column_to_rownames(gene_counts, var = colnames(gene_counts)[1])
  
  # --- Collapse duplicate symbols to the first occurence ---
  if (anyDuplicated(gene_counts$display_label)) {
    warning("Gene counts contain duplicate gene symbols; collapsing to first occurence")
    gene_counts <- gene_counts[!duplicated(gene_counts$display_label), , drop = FALSE]
  }

  # --- Extract counts from the metadata ---
  counts <- as.matrix(gene_counts[, 4:ncol(gene_counts)])
  rownames(counts) <- gene_counts$display_label

  # --- Sanity checks for data validity ---

  # 1. non‑empty
  if (nrow(counts) == 0 || ncol(counts) == 0) {
    stop("Gene counts are empty: ", gene_counts_path)
  }
  
  # 2. integer, non‑negative  
  if (!all(counts == round(counts), na.rm = TRUE)) {
    stop("Gene counts contain non‑integer values")
  }
  if (any(counts < 0, na.rm = TRUE)) {
    stop("Gene counts contain negative counts")
  }
  
  # 3. missing values  
  na_frac <- mean(is.na(counts))
  if (na_frac > 0) {
    if (na_frac > 0.01) {
      stop(sprintf("Gene counts contain too many missing values: %.1f%%", na_frac * 100))  
    } else {
      message(sprintf("Warning: %.1f%% missing counts; imputation will be performed", na_frac * 100))
    }
  }                          
  
  # 6. unique names  
  if (any(duplicated(colnames(counts)))) {
    stop("Gene counts contain duplicate sample names")
  }
  
  # --- Process the data into a SummarizedExperiment ---
  row_data = DataFrame(
    gene_id = rownames(gene_counts),
    gene_symbol = rownames(counts),
    biotype = gene_counts$biotype,
    row.names = rownames(counts)
  )
  col_data <- DataFrame(
    sample_id = colnames(counts),
    row.names = colnames(counts)
  )
  metadata <- list(
    reference_genome = "hg19",
    dataset_name = dataset_name,
    date_processed = Sys.Date()
  )
  
  expression_data <- SummarizedExperiment(
    assays = SimpleList(counts = counts),
    rowData = row_data,
    colData = col_data,
    metadata = metadata
  )
  
  # --- Save and return the processed data ---
  processed_path <- file.path(
    paths$data_processed, 
    paste0(dataset_name, "_adapted_expression.rds")
  )
  saveRDS(expression_data, processed_path)
  message("Wrote adapted expression data to: ", processed_path)
  
  return(methylation_data)
}




