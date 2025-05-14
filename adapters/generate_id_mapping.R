# ------------------------------------------------------------------------------
# Script:       generate_id_mapping.R
# Purpose:      Generates mappings for methylation/expression data and cohorts
# Author:       Sophia Li
# Affiliation:  Sadreyev lab
# Date:         2025-05-07
#
# Version:      R version 4.4.2 (2024-10-31)
# ------------------------------------------------------------------------------

# ===| A. Generate Mapping |====================================================

# --- Generates a mapping for DNA methylation data, RNA-seq data, and Cohort ---
generate_mapping <- function(config, metadata_paths, data_samples) {
  
  # --- Load metadata ---
  clinical <- read_xlsx(metadata_paths[1])
  sample_mapping <- read_xlsx(metadata_paths[2])
  
  # --- Extract the sample IDs for both cohorts ---
  id_column <- config$cohorts$id_column
  grp1 <- config$cohorts$group_1
  grp2 <- config$cohorts$group_2
  
  grp1_ids <- clinical[[id_column]][grepl(
      grp1$pattern, clinical[[grp1$column]], ignore.case = TRUE, perl= TRUE
    )]
  
  grp2_ids <- clinical[[id_column]][grepl(
    grp2$pattern, clinical[[grp2$column]], ignore.case = TRUE, perl= TRUE
  )]
  
  # --- Sanity checks for valid data ---
  if (length(grp1_ids) == 0) {
    stop("No samples were found for the group: ", grp1$name)
  }
  if (length(grp2_ids) == 0) {
    stop("No samples were found for the group: ", grp2$name)
  }
  
  # --- Subset the sample mapping for cohort IDs ---
  mapping <- sample_mapping[, c(2, 3, 4), drop = FALSE]
  mapping <- mapping[mapping[[id_column]] %in% c(grp1_ids, grp2_ids), ]
  
    
  # --- Filter for IDs with either methylation or expression data available ---
  meth_id <- config$cohorts$meth_id
  expr_id <- config$cohorts$expr_id
  
  # remove any IDs that don't exist in the provided data
  mapping[[meth_id]][!mapping[[meth_id]] %in% data_samples$beta_ids] <- NA_character_
  mapping[[expr_id]][!mapping[[expr_id]] %in% data_samples$counts_ids] <- NA_character_
  
  # remove rows that do not denote any valid corresponding data
  mapping <- mapping[!(is.na(mapping[[meth_id]]) & is.na(mapping[[expr_id]])), , drop = FALSE ]
  
  # --- Set cohort representations in the mapping ---
  mapping$cohort <- NA_character_
  mapping$cohort[mapping[[id_column]] %in% grp1_ids] <- grp1$name
  mapping$cohort[mapping[[id_column]] %in% grp2_ids] <- grp2$name

  # --- Save and return the processed data ---
  processed_path <- file.path(
    paths$data_processed, 
    paste0(config$general$dataset_name, "_mapping.rds")
  )
  saveRDS(mapping, processed_path)
  message("Wrote mapping to: ", processed_path)
  
  return(mapping)
}

# ===| B. Subset Cohorts |======================================================

# --- Subsets the given data into the defined cohorts ---
subset_by_cohort <- function(data, type, mapping, config) {
  
  message("Subsetting by cohort...")
  
  # --- Preprocess the IDs to align with the data representation
  if (type == "methylation") {
    
    meth_id <- config$cohorts$meth_id
    valid_mapping <- mapping[!is.na(mapping[[meth_id]]), ]
    grp1_ids <- valid_mapping[[meth_id]][valid_mapping$cohort == config$cohorts$group_1$name]
    grp1_ids <- paste0("X", gsub("-", ".", grp1_ids))
    grp2_ids <- valid_mapping[[meth_id]][valid_mapping$cohort == config$cohorts$group_2$name]
    grp2_ids <- paste0("X", gsub("-", ".", grp2_ids))
    
  } else {
    
    expr_id <- config$cohorts$expr_id
    grp1_ids <- mapping[[expr_id]][mapping$cohort == config$cohorts$group_1$name]
    grp2_ids <- mapping[[expr_id]][mapping$cohort == config$cohorts$group_2$name]
  }
  
  # --- Subset the data by cohort ---
  grp1_idx <- match(grp1_ids, colnames(data))
  grp1 <- data[, grp1_idx[!is.na(grp1_idx)]]
  grp2_idx <- match(grp2_ids, colnames(data))
  grp2 <- data[, grp2_idx[!is.na(grp2_idx)]]
  
  # --- Save and return the processed data ---
  cohorts <- list(grp1 = grp1, grp2 = grp2)
  
  processed_path <- file.path(
    paths$data_processed, 
    paste0(config$general$dataset_name, "_unprocessed_", type, ".rds")
  )
  saveRDS(cohorts, processed_path)
  message("Subsetted data saved to: ", processed_path)
  
  return(cohorts)
}

