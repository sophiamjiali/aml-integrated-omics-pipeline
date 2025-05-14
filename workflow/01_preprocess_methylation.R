# ------------------------------------------------------------------------------
# Script:       01_preprocess_methylation.R
# Purpose:      Preprocesses methylation data for differential analysis
# Author:       Sophia Li
# Affiliation:  Sadreyev lab
# Date:         2025-05-12
#
#
# Version:      R version 4.4.2 (2024-10-31)
# ------------------------------------------------------------------------------

# ===| A. Remove Low Quality Probes |===========================================

# --- Remove probes missing a threshold proportion of values ---
remove_low_quality_probes <- function(methylation_cohorts, config) {
  
  message("Removing low quality probes...")
  
  # --- Compute proportion of missing values across all samples ---
  all_beta <- cbind(assays(methylation_cohorts$grp1)$beta,
                    assays(methylation_cohorts$grp2)$beta)
  
  keep_probes <- rowSums(is.na(all_beta)) < ncol(all_beta) * 
    config$methylation$preproc$low_quality_threshold

  methylation_cohorts$grp1 <- methylation_cohorts$grp1[keep_probes, ]
  methylation_cohorts$grp2 <- methylation_cohorts$grp2[keep_probes, ]
  
  message("Successfully removed low quality probes")
  return(methylation_cohorts)
}

# ===| B. Remove Probes by P-value |============================================

# --- Remove probes by P-value threshold ---
remove_probes_by_pval <- function(methylation_cohorts, config) {
  
  message("Removing probes by detection p-value...")
  
  all_pval <- cbind(assays(methylation_cohorts$grp1)$detection_pval,
                    assays(methylation_cohorts$grp2)$detection_pval)
  
  pval_threshold <- config$methylation$preproc$pval_threshold
  sample_threshold <- config$methylation$preproc$proportion_samples
  
  keep_probes <- rowSums(all_pval > pval_threshold) < ncol(all_pval) * sample_threshold
  
  methylation_cohorts$grp1 <- methylation_cohorts$grp1[keep_probes, ]
  methylation_cohorts$grp2 <- methylation_cohorts$grp2[keep_probes, ]
  
  message("Successfully removed probes by detection p-value")
  return(methylation_cohorts)
}

# ===| A. Filter Probes |=======================================================

# --- Filter out probes set for removal in the configuration file ---
filter_probes <- function(methylation_cohorts, cpg_anno) {
  
  message("Filtering probes from methylation data...")
  keep_probes <- cpg_anno$probeID
  
  methylation_cohorts$grp1 <- methylation_cohorts$grp1[
    rownames(methylation_cohorts$grp1) %in% keep_probes, 
  ]
  methylation_cohorts$grp2 <- methylation_cohorts$grp2[
    rownames(methylation_cohorts$grp2) %in% keep_probes, 
  ]
  
  message("Probes successfully filtered from methylation data")
  return (methylation_cohorts)
}

# ===| B. Normalize Beta Values (BMIQ) |========================================

# --- Normalize beta value cohorts separately via BMIQ ---
normalize_beta_values <- function(methylation_cohorts, config) {
  
  if (config$methylation$preproc$normalize_beta == TRUE) {
    
    message("Normalizing beta values by the BMIQ algorithm...")
    
    assays(methylation_cohorts$grp1)$beta <- champ.norm(
      beta = assays(methylation_cohorts$grp1)$beta,
      method = "BMIQ"
    )
    assays(methylation_cohorts$grp2)$beta <- champ.norm(
      beta = assays(methylation_cohorts$grp2)$beta,
      method = "BMIQ"
    )
    
    message("Beta values successfully normalized")
    
  } else {
    message("Beta values are already normalized")
  }
  
  return(methylation_cohorts)
}

# ===| C. Impute Missing Beta values |==========================================

# --- Impute missing beta value cohorts separately, by chromosome ---
impute_beta_values <- function(methylation_cohorts, cpg_anno, config) {
  
  message("Imputing missing beta values per cohort per chromosome...")
  
  # divide probes by chromosome
  probe_by_chr <- split(
    cpg_anno$probeID, 
    factor(seqnames(cpg_anno), levels = paste0("chr", 1:22))
  )
  probe_by_chr <- probe_by_chr[sapply(probe_by_chr, length) > 0]
  
  # initialize parallel cores
  register(MulticoreParam(workers = 4))
  
  # impute via K-nearest cores; k specified in configuration file
  impute_by_chromosome <- function(beta, probe_by_chr, k) {
    
    imputed_chrs <- bplapply(probe_by_chr, function(probes) {
      
      if (length(probes) == 0) return(NULL)
      
      chr_data <- beta[probes, , drop = FALSE]
      
      if (nrow(chr_data) > 1 && anyNA(chr_data)) {
        impute.knn(chr_data, k = k)$data
      } else {
        chr_data
      }
    })
    
    imputed_chrs <- imputed_chrs[!sapply(imputed_chrs, is.null)]
    imputed_matrix <- do.call(rbind, imputed_chrs)
    imputed_matrix[rownames(beta), , drop = FALSE]
  }
  
  assays(methylation_cohorts$grp1)$beta <- impute_by_chromosome(
    assays(methylation_cohorts$grp1)$beta, 
    probe_by_chr, 
    config$methylation$preproc$num_neighbours
  )
  assays(methylation_cohorts$grp2)$beta <- impute_by_chromosome(
    assays(methylation_cohorts$grp2)$beta, 
    probe_by_chr,
    config$methylation$preproc$num_neighbours
  )
  
  message("Missing beta values successfully imputed")
  return(methylation_cohorts)
}

# ===| D. Save Processed Results |==============================================

# --- Saves processed cohorts ---
save_processed_methylation <- function(methylation_cohorts, paths) {
  
  message("Saving processed methylation data...")
  
  processed_path <- file.path(
    paths$data_processed, 
    paste0(config$general$dataset_name, "_processed_methylation.rds")
  )
  saveRDS(methylation_cohorts, processed_path)
  message("Wrote processed methylation data to: ", processed_path)
}

