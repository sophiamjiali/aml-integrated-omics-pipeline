# ------------------------------------------------------------------------------
# Script:       04_differential_expression.R
# Purpose:      Generates differentials for processed expression data
# Author:       Sophia Li
# Affiliation:  Sadreyev lab
# Date:         2025-05-13
#
# Version:      R version 4.4.2 (2024-10-31)
# ------------------------------------------------------------------------------

# ===| A. Generate Differentials |==============================================

# --- Generate differentials for expression data ---
generate_differentials <- function(expression_cohorts) {
  
  message("Generating differentials...")
  
  # -- Create a design matrix ---
  group <- factor(c(rep(config$cohorts$group_1$name, ncol(expression_cohorts$grp1)),
                    rep(config$cohorts$group_2$name, ncol(expression_cohorts$grp2))))
  group <- relevel(group, ref = config$cohorts$reference)
  design_matrix <- model.matrix(~ group)
  
  # --- Estimate dispersion for technical/biological variability ---
  dge <- estimateDisp(expression_cohorts$dge, design_matrix)
  
  # --- Fit the Quasi-Likelihood model ---
  fit <- glmQLFit(dge, design_matrix)
  
  # --- Test for differential expression
  qlf <- glmQLFTest(fit, coef = 2)
  results <- topTags(qlf, n = Inf)$table
  
  message("Successfully generated differentials")
  return(results)
}

# ===| B. Identify Differentially Expressed Genes |=============================

# --- Identify DEGs by significance thresholds ---
identify_degs <- function(expression_differentials, paths, config) {
  
  message("Identifying differentially expressed genes...")
  
  # --- Filter differentials for significance ---
  degs <- subset(expression_differentials,
                 abs(logFC) > config$expression$differential$logfc_threshold &
                 FDR < config$expression$differential$fdr_threshold)
  
  degs <- cbind(symbols = rownames(degs), degs)
  rownames(degs) <- NULL
  
  # --- Save and return the processed CpG probe annotation ---
  processed_path <- file.path(
    paths$results_expression, 
    paste0(config$general$dataset_name, "_degs.rds")
  )
  saveRDS(degs, processed_path)
  message("Wrote DEGs to: ", processed_path)
  processed_path <- file.path(
    paths$results_expression, 
    paste0(config$general$dataset_name, "_degs.xlsx")
  )
  write_xlsx(degs, processed_path)
  
  return(degs)
}

# ===| C. Split DEGs into Up/Downregulated |====================================

split_degs <- function(degs, paths, config) {
  
  message("Splitting DEGs into up/downregulated...")
  
  split_degs <- list(
    all = degs,
    upregulated = degs[degs$logFC > 0, ],
    downregulated = degs[degs$logFC < 0, ]
  )
  
  # --- Save and return the processed CpG probe annotation ---
  processed_path <- file.path(
    paths$results_expression, 
    paste0(config$general$dataset_name, "_degs.rds")
  )
  saveRDS(split_degs, processed_path)
  message("Wrote split DEGs to: ", processed_path)
  
  # save as a multi-sheet excel file
  processed_path <- file.path(
    paths$results_expression, 
    paste0(config$general$dataset_name, "_degs.xlsx")
  )
  write_xlsx(split_degs, processed_path)
  
  return(split_degs)
}



