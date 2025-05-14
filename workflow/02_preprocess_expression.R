# ------------------------------------------------------------------------------
# Script:       02_preprocess_expression.R
# Purpose:      Preprocesses expression data for differential analysis
# Author:       Sophia Li
# Affiliation:  Sadreyev lab
# Date:         2025-05-12
#
# Version:      R version 4.4.2 (2024-10-31)
# ------------------------------------------------------------------------------

# ===| A. Filter Genes |========================================================

# --- Filter for protein-coding genes specified in the configurations file ---
filter_genes <- function(expression_cohorts, promoters) {
  
  message("Filtering genes from expression data...")
  keep_genes <- promoters$symbol

  expression_cohorts$grp1 <- expression_cohorts$grp1[
    rownames(expression_cohorts$grp1) %in% keep_genes
  ]
  expression_cohorts$grp2 <- expression_cohorts$grp2[
    rownames(expression_cohorts$grp2) %in% keep_genes
  ]
  
  message("Genes successfully filtered from expression data")
  return (expression_cohorts)
}

# ===| B. Gene Quality Control |================================================

# --- Remove low quality genes with zero or low counts across all samples ---
remove_low_quality_genes <- function(expression_cohorts, config) {
  
  if (config$expression$preproc$remove_gene_qc == TRUE) {
    
    message("Removing low quality genes...")
    
    # compute quality across all samples
    all_counts <- cbind(assays(expression_cohorts$grp1)$counts, 
                        assays(expression_cohorts$grp2)$counts)
    
    # filter for CPM threshold defined in the configuration file
    cpm_counts <- cpm(all_counts)
    
    min_samples <- min(ncol(expression_cohorts$grp1), ncol(expression_cohorts$grp2))
    keep <- rowSums(cpm_counts > 1) >= min_samples
    counts <- all_counts[keep, ]
    
    # filter original expression data
    keep_genes <- rownames(counts)
    expression_cohorts$grp1 <- expression_cohorts$grp1[keep_genes, ]
    expression_cohorts$grp2 <- expression_cohorts$grp2[keep_genes, ]
    
    message("Low quality genes successfully removed")
    
  } else {
    message("Gene quality control disabled in configurations")
  }
  
  return(expression_cohorts)
}

# ===| C. Filter Lowly Expressed Genes |========================================

remove_lowly_expressed_genes <- function(expression_cohorts, config) {
  
  if (config$expression$preproc$remove_low_expr == TRUE) {
    
    message("Removing lowly expressed genes...")
    
    # normalize counts across all samples
    all_counts <- cbind(assays(expression_cohorts$grp1)$counts, 
                        assays(expression_cohorts$grp2)$counts)
    
    # define the design matrix with the specified reference group
    group <- factor(c(rep(config$cohorts$group_1$name, ncol(expression_cohorts$grp1)),
                      rep(config$cohorts$group_2$name, ncol(expression_cohorts$grp2))))
    group <- relevel(group, ref = config$cohorts$reference)
    design_matrix <- model.matrix(~ group)
    
    # convert raw counts to DGEList and remove lowly expressed genes
    dge <- DGEList(all_counts)
    keep_genes <- filterByExpr(dge, design_matrix)
    dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
    
    # add the DGEList as a third field in the expression object
    expression_cohorts$dge <- dge
    
    message("Successfully removed lowly expressed genes")
    
  } else {
    message("Lowly expressed gene removal disabled in configurations")
  }
  
  return(expression_cohorts)
}

# ===| D. Normalize Counts |====================================================

# --- Normalize counts to TMM ---
normalize_counts <- function(expression_cohorts, config) {
  
  if (config$expression$preproc$normalize_counts == TRUE) {
    
    message("Normalizing counts to TMM...")
    
    # calculate TMM normalization factors
    expression_cohorts$dge <- calcNormFactors(expression_cohorts$dge, method = "TMM")
    
    message("Successfully normalized counts")
    
  } else {
    message("Count normalization disabled in configurations")
  }
  
  return(expression_cohorts)
  
}

# ===| E. Save Processed Results |==============================================

# --- Saves processed cohorts ---
save_processed_expression <- function(expression_cohorts, paths) {
  
  message("Saving processed expression data...")
  
  processed_path <- file.path(
    paths$data_processed, 
    paste0(config$general$dataset_name, "_processed_expression.rds")
  )
  saveRDS(expression_cohorts, processed_path)
  message("Wrote processed expression data to: ", processed_path)
}