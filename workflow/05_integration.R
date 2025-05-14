# ------------------------------------------------------------------------------
# Script:       05_integration_expression.R
# Purpose:      Integrates DMG and DEG results
# Author:       Sophia Li
# Affiliation:  Sadreyev lab
# Date:         2025-05-13
#
# Version:      R version 4.4.2 (2024-10-31)
# ------------------------------------------------------------------------------

# ===| A. Find Overlapping DMGs and DEGs |======================================

find_common_genes <- function(dmgs, degs, paths, config) {
  
  # --- Identify the intersection of genes ---
  intersection <- dmgs$all$symbol[dmgs$all$symbol %in% degs$all$symbols]
  
  common_genes <- list(
    dmgs = dmgs$all[dmgs$all$symbol %in% intersection, ],
    degs = degs$all[degs$all$symbol %in% intersection, ]
  )
  
  # --- Save and return the processed CpG probe annotation ---
  processed_path <- file.path(
    paths$results_integration, 
    paste0(config$general$dataset_name, "_common_genes.rds")
  )
  saveRDS(common_genes, processed_path)
  message("Wrote the common set of DMGs and DEGs to: ", processed_path)
  processed_path <- file.path(
    paths$results_integration, 
    paste0(config$general$dataset_name, "_common_genes.xlsx")
  )
  write_xlsx(common_genes, processed_path)
  
  return (common_genes)
}

# ===| B. Aggregate Gene-Level Methylation |=================================================

# --- Process the methylation data of common genes for correlation analysis ---
aggregate_methylation <- function(common_genes, methylation_cohort, promoters, cpg_anno) {
  
  message("Processing methylation data of the common genes...")
  
  # --- Calculate gene-level methylation of all samples ---
  all_beta <- cbind(assays(methylation_cohorts$grp1)$beta, 
                    assays(methylation_cohorts$grp2)$beta)
  
  # --- Map CpGs to overlapping promoters ---
  common_promoters <- promoters[promoters$symbol %in% common_genes$dmgs$symbol, ]
  overlaps <- findOverlaps(cpg_anno, common_promoters, ignore.strand = TRUE)
  cpg_to_gene <- data.frame(
    probeID = cpg_anno$probeID[queryHits(overlaps)],
    gene = common_promoters$symbol[subjectHits(overlaps)]
  )
  
  # --- Map beta values to genes ---
  all_beta <- all_beta[rownames(all_beta) %in% cpg_to_gene$probeID, , drop = FALSE]
  beta_with_gene <- cbind(
    gene = cpg_to_gene$gene[match(rownames(all_beta), cpg_to_gene$probeID)],
    all_beta
  )
  
  beta_with_gene <- as.data.frame(beta_with_gene)
  beta_with_gene[, -1] <- lapply(beta_with_gene[, -1], as.numeric) 
  
  # --- Aggregate beta values per gene ---
  gene_beta <- beta_with_gene %>%
    group_by(gene) %>%
    summarize(across(everything(), function(x) mean(x, na.rm = TRUE)))
  gene_beta <- as.matrix(gene_beta)
  
  message("Successfully aggregated gene-level methylation")
  return (gene_beta)
}

# ===| B. Process Common Expression  |==========================================

# --- Normalizes expression of the common genes to logCPM ---
process_common_expression <- function(common_genes, expression_cohorts) {
  
  message("Processing expression of the common genes...")
  
  # --- Subset for the common genes ---
  common_dge <- expression_cohorts$dge[rownames(expression_cohorts$dge) %in% 
                                         common_genes$dmgs$symbol, ]
  
  # --- Calculate logCPM ---
  logCPM <- cpm(common_dge, log = TRUE)
  
  message("Successfully computed logCPM from the expression of common genes")
  return(logCPM)
}

# ===| C. Run Correlation Analysis |============================================

# --- Calculate Pearson's Correlation Coefficient ---
calculate_pearson_correlation <- function(common_methylation, common_expression,
                                          mapping, paths, config) {
  
  # --- Map expression data sample IDs to labIDs if available ---
  expr_id <- config$cohorts$expr_id
  meth_id <- config$cohorts$meth_id
  
  valid_mapping <- mapping[!is.na(mapping[[meth_id]]) & mapping[[expr_id]] 
                           %in% colnames(common_expression), ]
  expr_to_meth <- setNames(paste0("X", gsub("-", ".", valid_mapping[[meth_id]])), 
                           valid_mapping[[expr_id]])
  
  # subset the expression data and reassign colnames to methylation IDs
  expression <- common_expression[, names(expr_to_meth), drop = FALSE]
  colnames(expression) <- expr_to_meth
  
  # assign methylation symbols to rownames
  methylation <- common_methylation[, -1]
  rownames(methylation) <- common_methylation[, 1]
  
  # --- Subset the input data for the common set of samples ---
  common_samples <- intersect(colnames(methylation), colnames(expression))
  methylation <- methylation[, colnames(methylation) %in% common_samples]
  expression <- expression[, colnames(expression) %in% common_samples]
  
  # align rows and columns in case order was lost
  sample_order <- sort(intersect(colnames(expression), colnames(methylation)))
  expression  <- expression[, sample_order, drop = FALSE]
  methylation <- methylation[, sample_order, drop = FALSE]
  
  gene_order <- sort(intersect(rownames(expression), rownames(methylation)))
  expression  <- expression[gene_order, , drop = FALSE]
  methylation <- methylation[gene_order, , drop = FALSE]
  
  # calculate pearson correlation and p-value for each gene
  cor_results <- t(sapply(1:nrow(methylation), function(i) {
    ct <- cor.test(as.numeric(methylation[i, ]), as.numeric(expression[i, ]), method = "pearson")
    c(correlation = ct$estimate, p.value = ct$p.value)
  }))
  rownames(cor_results) <- rownames(methylation)
  
  # adjust P-values to FDR
  cor_results <- as.data.frame(cor_results)
  cor_results$FDR <- p.adjust(cor_results$p.value, method = "BH")
  
  # --- Assemble final results ---
  cor_results <- data.frame(
    gene = rownames(cor_results),
    r = cor_results$correlation.cor,
    p.value = cor_results$p.value,
    FDR = cor_results$FDR,
    row.names = NULL
  )
  
  # --- Save processed objects ---
  processed_path <- file.path(
    paths$results_integration,
    paste0(config$general$dataset_name, "_correlation.rds")
  )
  saveRDS(cor_results, processed_path)
  message("Correlation results saved to: ", processed_path)
  processed_path <- file.path(
    paths$results_integration,
    paste0(config$general$dataset_name, "_correlation.xlsx")
  )
  saveRDS(cor_results, processed_path)
  
  return(cor_results)
}


# ===| D. Identify Correlated Genes  |==========================================

# --- Identify correlated genes by threshold filtering ---
identify_correlated_genes <- function(correlation, paths, config) {
  
  message("Identifying correlated genes...")
  
  correlated_genes <- subset(correlation,
                             abs(r) > config$integration$r_threshold &
                               FDR < config$integration$fdr_threshold)
  
  processed_path <- file.path(
    paths$results_integration,
    paste0(config$general$dataset_name, "_correlated_genes.rds")
  )
  saveRDS(correlated_genes, processed_path)
  message("Correlated genes written to: ", processed_path)
  processed_path <- file.path(
    paths$results_integration,
    paste0(config$general$dataset_name, "_correlated_genes.xlsx")
  )
  write_xlsx(correlated_genes, processed_path)
  
  return (correlated_genes)
}
