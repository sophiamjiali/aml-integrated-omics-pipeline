# ------------------------------------------------------------------------------
# Script:       default_config.R
# Purpose:      Sets default configurations for the full integrated pipeline
# Author:       Sophia Li
# Affiliation:  Sadreyev lab
# Date:         2025-05-09
#
# Version:      R version 4.4.2 (2024-10-31)
# ------------------------------------------------------------------------------

config <- list(
  
  # --- General configurations ---
  general = list(
    dataset_name             = "default",         # Dataset name for directories
    n_cores                  = 8,                 # Parallel computation cores
    genome_build             = "hg19",            # Genome build of data
    seed                     = 123               # Seed for random generation
  ),
  
  # --- Cohort definitions ---
  cohorts = list(
    reference                = "group_1",               # Reference cohort
    id_column                = "dbgap_subject_id",      # Clinical data ID column
    meth_id                  = "labId",                 # Methylation data ID
    expr_id                  = "dbgap_rnaseq_sample",   # Expression data ID
    
    # A) Group 1
    group_1 = list(
      name                   = "Subtype 1",       # Name of cohort 1
      column                 = "column",          # Clinical data column
      pattern                = "mutation_name"    # Pattern name
    ),
    
    # B) Group 2
    group_2 = list(
      name                   = "Subtype 2",       # Name of cohort 2
      column                 = "column",          # Clinical data column
      pattern                = "mutation_name"    # Pattern name
    )
  ),
  
  # --- Methylation pipeline parameters ---
  methylation = list(
    
    assay                   = "Illumina EPIC",   # Methylation Assay
    
    # A) Preprocessing
    preproc = list(
      low_quality_threshold = 0.10,              # Missing values threshold
      pval_threshold        = 0.05,              # Detection p-value threshold
      proportion_samples    = 0.10,              # % samples for p-value threshold
      remove_non_cpg        = TRUE,              # Remove non-CpG sites
      remove_sex_chr        = TRUE,              # Remove sex chromosomes
      remove_cr_probes      = TRUE,              # Remove cross-reactive probes
      remove_snp_probes     = TRUE,              # Remove SNP-overlapped probes
      normalize_beta        = TRUE,              # Normalize beta values
      num_neighbours        = 10                 # Imputation via K-nearest neighbors
    ),
    
    # B) Differential Analysis
    diff = list(
      bumphunter = list(
        cutoff_q            = 0.99,              # DMR cutoff quantile
        max_gap             = 500,               # maximum gap b/w probes
        min_probes          = 3,                 # minimum probe count
        permutations        = 1000,              # permutations (min 1000)
      ),
      dmr_effect_threshold  = 0.10,              # Effect size threshold for sig. DMRs
      fdr_threshold         = 0.05               # FDR threshold for significant DMRs
    )
  ),
  
  expression = list(
    
    # A) Preprocessing
    preproc = list(
      remove_gene_qc        = TRUE,              # Remove low quality genes
      min_cpm               = 1,                 # Minimum CPM per gene
      remove_low_expr       = TRUE,              # Remove lowly expressed genes
      normalize_counts      = TRUE               # Normalize gene-level counts
    ),
    
    # B) Differential Analysis
    differential = list(
      logfc_threshold       = 1.0,               # LogFC threshold for significant DEGs
      fdr_threshold         = 0.05               # FDR threshold for significant DEGs
    )
  ),
  
  integration = list(
    r_threshold             = 0.3,               # Pearson threshold for anticorrelation
    fdr_threshold           = 0.05               # FDR threshold for sig. anticorrelation
  )
)
