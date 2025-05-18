# ------------------------------------------------------------------------------
# Script:       setup.R
# Purpose:      Initializes directories, packages, and dependencies
# Author:       Sophia Li
# Affiliation:  Sadreyev lab
# Date:         2025-05-09
#
# Input:        <dataset_config.R>: Configuration file for the chosen dataset
#
# Version:      R version 4.4.2 (2024-10-31)
# ------------------------------------------------------------------------------


# ===| 1. Load Libraries |======================================================

# --- Load all libraries needed for the full pipeline ---
setup_libraries <- function() {
  
  # --- Download all dependencies ---
  install_dependencies()
  
  message("[Setup] Loading required libraries...")

  # --- Load necessary libraries ---
  suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(minfi)
    library(limma)
    library(impute)
    library(ChAMP)
    library(edgeR)
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(EnsDb.Hsapiens.v75)
    library(Mfuzz)
    library(cluster)
    library(dplyr)
    library(tidyr)
    library(vroom)
    library(tibble)
    library(writexl)
    library(readxl)
    library(doRNG)
    library(here)
    library(BiocParallel)
    library(apeglm)
    library(preprocessCore)
    library(doParallel)
    library(factoextra)
    library(metap)
    library(matrixStats)
  })
  
  message("[Setup] Libraries loaded successfully.")
  
  # --- Return invisibly ---
  invisible(NULL)
}

# ===| 2. Install Dependencies |================================================

# --- Install package dependencies that are not already installed ---
install_dependencies <- function() {
  
  message("[Setup] Installing dependencies...")
  
  # --- List the required packages for the pipeline ---
  required_packages <- c(
    "SummarizedExperiment",
    "GenomicRanges",
    "GenomicFeatures",
    "minfi",
    "limma",
    "impute",
    "ChAMP",
    "edgeR",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "BSgenome.Hsapiens.UCSC.hg19",
    "EnsDb.Hsapiens.v75",
    "Mfuzz",
    "cluster",
    "dplyr",
    "tidyr",
    "vroom",
    "tibble",
    "writexl",
    "readxl",
    "doRNG",
    "here",
    "BiocParallel",
    "apeglm",
    "preprocessCore",
    "doParallel",
    "factoextra",
    "metap",
    "matrixStats"
    )
  
  # --- Install BiocManager if not already installed ---
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # --- Install all packages as necessary ---
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      
      tryCatch({
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      }, error = function(e) {
        install.packages(pkg)
      })
    }
  }
  
  message("[Setup] Dependencies successfully installed.")
  
  # --- Return invisibly ---
  invisible(NULL)
}

# ===| 3. Load all Sources |====================================================

# --- Loads all source scripts ---
load_sources <- function() {
  source(here::here("utilities/convert_txt_to_csv.R"))
  source(here::here("adapters/data_adapter.R"))
  source(here::here("adapters/generate_id_mapping.R"))
  source(here::here("adapters/generate_annotations.R"))
  source(here::here("workflow/01_preprocess_methylation.R"))
  source(here::here("workflow/02_preprocess_expression.R"))
  source(here::here("workflow/03_differential_methylation.R"))
  source(here::here("workflow/04_differential_expression.R"))
  source(here::here("workflow/05_integration.R"))
}

# ===| 4. Load Data Paths |=====================================================

# --- Set up and define file paths to main directories ---
setup_directories <- function(dataset_name) {
  
  message("[Setup] Initializing project directories for dataset: ", dataset_name)

  # --- Define base project directories accessed by the pipeline ---
  proj_dir    <- here::here()
  data_dir    <- file.path(proj_dir, "data", dataset_name)
  results_dir <- file.path(proj_dir, "results", dataset_name)
  
  # --- Ensure key directories exist ---
  paths <- list(
    data_default        = file.path(proj_dir, "data", "default"),
    data_raw            = file.path(data_dir, "raw"),
    data_processed      = file.path(data_dir, "processed"),
    templates           = file.path(proj_dir, "templates"),
    slurm_scripts       = file.path(results_dir, "01_methylation", "slurm", "scripts"),
    slurm_data          = file.path(results_dir, "01_methylation", "slurm", "data"),
    results_methylation = file.path(results_dir, "01_methylation"),
    results_expression  = file.path(results_dir, "02_expression"),
    results_integration = file.path(results_dir, "03_integration")
  )
    
  for (d in unique(paths)) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      message("[Setup] Created directory: ", d)
    }
  }
  
  message("[Setup] Directory initialization complete.")
  return(paths)
}


# ===| 5. Load Configurations |=================================================

# --- Fetch configurations for the current dataset ---
load_config <- function(config_path) {
  
  # --- Verify a valid configuration path was provided ---
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path)
  }
  
  # --- Fetch the configuration file as a source ---
  source(config_path)
  if (!exists("config")) {
    stop("Configuration file was improperly defined in ", config_path)
  }
  
  # --- Return the configuration object ---
  return(config)
}

