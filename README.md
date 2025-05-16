# Integrative Analysis of Molecular Subtypes in Acute Myeloid Leukemia

An R-based pipeline for integrated differential analysis of DNA methylation and RNA-seq data in acute myeloid leukemia (AML).

## Table of Contents
- [Project Overview](#project-overview)
  - [Features](#features)
    - [1. Preprocessing](#preprocessing)
    - [2. Differential Analysis](#differential-analysis)
    - [3. Integration](#integration)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Dependency Installation](#dependency-installation)
- [Usage](#usage)
- [Data Availability](#data-availability)
- [Methods](#methods)
- [Example Dataset](#example-dataset)
- [Contact]
- 
## Project Overview

This repository contains a pipeline for integrated differential analysis of DNA methylation and RNA-seq data for acute myeloid leukemia (AML) datasets. 

It assumes that DNA methylation data, detection p-values, and RNA-seq data are provided. Subsetting into cohorts require clinical data and a sample mapping file. The provided data must have non-zero sample overlap between DNA methylation and RNA-seq data availability. Any level of preprocessing is acceptable and is not required to be matched between the data. Thresholds and parameters are set in configuration files.

A list of differentially methylated genes, differentially expressed genes, and correlated genes are produced, along with intermediates at major processing points.

### Workflow

The pipeline consists of the following main steps:

**1. Configuration**

The user-specified configuration file is loaded and project directories are set up. 

**2. Data Import**

Raw data (beta values, detection p-values, gene-level counts) and metadata (clinical data, sample mapping) are loaded. Adapters standardize the raw data and metadata into compatible formats. 

**3. Preprocessing**

The standardized data is subsetted by cohort and preprocessed for quality control, filtering, imputation, and normalization as specified in the associated configuration file. 

**4. Differential Analysis**



**3. Differential Methylation**
**4. Differential Expression**
**5. Integration**

## Features

## Repository Structure
The following repository structure includes the provided BeatAML example dataset and all data provided for download.

```bash
.
├── config/
│   ├── default_config.R
│   └── beataml_config.R
|
├── templates/
│   ├── template_slurm_analysis.R
│   └── template_slurm_job.sh
│
├── utilities/
│   ├── convert_txt_to_csv.R
│   └── setup.R
│
├── adapters/
│   ├── data_adapter.R
│   ├── generate_id_mapping.R
│   └── generate_annotations.R
│
├── workflow/
│   ├── 01_preprocess_methylation.R
│   ├── 02_preprocess_expression.R
│   ├── 03_differential_methylation.R
│   ├── 04_differential_expression.R
│   └── 05_integration.R
│
├── data/
│   ├── default/
│   │   └── epic_cross_reactive_probes.csv
│   └── beataml/
│       ├── raw/
│       │   ├── beataml_beta_values.csv
│       │   ├── beataml_detection_pval.csv
│       │   ├── beataml_star_gene_counts.csv
│       │   ├── beataml_clinical.xlsx
│       │   └── beataml_sample_mapping.xlsx
│       └── processed/
│           ├── beataml_adapted_methylation.rds
│           ├── beataml_adapted_expression.rds
│           ├── beataml_adapted_detection_pval.rds
│           ├── beataml_unprocessed_methylation.rds
│           ├── beataml_unprocessed_expression.rds
│           ├── beataml_processed_methylation.rds
│           ├── beataml_processed_expression.rds
│           ├── beataml_mapping.rds
│           ├── beataml_cpg_annotation.rds
│           └── beataml_promoter_annotation.rds
│
├── results/
│   └── beataml/
│       ├── 01_methylation/
│       │   ├── slurm/
│       │   │   ├── scripts/
│       │   │   │   ├── beataml_slurm_analysis.R
│       │   │   │   └── beataml_slurm_job.sh
│       │   │   ├── data/
│       │   │   │   ├── beataml_slurm_data.sh
│       │   │   │   ├── beataml_slurm_info.sh
│       │   │   │   └── beataml_slurm_design.R
│       │   ├── beataml_raw_dmrss.rds
│       │   ├── beataml_significant_dmrs.rds
│       │   ├── beataml_dmgs.rds
│       │   └── beataml_dmgs.xlsx
│       ├── 02_expression/
│       │   ├── beataml_degs.rds
│       │   └── beataml_degs.xlsx
│       └── 03_integration
│           ├── beataml_background_genes.xlsx
│           ├── beataml_common_genes.rds
│           ├── beataml_common_genes.xlsx
│           ├── beataml_correlation.rds
│           ├── beataml_correlation_xlsx
│           ├── beataml_correlated_genes.rds
│           ├── beataml_correlated_genes.xlsx
│           └── beat_anticorrelated.xlsx
│
├── run_pipeline.R
├── aml-integrated-omics-pipeline.Rproj
└── README.md
```

## Installation
  ### Prerequisites
  ### Dependency Installation

## Usage

## Data Availability

Full data (both provided in this repository and omitted) are available here: https://utoronto-my.sharepoint.com/:f:/g/personal/sophiamjia_li_mail_utoronto_ca/ElHaP2t2nvpAi7xc0wz1zF8B1Em2vXN65rb-lLvH86eH1w?e=xNUiKx (expires June 12 2025)

## Configuration

## Methods

## Example Dataset
An example of the full pipeline on the Beat AML dataset is provided. The dataset contains DNA methylation data of 298 samples and RNA-seq data of 

## Contact

Contact Sophia Li at sophiamjia.li@utoronto.ca for questions regarding pipeline usage or data access.
