# Integrative Analysis of Molecular Subtypes in Acute Myeloid Leukemia

An R-based pipeline for integrated differential analysis of DNA methylation and RNA-seq data in acute myeloid leukemia (AML).

## Table of Contents
- [Project Overview](#project-overview)
  - [Features](#features)
    - [1. Preprocessing](#preprocessing)
    - [2. Differential Analysis](#differential-analysis)
    - [3. Integration](#integration)
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

## Installation
  ### Prerequisites
  ### Dependency Installation

## Usage

## Data Availability

Full data (both provided in this repository and omitted) are available here: https://utoronto-my.sharepoint.com/:f:/g/personal/sophiamjia_li_mail_utoronto_ca/ElHaP2t2nvpAi7xc0wz1zF8B1Em2vXN65rb-lLvH86eH1w?e=xNUiKx (expires June 12 2025)

### Dataset Compartments
**Configuration:** This compartment stores the configuration files for each dataset, including the default reference.
```
├── config/
│   ├── default_config.R
│   ├── ...
│   └── beataml_config.R
```
**Raw:** This compartment stores the raw data files used as input for the pipeline. It is necessary to provide methylation, expression, and metadata into this file.
```
├── raw/
│   ├── beataml_beta_values.csv
│   ├── beataml_detection_pval.csv
│   ├── beataml_clinical.xlsx
│   └── beataml_sample_mapping.xlsx
```
**Processed:** This compartment stores the processed input data outputted from the adapter and preprocessing files. They are standardized to enter the pipeline, where intermediate files are generated post-adaption, pre-processing, and post-processing. The annotation files generated for the dataset are additionally stored here. 
```
├── processed/
│   ├── beataml_adapted_methylation.rds
│   ├── beataml_adapted_expression.rds
│   ├── beataml_adapted_detection_pval.rds
│   ├── beataml_unprocessed_methylation.rds
│   ├── beataml_unprocessed_expression.rds
│   ├── beataml_processed_methylation.rds
│   ├── beataml_processed_expression.rds
│   ├── beataml_mapping.rds
│   ├── beataml_cpg_annotation.rds
│   └── beataml_promoter_annotation.rds
```
**Slurm:** This compartment stores the data and script files needed to submit the differentially methlyated regions (DMRs) using the Bumphunter function from the minfi package. Input data, R scrips, and slurm submission scripts are generated during differential analysis of the provided methylation data to be submitted to a high-performance computing cluster (HPC).
```
├── slurm/
│   ├── scripts/
│   │   ├── beataml_slurm_analysis.R
│   │   └── beataml_slurm_job.sh
│   └── data/
│       ├── beataml_slurm_data.sh
│       ├── beataml_slurm_info.sh
│       └── beataml_slurm_design.R
```
**01. Methylation:** This compartment stores the slurm compartment and results of differential methylation analysis. Intermediate and final result files are stored here in both R-native (RDS) and excel file (XLSX) formats.
```
├── 01_methylation/
│   ├── slurm/
│   |   └── ...
│   ├── beataml_raw_dmrs.rds
│   ├── beataml_significant_dmrs.rds
│   ├── beataml_dmgs.rds
│   └── beataml_dmgs.xlsx
```
**02. Expression:** This compartment stores the results of differential expression analysis. Final result files are stored here in both R-native (RDS) and excel file (XLSX) formats.
```
├── 02_expression/
│   ├── beataml_degs.rds
│   └── beataml_degs.xlsx
```
**03. Integration:** This compartment stores the results of integrative analysis of the differential methylation and expression analysis results. Correlation results of all identified genes and genes filtered for significance are stored here in both R-native (RDS) and excel file (XLSX) formats.
```
└── 03_integration
│   ├── beataml_background_genes.xlsx
│   ├── beataml_common_genes.rds
│   ├── beataml_common_genes.xlsx
│   ├── beataml_correlation.rds
│   ├── beataml_correlation_xlsx
│   ├── beataml_correlated_genes.rds
│   ├── beataml_correlated_genes.xlsx
│   └── beataml_anticorrelated.xlsx
```

## Configuration

## Methods

## Example Dataset
An example of the full pipeline on the Beat AML dataset is provided. The dataset contains DNA methylation data of 298 samples and RNA-seq data of 

## Contact

Contact Sophia Li at sophiamjia.li@utoronto.ca for questions regarding pipeline usage or data access.
