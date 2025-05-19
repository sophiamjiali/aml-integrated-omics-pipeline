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
  - [Input Data Format](#input-data-format)
- [Data Availability](#data-availability)
- [File Descriptions](#file-descriptions)
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

**1. Configuration:** The user-specified configuration file is loaded and project directories are set up. The format must follow the default configuration file provided as example, following EWAS standards for significance.

**2. Data Import:** Raw data (beta values, detection p-values, gene-level counts) and metadata (clinical data, sample mapping) are loaded. Adapters standardize the raw data and metadata into compatible formats. 

**3. Preprocessing:** The standardized data is subsetted by cohort and preprocessed for quality control, filtering, imputation, and normalization as specified in the associated configuration file. 

**4. Differential Analysis:**



**3. Differential Methylation**
**4. Differential Expression**
**5. Integration**

## Features

## Installation
  ### Prerequisites
  ### Dependency Installation

## Usage

This section provides step-by-step instructions for running the integrated analysis pipeline on your own acute myeloid leukemia (AML) datasets. Refer to the [Beat AML example dataset](#example-dataset) for a full run-through of the pipeline, including all expected input and output files. In the following instructions, replace `<dataset_name` with your actual dataset name (e.g. `beataml`).

### 1. Prepare Input Data

Ensure your DNA methylation, detection p-values, RNA-seq, clinical data, and sample mapping are formatted as described in the [Input Data Format](#input-data-format) section. Create a directory `data/<dataset_name>/raw/` and place all input data inside.

### 2. Prepare Configuration File

Copy the default configuration file (`config/default_config.R`) to a new file (`config/<dataset_name>_config.R`). Edit this file to specify metadata and analysis parameters. Note that the default configuration file provides EWAS standards for statistical thresholds, quality control, normalization, and reporting.

### 3. Run the Pipeline

  **A. Running in RStudio:** Open R/RStudio in the project directory (`aml-integrated-omics-pipeline.Rproj`). Open
  `run_pipeline.R` and provide the local path from the folder directory to the dataset's configuration file as `config_path <- "config/<dataset_name>_config.R"`. You may now run the pipeline line-by-line, omitting lines 36-44.
  
  **B. Running in bash:** Open a command-line program (e.g. terminal) and navigate to the project directory (`cd <path_to_directory>/aml-integrated-omics-pipeline/`). Execute the script with the dataset's configuration file provided as the first command-line argument (`./run_pipeline.R <dataset_name>_config.R`).

### 4. Review Output

Observe status messages sent by the pipeline for any errors encountered during pipeline execution that may indicate failure (e.g. "no input data provided"). After the pipeline completes successfully, results will be saved in the output directory specified (`results/<dataset_name>/`). You will find:

- Differentially methylated genes
- Differentially expressed genes
- Integrated results
- Intermediate files for quality control and troubleshooting

Refer to [Dataset Compartments](#dataset-compartments) for the general layout of what data is generated and where it is stored. Refer to [File Descriptions](#file-descriptions) for full details of each file.

  ### Input Data Format

  For beta values, detection p-values, and star gene counts, the data is expected to provide rownames as the first data column, where subsequent columns are data points with sample ID as the column name. Column names do not need to be standardized in the raw input data (e.g. `methylationID` in the sample mapping), they are to be provided in the configuration file.  

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

## File Descriptions

Refer to [Dataset Compartments](#dataset-compartments) for the detailed repository layout. Let `<dsn>` represent `<dataset_name>` in the following tables. Refer to [Data Availability](#data-availability) for example data.

### Input Files

**A. Raw Files:** This compartment stores raw data expected to be provided as input by the user. 
| File Name                     | Description                                 | Required Columns              | Format |
|-------------------------------|---------------------------------------------|-------------------------------|--------|
| `<dsn>_beta_values.csv`       | DNA methylation beta values                 | ProbeID, SampleIDs...         | CSV    |
| `<dsn>_detection_pval.csv`    | Detection p-values for methylation probes   | ProbeID, SampleIDs...         | CSV    |
| `<dsn>_star_gene_counts.csv`  | Gene-level RNA-seq counts                   | Gene, SampleIDs...            | CSV    |
| `<dsn>_clinical.xlsx`         | Clinical metadata for each sample           | SampleID, ...                 | XLSX   |
| `<dsn>_sample_mapping.xlsx`   | Mapping between sample IDs across datasets  | MethylationID, RNAseqID, ...  | XLSX   |

Refer to [Input Data Format](#input-data-format) for specific requirements of the raw input data.

**B. Processed Files:** This compartment stores processed data at various intermediate stages of adaption and preprocessing.
| File Name                           | Description                              | Format |
|-------------------------------------|------------------------------------------|--------|
| `<dsn>_adapted_methylation.rds`     | Adapted beta values                      | RDS    |
| `<dsn>_adapted_detection_pval.rds`  | Adapted detection p-values               | RDS    |
| `<dsn>_adapted_expression.csv`      | Adapted gene-level counts                | RDS    |
| `<dsn>_unprocessed_methylation.rds` | Methylation data split by cohort         | RDS   |
| `<dsn>_unprocessed_expression.rds`  | Expression data split by cohort          | RDS   |
| `<dsn>_processed_methylation.rds`   | Preprocessed methylation split by cohort | RDS   |
| `<dsn>_processed_expression.rds`    | Preprocessed Expression split by cohort  | RDS   |
| `<dsn>_mapping.rds`                 | Adapted sample mapping                   | RDS   |
| `<dsn>_cpg_annotation.rds`          | CpG probe annotations for available data | RDS   |
| `<dsn>_promoter_annotation.rds`     | Promoter annotations for available data  | RDS   |

Note that CpG probe annotations are filtered by the same preprocessing steps performed on the methylation data, such that it only provides annotations probes remaining in the preprocessed data. Promoter annotations are generated for protein-coding genes defined in the ENSEMBL gene annotation corresponding to the dataset's genomic build (as provided in the configuration file). Promoters are defined as ±3kB TSS.

### Output Files

**A. Methylation Results:** This compartment stores intermediate stages and results of differential methylation analysis.
| File Name                     | Description                                    | Format    |
|-------------------------------|------------------------------------------------|-----------|
| `<dsn>_raw_dmrs.rds`          | DMRs generated by Bumphunter                   | RDS       |
| `<dsn>_significant_dmrs.rds`  | DMRs filtered for statistical significance     | RDS       |
| `<dsn>_dmgs.rds/xlsx`         | DMGs whose promoters overlap significant DMRs  | RDS, XLSX |
| `<dsn>_slurm_design.rds`      | Bumphunter design input for SLURM              | RDS       |
| `<dsn>_slurm_info.rds`        | Bumphunter genomic info input for SLURM        | RDS       |
| `<dsn>_slurm_data.rds`        | Bumphunter beta value input for SLURM          | RDS       |
| `<dsn>_slurm_analysis.R`      | Bumphunter analysis script for SLURM           | R         |
| `<dsn>_slurm_job.sh`          | Bumphunter job submission script for SLURM     | sh        |

**B. Expression Results:** This compartment stores intermediate stages and results of differential expression analysis.
| File Name                     | Description                                    | Format    |
|-------------------------------|------------------------------------------------|-----------|
| `<dsn>_degs.rds/xlsx`         | DEGs generated by edgeR and Limma              | RDS       |

**C. Integration Results:** This compartment stores intermediate stages and results of integrative analysis.
| File Name                          | Description                                     | Format    |
|------------------------------------|-------------------------------------------------|-----------|
| `<dsn>_background_genes.xlsx`      | All possible genes identifiable by current data | XLSX      |
| `<dsn>_common_genes.rds/xlsx`      | Genes intersecting both DMGs and DEGs           | RDS, XLSX |
| `<dsn>_correlation.rds/xlsx`       | Pearson's Correlation of all common genes       | RDS, XLSX |
| `<dsn>_correlated_genes.rds/xlsx`  | Significantly correlated genes                  | RDS, XLSX |

Note that background genes are defined as protein-coding genes whose promoters, as defined in `<dsn>_promoter_annotation.rds`, has at least one overlapping CpG probe present in the preprocessed methylation data.

## Configuration

## Methods

## Example Dataset
An example of the full pipeline on the Beat AML dataset is provided. The dataset contains DNA methylation data of 298 samples and RNA-seq data of 

## Contact

Contact Sophia Li at sophiamjia.li@utoronto.ca for questions regarding pipeline usage or data access.
