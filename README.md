# Integrative Analysis of Molecular Subtypes in Acute Myeloid Leukemia

An R-based pipeline for integrated differential analysis of DNA methylation and RNA-seq data in acute myeloid leukemia (AML)

Full data (both provided in this repository and omitted) are available here: https://utoronto-my.sharepoint.com/:f:/g/personal/sophiamjia_li_mail_utoronto_ca/ElHaP2t2nvpAi7xc0wz1zF8B1Em2vXN65rb-lLvH86eH1w?e=xNUiKx (expires June 12 2025)

## Table of Contents
- [Repository Layout](#repository-layout)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Dependency Installation](#dependency-installation)
- [Usage](#usage)
- [Configuration](#configuration)
- [Inputs and Outputs](#inputs-and-outputs)
- [Pipeline Adaption](#pipeline-adaption)
  - [Adapting to other Subtypes](#adapting-to-other-subtypes)
  - [Adapting to other Datasets](#adapting-to-other-datasets)
- [Methods](#methods)
- [Assays Supported](#assays-supported)


## Repository Layout
The following repository layout includes the provided BeatAML example dataset and all data provided for download.

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

## Configuration

## Inputs and Outputs

## Pipeline Adaption
  ### Adapting to other Subtypes
  ### Adapting to other Datasets

## Methods

## Assays Supported
- Cross-reactive probes for EPIC assay fetched from Pidsley et al., 2016
