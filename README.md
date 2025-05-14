# Integrative Analysis of Molecular Subtypes in Acute Myeloid Leukemia

_Repository: aml-epi-transcriptome-integration_  

An R-based, automated, and reproducible workflow for integrated DNA methylation and RNA-seq analysis in acute myeloid leukemia. Implements differential methylation (DMG) and expression (DEG) comparisons between NPM1-mutant and PML-RARA fusion subtypes from the Beat AML cohort, with modular design such that adapter files can be added to extend the pipeline to other AML subtypes or datasets.

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
The following repository layout includes the provided BeatAML example dataset. 

```bash
.
├── config/
│   ├── default_config.R
│   └── beataml_config.R
│
├── utilities/
│   └── setup.R
│
├── adapters/
│   ├── methylation_adapter.R
│   ├── expression_adapter.R
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
│   │   ├── raw/
│   │   └── processed/
│   └── beataml/
│       ├── raw/
│       │   ├── beataml_beta_values.csv
│       │   └── beataml_star_gene_counts.csv
│       └── processed/
│           ├── beataml_methylation_data.rds
│           └── beataml_expression_data.rds
│
├── results/
│   ├── default/
│   │   ├── 01_methylation/
│   │   │   └── slurm/
│   │   │   │   ├── scripts/
│   │   │   │   └── data/
│   │   ├── 02_expression/
│   │   └── 03_integration
│   └── beataml/
│       ├── 01_methylation/
│       │   ├── slurm/
│       │   │   ├── scripts/
│       │   │   │   ├── beataml_slurm.sh
│       │   │   │   └── beataml_slurm.R
│       │   │   ├── data/
│       │   │   │   ├── beataml_slurm_data.sh
│       │   │   │   ├── beataml_slurm_info.sh
│       │   │   │   └── beataml_slurm_design.R
│       │   ├── beat_dmgs.rds
│       │   ├── beat_dmgs.xlsx
│       │   ├── beat_methylation_summary.rds
│       │   └── beat_methylation_summary.xlsx
│       ├── 02_expression/
│       │   ├── beat_degs.rds
│       │   └── beat_degs.xlsx
│       └── 03_integration
│           ├── beat_anticorrelated.rds
│           └── beat_anticorrelated.xlsx
│
├── run_pipeline.R
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
