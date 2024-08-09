# Salmon quantification data processing with the evaluation of multiple ML classification algorithms

## Overview

GEMPD is a comprehensive pipeline for processing gene expression data stored as salmon quantification files and applying multiple machine learning classification algorithms to RNA-seq data. This pipeline is designed to streamline the analysis of transcriptomic data.

## Table of Contents

1. [Features](#features)
2. [Prerequisites](#prerequisites)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Pipeline Steps](#pipeline-steps)
6. [Configuration](#configuration)
7. [Output](#output)
8. [Troubleshooting](#troubleshooting)

## Features

- Process salmon quantification files into a consolidated format
- Filter and align gene expression data
- Apply multiple machine learning models for classification
- Parallel processing for improved performance
- Configurable parameters for flexibility
- Comprehensive output including model performance metrics

## Prerequisites

- Python 3.7+
- R 4.0+
- Jupyter Notebook
- Required Python libraries: pandas, numpy, tqdm, yaml
- Required R libraries: MLSeq, DESeq2, S4Vectors, kernlab, parallel, doParallel, foreach

## Installation

1. Clone the repository:
[git clone https://github.com/your-username/GEMPD.git](https://github.com/ensiferum877/salmon-rna-ml-pipeline)

2. Install required Python libraries:

pip install pandas numpy tqdm pyyaml

3. Install required R libraries:
```R
install.packages(c("MLSeq", "DESeq2", "S4Vectors", "kernlab", "parallel", "doParallel", "foreach"))
```

## Usage

1. Update the config.yaml file with your specific parameters and file paths.
2. Run the Jupyter notebook **PipelineNotebook.ipynb** to execute the entire pipeline.

## Pipeline Steps

Load Configuration: The pipeline starts by loading parameters from config.yaml.
Process Salmon Data:

Consolidates salmon quantification files into a tabular format
Filters data based on specified genes and patients
Aligns data with metadata


Machine Learning Evaluation:

Prepares data for machine learning
Creates DESeq objects
Runs multiple machine learning models in parallel
Evaluates model performance


## Configuration

Key parameters in config.yaml:

salmon_input_folder: Path to the folder containing salmon quantification files
output_folder: Path for output files
metadata_file: Path to the metadata CSV file
patients_file: Path to the file containing selected patient IDs
columns_to_extract: Columns to extract from salmon files (e.g., "TPM")
top_genes_count: Number of top genes to include in the analysis
reference_class: Reference class for classification
max_cores: Maximum number of cores to use for parallel processing
max_models: Maximum number of machine learning models to run

## Output

Consolidated and filtered gene expression data: filtered_data.tsv
Machine learning model results: ml_model_results.csv

The ml_model_results.csv file contains:

Model name
Method
Transformation/Normalization technique
Accuracy
Kappa statistic
Execution time

## Troubleshooting

Ensure all file paths in config.yaml are correct and accessible.
Check that all required R and Python libraries are installed.
