# Gene Expression Machine Learning (GEML) pipeline 

## Overview

GEML is a comprehensive pipeline for processing gene expression data stored as salmon quantification files and applying multiple machine learning classification algorithms to RNA-seq data. This pipeline is designed to streamline the analysis of transcriptomic data.

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
```bash
pip install pandas numpy tqdm pyyaml
```

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

The `config.yaml` file is the cornerstone of the GEMPD pipeline, allowing users to customize various aspects of the data processing and analysis. Below is a detailed description of each parameter:

### Input/Output Paths
- `salmon_input_folder`: Path to the directory containing salmon quantification files (.sf)
  - Example: "/home/user/data/salmon_quant/"
- `output_folder`: Directory where all output files will be saved
  - Example: "/home/user/results/GEMPD_output/"
- `metadata_file`: Path to the CSV file containing metadata for the samples
  - Example: "/home/user/data/metadata.csv"
- `patients_file`: Path to the TSV file listing selected patient IDs
  - Example: "/home/user/data/selected_patients.tsv"

### Data Processing Parameters
- `columns_to_extract`: List of columns to extract from salmon files
  - Example: ["TPM"]  # Transcript Per Million
- `batch_size`: Number of files to process in each batch (for memory management)
  - Example: 100
- `sample_id_pattern`: Regular expression to extract sample IDs from filenames
  - Example: '(?P<sample_id>\d{4}-SL-\d{4})'

### Filtering Parameters
- `patients_id_column`: Column name in the patients file containing IDs
  - Example: "HudAlphaID"
- `meta_analysis_file`: Path to file containing gene importance information
  - Example: "/home/user/data/gene_importance.csv"
- `feature_importance_column`: Column name for feature importance in meta-analysis file
  - Example: "Feature_Importance_1"
- `top_genes_count`: Number of top genes to include based on importance
  - Example: 151
- `filtered_data_filename`: Name of the output file for filtered data
  - Example: "filtered_data.tsv"
- `proceed_with_filtering`: Boolean to control whether to apply filtering
  - Example: true

### Machine Learning Parameters
- `reference_class`: The reference class for classification tasks
  - Example: "Control"
- `max_cores`: Maximum number of CPU cores to use (null for all available)
  - Example: 4
- `max_models`: Maximum number of ML models to run (null for all available)
  - Example: null

### Example Configuration
```yaml
salmon_input_folder: "/home/user/data/salmon_quant/"
output_folder: "/home/user/results/GEMPD_output/"
metadata_file: "/home/user/data/metadata.csv"
patients_file: "/home/user/data/selected_patients.tsv"

columns_to_extract: ["TPM"]
batch_size: 100
sample_id_pattern: '(?P<sample_id>\d{4}-SL-\d{4})'

patients_id_column: "HudAlphaID"
meta_analysis_file: "/home/user/data/gene_importance.csv"
feature_importance_column: "Feature_Importance_1"
top_genes_count: 151
filtered_data_filename: "filtered_data.tsv"
proceed_with_filtering: true

reference_class: "Control"
max_cores: 4
max_models: null
```

## Output

The GEMPD pipeline produces two main outputs:

1. Filtered Gene Expression Matrix
2. Machine Learning Evaluation Metrics

### 1. Filtered Gene Expression Matrix

Filename: `filtered_data.tsv`

This file contains the processed and filtered gene expression data in a tabular format.

Format:
- Tab-separated values (TSV)
- Rows represent genes
- Columns represent samples

Content:
- First column: Gene identifiers (e.g., ENSEMBL IDs or gene symbols)
- Subsequent columns: Expression values for each sample (e.g., TPM values)


### 2. Machine Learning Evaluation Metrics

Filename: `ml_model_results.csv`.  This file contains:

This file contains the performance metrics and details for each machine learning model applied to the gene expression data.

Format:
- Comma-separated values (CSV)
- Each row represents a different machine learning model or configuration

Columns:
1. `Model`: Name of the machine learning model
2. `Method`: Specific algorithm or method used
3. `Transformation`: Data transformation technique applied (if any)
4. `Normalization`: Data normalization method used
5. `PreProcessing`: Any preprocessing steps applied to the data
6. `Reference`: Reference class used for classification
7. `Accuracy`: Overall accuracy of the model
8. `Kappa`: Cohen's Kappa statistic (measure of agreement)
9. `Time`: Execution time for the model (in minutes)
10. `Error`: Any error messages encountered during model execution (if applicable)

## Troubleshooting

Ensure all file paths in config.yaml are correct and accessible.
Check that all required R and Python libraries are installed.
