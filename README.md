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
9. [Contributing](#contributing)
10. [License](#license)

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
```python
pip install pandas numpy tqdm pyyaml

4. 3. Install required R libraries:
```R
install.packages(c("MLSeq", "DESeq2", "S4Vectors", "kernlab", "parallel", "doParallel", "foreach"))
