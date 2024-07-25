# Salmon Data Extractor

## Description
Salmon Data Extractor is a Python-based tool designed to consolidate and process data from multiple Salmon quantification files. It's particularly useful for bioinformaticians and researchers working with RNA-seq data, offering functionality similar to the `tximport` package in R, but implemented in Python.

## Features
- Consolidates multiple Salmon quantification files into a single comprehensive table
- Flexible data extraction (e.g., TPM, NumReads, or other columns from Salmon output)
- Parallel processing for efficient handling of large datasets
- Custom filtering options for both rows (genes) and columns (samples)
- Designed to handle large-scale RNA-seq datasets, such as those from the Parkinson's Progression Markers Initiative (PPMI)

## Requirements
- Python 3.7+
- pandas
- tqdm
- pyyaml

## Installation
1. Clone this repository:
git clone https://github.com/ensiferum877/salmon-data-extractor.git
Copy2. Navigate to the project directory:
cd salmon-data-extractor
Copy3. Install required packages:
pip install -r requirements.txt


## Usage
1. Edit the `config.yaml` file to specify your input and output settings:
```yaml

1 Configuration Options

input_folder: 'path/to/your/salmon/quant/files'
output_file: 'output_filename.tsv'
columns_to_extract: 
  - 'TPM'
  - 'NumReads'
batch_size: Number of files to process in each batch (default: 100)
sample_id_pattern: Regular expression pattern to extract sample IDs from filenames

2. Run the main script:

python salmon_processing.py

The consolidated data will be saved in the specified output file.

Filtering
The tool also includes options for filtering the consolidated dataset:

Row filtering: Specify genes of interest
Column filtering: Select specific samples

Edit the filtering section in the script to customize your data selection.
Examples

# Example of running the script with custom filters
filter_dataset(input_file='consolidated_data.tsv', 
               output_file='filtered_data.tsv', 
               row_condition=gene_list, 
               column_condition=sample_list)

Contributing
Contributions to improve the Salmon Data Extractor are welcome. Please feel free to submit a Pull Request.
License
This project is licensed under the MIT License - see the LICENSE file for details.

Acknowledgments

Salmon developers for the quantification tool
Inspired by the functionality of tximport in R
