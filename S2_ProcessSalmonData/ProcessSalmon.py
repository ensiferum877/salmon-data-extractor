import os
import pandas as pd
import numpy as np
import re
from tqdm import tqdm

def process_and_filter_salmon_data(config):
    input_folder = config['salmon_input_folder']
    output_folder = config['output_folder']
    output_file = os.path.join(output_folder, config['filtered_data_filename'])
    
    # Load patient IDs
    patients = pd.read_csv(config['patients_file'], sep='\t')
    patient_ids = set(patients[config['patients_id_column']].astype(str).tolist())
    
    print(f"Number of patient IDs: {len(patient_ids)}")
    
    # Load top genes
    meta_analysis = pd.read_csv(config['meta_analysis_file'])
    top_genes = set(meta_analysis.nlargest(config['top_genes_count'], config['feature_importance_column'])['Gene'].tolist())
    
    print(f"Number of top genes: {len(top_genes)}")
    
    value_column = config['columns_to_extract'][0]
    print(f"Using value column: {value_column}")
    
    # Process files
    files = [f for f in os.listdir(input_folder) if f.endswith('.sf')]
    print(f"Number of .sf files found: {len(files)}")
    
    gene_data = {gene: {} for gene in top_genes}
    processed_patients = set()
    
    for file in tqdm(files, desc="Processing files"):
        sample_id_match = re.search(config['sample_id_pattern'], file)
        if sample_id_match and sample_id_match.group('sample_id') in patient_ids:
            sample_id = sample_id_match.group('sample_id')
            file_path = os.path.join(input_folder, file)
            
            df = pd.read_csv(file_path, sep='\t', usecols=['Name', value_column])
            df = df[df['Name'].isin(top_genes)].set_index('Name')
            
            for gene in top_genes:
                if gene in df.index and not pd.isna(df.loc[gene, value_column]):
                    gene_data[gene][sample_id] = str(df.loc[gene, value_column])
            
            processed_patients.add(sample_id)
    
    # Remove genes with no data
    genes_with_data = [gene for gene in top_genes if gene_data[gene]]
    genes_removed = top_genes - set(genes_with_data)
    
    # Write consolidated data to file
    with open(output_file, 'w') as f:
        f.write("Gene\t" + "\t".join(sorted(processed_patients)) + "\n")
        for gene in genes_with_data:
            values = [gene_data[gene].get(patient, 'NA') for patient in sorted(processed_patients)]
            f.write(gene + "\t" + "\t".join(values) + "\n")
    
    # Check for missing patients
    missing_patients = patient_ids - processed_patients
    
    print(f"\nFiltered data saved to {output_file}")
    print(f"Processed {len(processed_patients)} patients out of {len(patient_ids)} total patients")
    print(f"Number of missing patients: {len(missing_patients)}")
    if missing_patients:
        print(f"Missing patients: {', '.join(sorted(missing_patients))}")
    
    print(f"\nNumber of genes included in the output: {len(genes_with_data)}")
    print(f"Number of genes removed due to no data: {len(genes_removed)}")
    if genes_removed:
        print(f"Removed genes: {', '.join(sorted(genes_removed))}")
    
    # Load the written file to check for missing values
    df_check = pd.read_csv(output_file, sep='\t', index_col='Gene')
    missing_count = (df_check == 'NA').sum().sum()
    print(f"\nNumber of 'NA' values in the final data: {missing_count}")
    
    return output_file