# import packages
from pathlib import Path, PurePath
import pandas as pd


# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/" # CHANGE HERE

# Step 1: Load Ensembl gene IDs from TXT file
filename_txt = "Cosmic_CancerGeneCensus_v99_GRCh38.ensemblid.txt"  # Adjust the path if needed
id_file = Path(PurePath(base_path, filename_txt))
with open(id_file, 'r') as file:
    gene_ids = {line.strip() for line in file}

# Step 2: Parse the GTF file and filter records
filename_gtf = "gencode.v40.basic.annotation.gtf"  # Adjust the path if needed
gtf_file = Path(PurePath(base_path,filename_gtf ))
filtered_records = []

with open(gtf_file, 'r') as file:
    for line in file:
        if line.startswith('#'):
            continue  # Skip header lines
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue  # Skip incomplete lines
        attributes = fields[8]
        gene_id_field = next((attr for attr in attributes.split(';') if 'gene_id' in attr), None)
        if gene_id_field:
            gene_id = gene_id_field.split('"')[1].split('.')[0]  # Extract the gene ID, remove version if present
            if gene_id in gene_ids and '_PAR_Y' not in gene_id_field:
                filtered_records.append(fields)  # Store the fields as a list for easier DataFrame creation

# Convert filtered records to a DataFrame
df = pd.DataFrame(filtered_records, columns=['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

# Save to a TSV file
output_filename = 'gencode.v40.basic.annotation.CancerGeneCensus.gtf'
df.to_csv(output_filename, sep='\t', index=False)

print(f'Data saved to {output_filename}')
