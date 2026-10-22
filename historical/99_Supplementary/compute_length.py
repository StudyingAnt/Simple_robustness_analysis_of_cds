# import packages
from pathlib import Path, PurePath

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/" # CHANGE HERE

# GENCODE release number
release_num = 40


# make files' names and files' pathes
input_file_name = f"gencode.v{release_num}.pc_transcripts.nopary.cdsplus.fa"
input_file = Path(PurePath(base_path, input_file_name))


# for every record in the file
seq_name = []
seq_length = []
for cnt, seq_record in enumerate(SeqIO.parse(input_file, "fasta")):   
    seq_id = seq_record.id
    seq_id_list = seq_id.split("|")

    print(f"{seq_id_list[4]}\t{len(seq_record.seq)-2}")
    seq_length.append(len(seq_record.seq)-2)
    seq_name.append(seq_id_list[4])

# Create DataFrame
df = pd.DataFrame({
    'transcript_name': seq_name,
    'length': seq_length
})

output_file_name = f"gencode.v{release_num}.pc_transcripts.nopary.cdsplus.length.csv"
output_file = Path(PurePath(base_path, output_file_name))

# Save to CSV without index
df.to_csv(output_file, index=False)
