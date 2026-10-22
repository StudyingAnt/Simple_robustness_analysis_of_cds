# import packages
from pathlib import Path, PurePath

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/" # CHANGE HERE
bat_genome_path = Path(PurePath(base_path, "Bat_genome"))

# test bed
#molMol_bed = Path(PurePath(bat_genome_path, "test.bed"))
molMol_bed = Path(PurePath(bat_genome_path, "HLmolMol2.hg38.toga.bed"))
molMol_fasta = Path(PurePath(bat_genome_path, "HLmolMol2.fa"))

# output files
# output for Sinabro
output_file_cdsplus_name = f"HLmolMol2.pc_transcripts.cdsplus.fa"
output_file_cdsplus = Path(PurePath(bat_genome_path, output_file_cdsplus_name))

# output for sanity check
nopass_names = []
for i in range(1,8,1):
    nopass_names.append(f"HLmolMol2.pc_transcripts.parity{i}.fa")

output_file_nopasses = []
for i in range(1,8,1):
    output_file_nopass = Path(
        PurePath(bat_genome_path, nopass_names[i-1])
    )
    output_file_nopasses.append(output_file_nopass)

# Read in BED file
bed_df = pd.read_csv(molMol_bed, sep='\t', header=None, 
                     names=['chrom', 'chromStart', 'chromEnd', 'name', 
                            'score', 'strand', 'thickStart', 'thickEnd', 
                            'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'])

# print(bed_df)

# Read in genome fasta save into dict
molMol_genome = {}
for seq_record in SeqIO.parse(molMol_fasta, "fasta"):
    molMol_genome[seq_record.id] = seq_record

# make empty lists to store new record
records_cdsplus = []
records_nopass_dict = {}
for i in range(1,8,1):
    records_nopass_dict[f"parity{i}"] = []

for index, row in bed_df.iterrows():
    chrom = row['chrom']
    chromStart = int(row['chromStart'])
    chromEnd = int(row['chromEnd'])
    name = row['name']
    strand = row['strand']
    blockCount = int(row['blockCount'])
    blockSizes = row['blockSizes'].split(",")[:-1]
    blockSizes = [int(char) for char in blockSizes]
    blockStarts = row['blockStarts'].split(",")[:-1]
    blockStarts = [int(char) for char in blockStarts]

    parity = 0

    seq_str = ""
    for i in range(blockCount):
        blockStart = chromStart+blockStarts[i]
        blockEnd = chromStart+blockStarts[i]+blockSizes[i]

        seq_str = seq_str + molMol_genome[chrom].seq[blockStart:blockEnd].upper()

    five_prime = molMol_genome[chrom].seq[chromStart-1]
    three_prime = molMol_genome[chrom].seq[chromEnd]

    seq = Seq(five_prime+seq_str+three_prime)

    if strand == "-":
        seq = seq.reverse_complement()

    # Filtering
    # 1. Is the length multiple of 3?
    if (len(seq)-2)%3 != 0:
        parity = parity + 1
    
    # 2. Does it start with start codon?
    if seq[1:4] != "ATG":
        parity = parity + 2
        
    # 3. Does it ends with stop codon?
    if seq[-4:-1] not in ["TAA", "TAG", "TGA"]:
        parity = parity + 4

    # collect sequences that passes all check
    record_cdsplus = SeqRecord(seq, id=name ,description="")
    if parity == 0:
        records_cdsplus.append(record_cdsplus)
    else:
        records_nopass_dict[f"parity{parity}"].append(record_cdsplus)

# save passed sequences
SeqIO.write(records_cdsplus, output_file_cdsplus, "fasta")

# save sequences that are not passed by parity
for i in range(1,8,1):
    SeqIO.write(records_nopass_dict[f"parity{i}"], output_file_nopasses[i-1], "fasta")
        