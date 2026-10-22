import pandas as pd
import numpy as np

from Bio import SeqIO
from pathlib import Path, PurePath

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/" # CHANGE HERE
sequences_file = PurePath(base_path, f"gencode.v40.pc_transcripts.nopary.cdsplus.fa")

transcript_names = []
gc_totals = []
gc_wobbles = []
gc_ratios = []
for seq_record in SeqIO.parse(sequences_file, "fasta"):
    transcript_name = seq_record.id.split("|")[4]
    seq = seq_record.seq

    third_pos = {"A": 0, "C": 0, "G": 0, "T": 0}
    all_pos = {"A": 0, "C": 0, "G": 0, "T": 0}
    
    for i in range(1, len(seq)-1, 1):
        if i%3 == 0:
            if seq[i] == "A":
                all_pos["A"] = all_pos["A"]+1
                third_pos["A"] = third_pos["A"]+1
            elif seq[i] == "C":
                all_pos["C"] = all_pos["C"]+1
                third_pos["C"] = third_pos["C"]+1
            elif seq[i] == "G":
                all_pos["G"] = all_pos["G"]+1
                third_pos["G"] = third_pos["G"]+1
            elif seq[i] == "T":
                all_pos["T"] = all_pos["T"]+1
                third_pos["T"] = third_pos["T"]+1
        else:
            if seq[i] == "A":
                all_pos["A"] = all_pos["A"]+1
            elif seq[i] == "C":
                all_pos["C"] = all_pos["C"]+1
            elif seq[i] == "G":
                all_pos["G"] = all_pos["G"]+1
            elif seq[i] == "T":
                all_pos["T"] = all_pos["T"]+1

    total_length = len(seq)-2
    wobble_num = total_length/3
    all_freq = {key: value/total_length for key, value in all_pos.items()}
    third_freq = {key: value/wobble_num for key, value in third_pos.items()}
    gc_total = (all_pos["C"]+all_pos["G"])/total_length
    gc_wobble = (third_pos["C"]+third_pos["G"])/wobble_num

    transcript_names.append(transcript_name)
    gc_totals.append(gc_total)
    gc_wobbles.append(gc_wobble)
    gc_ratios.append(gc_wobble/gc_total)

output = pd.DataFrame.from_dict(
    {
        "Transcript_name": transcript_names,
        "GC_total": gc_totals,
        "GC_wobble": gc_wobbles,
        "GC_ratio": gc_ratios
    }
)

output_file = PurePath(base_path, f"gencode.v40.pc_transcripts.nopary.cdsplus.gc_contents.csv")
output.to_csv(output_file, index=False)