# import packages
from pathlib import Path, PurePath

import pandas as pd

# Bat species
bat = "myoMyo6"

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/" # CHANGE HERE
bat_genome_path = Path(PurePath(base_path, "Bat_genome", "results", bat+"_apobec"))

signatures = [
    "APOBEC3A", "APOBEC3C", "APOBEC1"
]

with open(PurePath(base_path, bat+"_rums_apobec.csv"), "w") as f:
    line = ",".join(signatures)
    f.write(f"Transcript,{line}\n")
    for child in Path(bat_genome_path).iterdir():
        transcript_name = child.name.split("_rus_")[0]
        print(transcript_name)
        cnts_transcripts = pd.read_csv(child)
        line_list = [str(item) for item in cnts_transcripts.mean().to_list()]
        line = ",".join(line_list)
        f.write(f"{transcript_name},{line}\n")


    