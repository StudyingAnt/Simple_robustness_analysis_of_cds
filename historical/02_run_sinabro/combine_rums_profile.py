import pathlib

import pandas as pd

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/" # CHANGE HERE
input_dir = pathlib.PurePath(base_path, "rums_analysis")

signatures = [
        "SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS6", "SBS7a", "SBS7b", 
        "SBS7c", "SBS7d", "SBS8", "SBS9", "SBS10a", "SBS10b", "SBS10c", 
        "SBS10d", "SBS11", "SBS12", "SBS13", "SBS14", "SBS15", "SBS16", 
        "SBS17a", "SBS17b", "SBS18", "SBS19", "SBS20", "SBS21", "SBS22", 
        "SBS23", "SBS24", "SBS25", "SBS26", "SBS27", "SBS28", "SBS29", 
        "SBS30", "SBS31", "SBS32", "SBS33", "SBS34", "SBS35", "SBS36", 
        "SBS37", "SBS38", "SBS39", "SBS40", "SBS41", "SBS42", "SBS43", 
        "SBS44", "SBS45", "SBS46", "SBS47", "SBS48", "SBS49", "SBS50", 
        "SBS51", "SBS52", "SBS53", "SBS54", "SBS55", "SBS56", "SBS57",
        "SBS58", "SBS59", "SBS60", "SBS84", "SBS85", "SBS86", "SBS87", 
        "SBS88", "SBS89", "SBS90", "SBS91", "SBS92", "SBS93", "SBS94", 
        "SBS95"
        ]

with open(pathlib.PurePath(base_path, "all_gencode_rums_profile.csv"), "w") as f:
    line = ",".join(signatures)
    f.write(f"Transcript,{line}\n")
    for child in pathlib.Path(input_dir).iterdir():
        transcript_name = child.name.split("_rums_")[0]
        print(transcript_name)
        cnts_transcripts = pd.read_csv(child)
        line_list = [str(item) for item in cnts_transcripts.mean().to_list()]
        line = ",".join(line_list)
        f.write(f"{transcript_name},{line}\n")


    