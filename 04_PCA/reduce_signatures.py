import pandas as pd

from pathlib import Path, PurePath

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/" # CHANGE HERE
all_rums_file = PurePath(base_path, "all_gencode_rums_profile.csv")

all_rums = pd.read_csv(all_rums_file)

signatures_noseqerr = [
    "Transcript",
    "SBS1", "SBS2", "SBS3", "SBS4", "SBS5", 
    "SBS6", "SBS7a", "SBS7b", "SBS7c", "SBS7d",   
    "SBS8", "SBS9", "SBS10a", "SBS10b", "SBS10c", 
    "SBS10d", "SBS11", "SBS12", "SBS13", "SBS14",   
    "SBS15", "SBS16", "SBS17a", "SBS17b", "SBS18", 
    "SBS19", "SBS20", "SBS21", "SBS22", "SBS23",   
    "SBS24", "SBS25", "SBS26",  "SBS28", "SBS29", 
    "SBS30", "SBS31", "SBS32", "SBS33", "SBS34",   
    "SBS35", "SBS36", "SBS37", "SBS38", "SBS39", 
    "SBS40", "SBS41", "SBS42", "SBS44", "SBS84",   
    "SBS85", "SBS86", "SBS87", "SBS88", "SBS89", 
    "SBS90", "SBS91", "SBS92", "SBS93", "SBS94"
]
  
signatures_aetiology = [
    "Transcript",
    "SBS1", "SBS2", "SBS13", "SBS84", "SBS85",     # Deamination activity
    "SBS4", "SBS29", "SBS92",                      # Tobacco exposure
    "SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS38",   # UV exposure
    "SBS18",                                       # ROS exposure
    "SBS22", "SBS24", "SBS42", "SBS88", "SBS90",   # Chemical exposure
    "SBS11", "SBS25", "SBS31", "SBS32", "SBS35",   # Anti-cancer treatment
    "SBS86", "SBS87",
    "SBS3", "SBS6", "SBS9", "SBS10a", "SBS10b",    # Defective DNA repair and
    "SBS10c", "SBS10d", "SBS14", "SBS15", "SBS20", # hypermutations 
    "SBS21", "SBS26", "SBS30", "SBS36", "SBS44"
]

noseqerr_rums = all_rums.loc[:, signatures_noseqerr]
aetiology_rums = all_rums.loc[:, signatures_aetiology]

file_name = PurePath(base_path, "all_gencode_rums_profile_noseqerr.csv")
noseqerr_rums.to_csv(file_name, index=False)
file_name = PurePath(base_path, "all_gencode_rums_profile_aetiology.csv")
aetiology_rums.to_csv(file_name, index=False)
