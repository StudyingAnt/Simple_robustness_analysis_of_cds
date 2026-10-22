import pandas as pd
import numpy as np

from pathlib import PurePath, Path

def entropy(p):
    return -1*np.dot(p, np.log2(p))

def roughness(p):
    return -1*np.log2(1/len(p)) - entropy(p)

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/" # CHANGE HERE
signatures_file = PurePath(base_path, "COSMIC_v3.3.1_SBS_GRCh37.txt")

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

original_sbs = pd.read_csv(signatures_file, sep="\t", index_col=0)

mut_types = []
for mut in ["[C>A]","[C>G]","[C>T]","[T>A]","[T>C]","[T>G]"]:
    for nt1 in ["A", "C", "G", "T"]:
        for nt2 in ["A", "C", "G", "T"]:
            mut_types.append(f"{nt1}{mut}{nt2}")

sbs_dict = {}
for signature in signatures:
    sbs_dict[signature] = []
    for mut_type in mut_types:
        sbs_dict[signature].append(original_sbs[signature][mut_type]) 

sbs_transversion = []
sbs_transition = []
sbs_trv_tr = []
sbs_trs_tr = []
for signature in signatures:
    val_trv = sum(sbs_dict[signature][0:32]) + sum(sbs_dict[signature][48:64]) + sum(sbs_dict[signature][80:])
    val_trs = sum(sbs_dict[signature][32:48]) + sum(sbs_dict[signature][64:80])
    
    sbs_transversion.append(val_trv)
    sbs_transition.append(val_trs)
    sbs_trv_tr.append(val_trv-0.6)
    sbs_trs_tr.append(val_trs-0.4)


data = pd.DataFrame(
        {
            "Signature": signatures,
            "SBS_transversion": sbs_transversion,
            "SBS_transition": sbs_transition,
            "SBS_trv_tr": sbs_trv_tr,
            "SBS_trs_tr": sbs_trs_tr
        }
    )

data.to_csv(PurePath(base_path, f"SBS_transversion_transition.csv"), index=False)
