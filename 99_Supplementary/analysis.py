import pathlib
import re

import pandas as pd

from pathlib import PurePath, Path



# import required libraries
import numpy as np
from numpy.linalg import norm

def cosine_similarity(a, b):
    a = np.array(a)
    b = np.array(b)
    return np.dot(a,b)/(norm(a)*norm(b))


proj_data_dir = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/signature_interaction/original_data"
exp_dir = pathlib.Path(__file__).parent.parent.resolve()
input_dir = PurePath(exp_dir, "in")
output_dir = PurePath(exp_dir, "out")

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

mut_types = []
for mut in ["[C>A]","[C>G]","[C>T]","[T>A]","[T>C]","[T>G]"]:
    for nt1 in ["A", "C", "G", "T"]:
        for nt2 in ["A", "C", "G", "T"]:
            mut_types.append(f"{nt1}{mut}{nt2}")



original_sbs_file = PurePath(output_dir, f"COSMIC_v3.3.1_SBS_GRCh37.txt")

original_sbs = pd.read_csv(original_sbs_file, sep="\t", index_col=0)

for signature in signatures:
    mut_type_cnts = {}
    for mut_type in mut_types:
        mut_type_cnts[mut_type] = original_sbs[signature][mut_type]
    
    data = pd.DataFrame(
                    {
                        "mut_type": mut_type_cnts.keys(),
                        "percentage": mut_type_cnts.values()
                    }
                )
    output_path = Path(f"/mnt/c/Users/CEEL-PC-005/Desktop/Joon/signature_interaction/exp/20221207-validation/out/{signature.lower()}")
    output_path.mkdir(parents=True, exist_ok=True)

    data.to_csv(f"/mnt/c/Users/CEEL-PC-005/Desktop/Joon/signature_interaction/exp/20221207-validation/out/{signature.lower()}/{signature}.csv", index=False)



