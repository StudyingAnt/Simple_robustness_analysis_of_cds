import pandas as pd
import numpy as np
from scipy import stats

from pathlib import Path, PurePath
from sklearn.decomposition import PCA

# directory to input and output files
base_path = "/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Simple_robustness_analysis_of_cds_files/" # CHANGE HERE

all_rums_file = PurePath(base_path, "all_gencode_rums_profile.csv")
noseqerr_rums_file = PurePath(base_path, "all_gencode_rums_profile_noseqerr.csv")
aetiology_rums_file = PurePath(base_path, "all_gencode_rums_profile_aetiology.csv")

all_rums = pd.read_csv(all_rums_file)
noseqerr_rums = pd.read_csv(noseqerr_rums_file)
aetiology_rums = pd.read_csv(aetiology_rums_file)

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

signatures_noseqerr = [
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

# run PCA
# for all signatures
data = all_rums.loc[:, all_rums.columns != "Transcript"].to_numpy()
pca = PCA()

pcm = pca.fit(data).components_

pcs = []
sbss = []
values = []
for i in range(0,len(signatures)):
    for j, sbs in enumerate(signatures):
        pcs.append(f"PC{i+1}")
        sbss.append(sbs)
        values.append(pcm[i,j])

pcm_out = pd.DataFrame.from_dict(
    {
        "Principal_component": pcs,
        "Signature": sbss,
        "Coefficient": values
    }
)

pcm_out_file = PurePath(base_path, f"all_gencode_signatures_PCA_matrix.csv")
pcm_out.to_csv(pcm_out_file, index=False)

x = pca.fit_transform(data)

output = pd.concat([all_rums["Transcript"], pd.DataFrame(x, columns=[f"PC{i+1}" for i in range(0,len(signatures))])], axis=1)

output_file = PurePath(base_path, f"all_gencode_signatures_PCA_transformed.csv")
output.to_csv(output_file, index=False)


# for sequencing errors removed signatures
data = noseqerr_rums.loc[:, noseqerr_rums.columns != "Transcript"].to_numpy()
pca = PCA()

pcm = pca.fit(data).components_

pcs = []
sbss = []
values = []
for i in range(0,len(signatures_noseqerr)):
    for j, sbs in enumerate(signatures_noseqerr):
        pcs.append(f"PC{i+1}")
        sbss.append(sbs)
        values.append(pcm[i,j])

pcm_out = pd.DataFrame.from_dict(
    {
        "Principal_component": pcs,
        "Signature": sbss,
        "Coefficient": values
    }
)

pcm_out_file = PurePath(base_path, f"all_gencode_noseqerr_signatures_PCA_matrix.csv")
pcm_out.to_csv(pcm_out_file, index=False)

x = pca.fit_transform(data)

output = pd.concat([noseqerr_rums["Transcript"], pd.DataFrame(x, columns=[f"PC{i+1}" for i in range(0,len(signatures_noseqerr))])], axis=1)

output_file = PurePath(base_path, f"all_gencode_noseqerr_signatures_PCA_transformed.csv")
output.to_csv(output_file, index=False)


# for known aetiology signatures
data = aetiology_rums.loc[:, aetiology_rums.columns != "Transcript"].to_numpy()
pca = PCA()

pcm = pca.fit(data).components_

pcs = []
sbss = []
values = []
for i in range(0,len(signatures_aetiology)):
    for j, sbs in enumerate(signatures_aetiology):
        pcs.append(f"PC{i+1}")
        sbss.append(sbs)
        values.append(pcm[i,j])

pcm_out = pd.DataFrame.from_dict(
    {
        "Principal_component": pcs,
        "Signature": sbss,
        "Coefficient": values
    }
)

pcm_out_file = PurePath(base_path, f"all_gencode_aetiology_signatures_PCA_matrix.csv")
pcm_out.to_csv(pcm_out_file, index=False)

x = pca.fit_transform(data)

output = pd.concat([aetiology_rums["Transcript"], pd.DataFrame(x, columns=[f"PC{i+1}" for i in range(0,len(signatures_aetiology))])], axis=1)

output_file = PurePath(base_path, f"all_gencode_aetiology_signatures_PCA_transformed.csv")
output.to_csv(output_file, index=False)

"""
# Save PC_i PC_j separatedly

# run PCA
# for all signatures
out_path = PurePath(base_path, f"PCA_all_gencode_signatures")
Path(out_path).mkdir(parents=True, exist_ok=True)
data = all_rums.loc[:, all_rums.columns != "Transcript"].to_numpy()
pca = PCA()
x = pca.fit_transform(data)

for i in range(0,len(signatures)):
    for j in range(0,len(signatures)):
        if i != j:
            output = pd.DataFrame.from_dict({
                "Transcript": all_rums["Transcript"].tolist(),
                "X": x[:, i],
                "Y": x[:, j]
                })
            
            output_file = PurePath(out_path, f"all_gencode_signatures_PC{i+1}_PC{j+1}.csv")
            output.to_csv(output_file, index=False)

# for sequencing errors removed signatures
out_path = PurePath(base_path, f"PCA_all_gencode_noseqerr_signatures")
Path(out_path).mkdir(parents=True, exist_ok=True)
data = noseqerr_rums.loc[:, noseqerr_rums.columns != "Transcript"].to_numpy()
pca = PCA()
x = pca.fit_transform(data)

for i in range(0,len(signatures_noseqerr)):
    for j in range(0,len(signatures_noseqerr)):
        if i != j:
            output = pd.DataFrame.from_dict({
                "Transcript": noseqerr_rums["Transcript"].tolist(),
                "X": x[:, i],
                "Y": x[:, j]
                })
            
            output_file = PurePath(base_path, f"all_gencode_noseqerr_signatures_PC{i+1}_PC{j+1}.csv")
            output.to_csv(output_file, index=False)

# for known aetiology signatures
out_path = PurePath(base_path, f"PCA_all_gencode_aetiology_signatures")
Path(out_path).mkdir(parents=True, exist_ok=True)
data = aetiology_rums.loc[:, aetiology_rums.columns != "Transcript"].to_numpy()
pca = PCA()
x = pca.fit_transform(data)

for i in range(0,len(signatures_aetiology)):
    for j in range(0,len(signatures_aetiology)):
        if i != j:
            output = pd.DataFrame.from_dict({
                "Transcript": aetiology_rums["Transcript"].tolist(),
                "X": x[:, i],
                "Y": x[:, j]
                })
            
            output_file = PurePath(base_path, f"all_gencode_aetiology_signatures_PC{i+1}_PC{j+1}.csv")
            output.to_csv(output_file, index=False)
"""