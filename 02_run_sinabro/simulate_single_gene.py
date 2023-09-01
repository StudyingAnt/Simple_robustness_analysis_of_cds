import sys
import pathlib

import pandas as pd

from Bio import SeqIO

# import sinabro
sys.path.insert(0, "/home/jsong/sinabro")
import sinabro.sinabro as snbr

gene_idx = int(sys.argv[1])-1

# CHANGE DIRECTORIES ACCORDINGLY
exp_dir = pathlib.Path(__file__).parent.parent.resolve()
in_dir = pathlib.PurePath(exp_dir, "in")
out_dir = pathlib.PurePath(exp_dir, "out")

genes_fasta_file = pathlib.PurePath(in_dir, "gencode.v40.pc_transcripts.nopary.cdsplus.fa")

seq_records = list(SeqIO.parse(genes_fasta_file, "fasta"))

transcript_name = seq_records[gene_idx].id.split("|")[4]
seq = seq_records[gene_idx].seq

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

sbs_dict = {}
for sbs in signatures:
    sbs_dict[sbs] = []

n_traj = 1000
for sbs in signatures:
    for _ in range(n_traj):
        traj = snbr.Trajectory(transcript_name, seq)
        traj.autofill(condition="nonsynonymous", method="signature", mutational_signature=sbs, note=f"{sbs}")

        sbs_dict[sbs].append(traj._length-1)

out_data = pd.DataFrame(sbs_dict)
out_file = pathlib.PurePath(out_dir, f"{transcript_name}_rus_profile.csv")
out_data.to_csv(out_file, index=False)