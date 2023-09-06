# Principal component anlysis (PCA): GC contents of wobble postion

## Reduce columns
There are total 79 mutational signatures released in COSMIC. However, 19 of them are potential sequencing artifacts and another 19 of them do not have known aetiology. We removed those signatures and generated two processed files: `all_gencode_rums_profile_noseqerr.csv` and `all_gencode_rums_profile_aetiology.csv`.
```console
python3 reduce_signatures.py
```

## Run PCA
We ran PCA on data with all signatures, without sequencing artifacts, and with only known aetiology. Further, analysis was performed mostly with PCA without sequencing artifacts.
```console
python3 run_pca.py
```

## Draw Figure
One can use R script `draw_violin_plots_by_rank.R` to draw **Figure XX**.