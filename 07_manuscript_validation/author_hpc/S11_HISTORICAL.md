# S11: matched-transcript analysis of recovered historical estimates

> Author execution record: these commands and paths describe the authors’ SeaWulf setup. SeaWulf requires institutional access. For ordinary terminal/Python commands, use the stage README one level above. Researchers using an HPC system must adapt these examples to their own site.

This stage regenerates the criterion comparison using the same 3,740 transcripts on both sides. It uses recovered historical simulation estimates and does not launch new Sinabro paths. The bundled gzip CSV contains the transcript IDs/names and the two historical estimates, with source fingerprints and the original ordering recorded in its manifest. Original HDF5 files and raw COSMIC data are not required on the target computer. The pinned human FASTA is still required to verify transcript identity.

## Completed full HPC verification

Job **2345768** on **dg026** passed the full matched-transcript analysis: all **3,740** sequences verified, both populations **3,740**, and the source manifest, value archive and executable script hashes match the checked local run. Mean estimates are **0.8116874332** (identity) and **3.8720010695** (BLOSUM). Reported elapsed time was **6.316 seconds**. The report is at `PROJECT/results/s11_historical_2345768/s11_historical_summary.json`. This confirms historical postprocessing, not a fresh BLOSUM calculation under the corrected package. Existing successful output does not need resubmission.

## SeaWulf: one job

From the project root:

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
git -C Simple_robustness_analysis_of_cds pull --ff-only
mkdir -p logs
sbatch Simple_robustness_analysis_of_cds/07_manuscript_validation/author_hpc/analyze_s11_historical_seawulf.slurm
```

This job reuses the existing `sinabro-revision` environment and `data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa`. It requests one shared CPU, 4 GB and a 30-minute time limit; the time limit is not a runtime estimate. No new environment, HDF5 library, data transfer or mutation simulation is required. The two completed core validation jobs do not need resubmission.

Output directory: `PROJECT/results/s11_historical_JOBID/`.

- `s11_historical_summary.json`: checksums, completion status, scope and numerical summaries.
- `s11_matched_transcript_values.csv`: both archived estimates and their difference for each verified transcript, with CDS lengths and sequence hashes.
- `figureS11_historical_matched.pdf/.png`: two criterion distributions on the same transcript population.
- `s11_historical_criterion_scatter.pdf/.png`: matched-transcript scatter plot.

The Slurm log names are `logs/s11-historical-JOBID.out` and `.err`. After completion, inspect or share `s11_historical_summary.json`. Expected values include:

```json
{
  "analysis_passed": true,
  "historical_results_only": true,
  "fresh_simulation": false,
  "corrected_package_blosum_recalculation": false,
  "n_transcripts": 3740,
  "input_sequences_verified": 3740
}
```

The historical matched-population means are approximately **0.8116874332** under protein identity and **3.8720010695** under BLOSUM. These reference numbers describe the recovered estimates; they are not substituted values from a new simulation.

## Other computers

With the analysis environment active, from the repository root:

```bash
python 07_manuscript_validation/analyze_s11_historical.py \
  --input /PATH/TO/gencode.v40.pc_transcripts.nopary.cdsplus.fa \
  --out-dir /PATH/TO/PROJECT/results/s11_historical
```

Choose a new output directory. The program checks the bundled values, input FASTA and transcript ID/name pairing, and publishes outputs atomically only after all checks succeed. It refuses outputs nested inside an input directory. `SINABRO_ENV`, `SINABRO_FASTA`, `ANALYSIS_DIR`, `PROJECT` and `SINABRO_S11_OUTPUT` can override the Slurm defaults.

## Interpretation and relationship to the original figure

The recovered original SVG contains 73,214 identity points and 3,740 BLOSUM points. This script deliberately uses the identity values for those same 3,740 transcripts, removing that population mismatch. The historical script specifies 1,000 paths per estimate, but exact historical random seeds and per-path counts for this result were not recovered. Matching transcript identities is not matching random paths; no Monte Carlo confidence interval or paired-path test is claimed.

This establishes reproducible postprocessing of the historical numerical source. A fresh BLOSUM sensitivity run under the corrected package would be a separate calculation. The twelve-case illustrative reconstruction is also separate and is not substituted for this recovered original population. BLOSUM distance is an amino-acid alignment proxy, not measured biological functional loss. See [historical input recovery](records/SUPPLEMENTARY_SOURCE_AUDIT.md) for provenance.
