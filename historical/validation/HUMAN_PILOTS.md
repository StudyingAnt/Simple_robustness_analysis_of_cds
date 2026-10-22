# Small timing pilot with the refactored package

For the current full-run workflow, use [portable execution and audit](../../02_run_sinabro/REPRODUCTION.md), [SeaWulf submission](../../02_run_sinabro/HPC_SUBMISSION.md), and [audited NPZ-to-table export](../../02_run_sinabro/AGGREGATION.md). The user reported that the full human saved-output audit passed on 2026-10-16: 73,214 transcripts, 79 signatures and 1,000 trajectories per pair. The pilot descriptions and historical instructions below are retained as earlier stages.

The checked pilots use Sinabro 3.0.0, the recorded scientific-library environment, and the freshly reproduced GENCODE v40 CDS-plus file. Both 20 and 100 trajectories per combination have been measured. The 1,000-trajectory manuscript setting has not been run in this pilot.

## Preserved 20-path reference

```bash
python 02_run_sinabro/benchmark_pilot.py \
  --input data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa \
  --sinabro-source ../sinabro \
  --out-dir results/benchmark20_2026-10-14 \
  --n-traj 20 --seed 42 --max-steps 10000 --case-timeout 180

python 02_run_sinabro/evaluate_benchmark_pilot.py \
  --results results/benchmark20_2026-10-14
```

The benchmark selects three distinct genes per CDS length bin (up to 1,000; 1,001–3,000; 3,001–10,000; above 10,000 nt), including the longest transcript. Internal-stop records are excluded only from this timing sample. The unchanged complete preprocessing data remains available.

There are 12 transcripts × SBS1/SBS2/SBS13 × 20 trajectories = **720 trajectories**, with one measured batch per combination and no extra trajectory warm-up. Each combination runs in an isolated serial worker with numerical-library threads limited to one. Profile loading is prepared before timing and included in setup and end-to-end costs. The current memory instrumentation uses Linux/WSL peak worker RSS.

Outputs include selection CSV/FASTA, input and source hashes, seeds, environment, per-combination wall/CPU time, peak worker RSS, per-path stopping reasons, mutation counts, and trajectory hashes. Every path is checked for a single substitution per step, correct CDS coordinates, preserved flanks, protein identity before the terminal event, and correct accepted-step counts. The output directory must be new to protect previous measurements.

The evaluator audits a completed result and writes CSV/JSON summaries, a report and time/memory plots. It **does not simulate any trajectories**. The original 20-path report remains preserved in its output directory.

## 100 paths per combination

```bash
python 02_run_sinabro/benchmark_pilot.py \
  --input data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa \
  --sinabro-source ../sinabro \
  --out-dir results/benchmark100_2026-10-14 \
  --n-traj 100 --seed 42 --max-steps 10000 --case-timeout 180

python 02_run_sinabro/evaluate_benchmark_pilot.py \
  --results results/benchmark100_2026-10-14 \
  --baseline results/benchmark20_2026-10-14
```

The same 12 transcripts and 3 signatures produce 3,600 trajectories. The case seed does not depend on the requested count, allowing the first 20 paths in each 100-path case to be compared with the earlier measurement. All 720 prefix paths matched in this run.

The historical simulate_single_gene.py sets n_traj = 1000. This parameter is the number of independent trajectories per transcript/signature pair; each trajectory stops at the first protein-sequence change. The pilot reduces that repeat count to 100. A 1,000-path estimate in the 100-path report applies only to the same 12 transcripts and 3 signatures, and is not a measured result or a whole-genome projection.

## Preserved historical full-analysis instructions

The instructions below describe the original cluster workflow. Their import paths and fixed paths must be reviewed before use with Sinabro 3.0.0. For current audited NPZ results, use aggregate_results.py rather than the historical CSV combiner below.

# Run Sinabro

## Sinabro
Here we utilized Python package [Sinabro](https://github.com/StudyingAnt/sinabro) we developed to measure the robustness under mutational signature (RUMS).

One can import **Sinabro** by adding into `PATH` using `sys` module as in the following code:
```python
# import sinabro
sys.path.insert(0, "/home/jsong/sinabro")
import sinabro.sinabro as snbr
```

Or, use `conda develop /PATH/TO/SINABRO/` if you are using `conda` environment.
```console
conda develop /PATH/TO/SINABRO/
```

## Run Sinbaro on SGE cluster
We run total of 1,000 simulations per sequence for more than 70,000 sequences using SGE based computing cluster at Laufer center. Following is the job script to run all simulations using `simulate_single_gene.py`. Parameters should be modified according to the cluster you want to use.
```console
#!/bin/bash
#$ -N NAME
#$ -S /bin/bash
#$ -o /PATH/TO/STDOUTS
#$ -e /PATH/TO/STDERRS
#$ -j y
#$ -q cpu_short
#$ -P PROJECT_NAME
#$ -t 1-73214

conda activate py39

python3 /titan/jsong/robustness_analysis_all_genes/scripts/simulate_single_gene.py ${SGE_TASK_ID}

conda deactivate
```

## Combine results into single file
We used `combine_rums_profile.py` to combine results into single `all_gencode_rums_profile.csv` file.
```console
python3 combine_rums_profile.py
```
