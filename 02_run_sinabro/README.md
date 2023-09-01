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