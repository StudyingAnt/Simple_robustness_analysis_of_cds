#!/bin/bash
#$ -N rusAnal
#$ -S /bin/bash
#$ -o /home/jsong/stdouts
#$ -e /home/jsong/stderrs
#$ -j y
#$ -q cpu_short
#$ -P tomprj
#$ -t 1-73214

conda activate py39

python3 /titan/jsong/robustness_analysis_all_genes/scripts/simulate_single_gene.py ${SGE_TASK_ID}

conda deactivate
