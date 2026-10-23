"""Read-only access to the previously audited human inputs."""
import csv
import json
import os
from pathlib import Path
import sys

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / '01_statistics_and_sequence_features'))
import run_analysis as prior
import numpy as np
from scipy import stats

require, digest, write_csv = prior.require, prior.digest, prior.write_csv


def load_human(args):
    data = prior.load_dataset(args.human_export, args.human_fasta)
    cfg = data['manifest']['configuration']
    require(cfg['n_trajectories'] == args.n_traj, 'Unexpected saved trajectory count')
    require(cfg.get('runner') == 'rums-v3' and cfg.get('translation_table') == 1
            and cfg.get('cds_start') == 1 and cfg.get('cds_end') == 'len(sequence)-1'
            and cfg.get('sinabro') == '3.0.0', 'Unsupported simulation/translation/coordinate policy')
    if not args.allow_subset:
        require(len(data['ids']) == 73214 and len(data['signatures']) == 79
                and args.n_traj == 1000, 'Expected the full 73,214 x 79 x 1,000 human export')
    selected = prior.representatives(data)
    require(args.allow_subset or len(selected) == 19309, 'Unexpected representative gene count')
    return data, selected


def check_output(output, inputs):
    output = output.resolve()
    require(not output.exists(), 'Output already exists; choose a new directory')
    require(all(output != p.resolve() and p.resolve() not in output.parents
                and output not in p.resolve().parents for p in inputs), 'Output overlaps inputs')
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def correlation(x, y):
    return None if np.ptp(x) == 0 or np.ptp(y) == 0 else float(stats.spearmanr(x, y)[0])


def finish(stage, output, report):
    report['outputs'] = [dict(filename=p.name, sha256=digest(p), bytes=p.stat().st_size)
                         for p in sorted(stage.iterdir()) if p.is_file()]
    (stage / 'analysis_summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    require(not output.exists(), 'Output appeared during analysis')
    stage.rename(output)
    print(json.dumps(report, indent=2), flush=True)
    print('Report:', output / 'analysis_summary.json', flush=True)


def add_arguments(parser):
    parser.add_argument('--human-export', type=Path, required=True)
    parser.add_argument('--human-fasta', type=Path)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--n-traj', type=int, default=1000)
    parser.add_argument('--allow-subset', action='store_true', help='Explicit software validation only')
