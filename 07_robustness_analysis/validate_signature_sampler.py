"""Recreate original S1/S2 random-sequence validation with a conditional oracle."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import sys
import tempfile
import time

for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / '02_run_sinabro'))
import run_rums as package_helpers
from analyze_robustness_results import digest, require, write_csv, save_plot
import matplotlib.pyplot as plt


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT','TGCA'))[::-1]


def conditional_channel_probabilities(sequence, labels, profile):
    # Independent event enumeration by counts of the 64 trinucleotide contexts.
    counts = Counter(sequence[i:i+3] for i in range(len(sequence)-2))
    opportunity = np.array([counts[label[0]+label[2]+label[-1]] +
                            counts[reverse_complement(label[0]+label[2]+label[-1])] for label in labels])
    weights = profile * opportunity
    require(weights.sum() > 0, 'No eligible signature event in sequence')
    return weights / weights.sum()


def channel_checks(observed, expected, variance, n_sequences):
    mask = variance > 0
    require(np.all(observed[~mask] == expected[~mask]), 'Impossible channel observed')
    residual = np.zeros(len(observed))
    residual[mask] = (observed[mask]-expected[mask])/np.sqrt(variance[mask])
    tested = mask & (expected >= 10) & (n_sequences-expected >= 10)
    maximum = float(np.abs(residual[tested]).max()) if tested.any() else 0.
    # Poisson-binomial Chernoff bound also checks unexpectedly frequent rare bins.
    log_bound = -expected.copy()
    positive = (observed > 0) & (expected > 0)
    log_bound[positive] = observed[positive] - expected[positive] - observed[positive]*np.log(observed[positive]/expected[positive])
    threshold = np.log(1e-8 / (2 * 96 * 79))
    return dict(passed=bool(maximum < 8. and float(log_bound.min()) >= threshold),
                largest_standardized_residual=maximum, tested_channels=int(tested.sum()),
                smallest_log_chernoff_bound=float(log_bound.min()), log_chernoff_threshold=float(threshold)), residual


def sample_signature(task):
    settings, signature = task
    package_helpers.load_package(Path(settings['sinabro_source']))
    from sinabro import core
    path = Path(settings['signature_path'])
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        channel = reader.fieldnames[0]
        rows = list(reader)
    labels = [row[channel] for row in rows]
    profile = np.array([float(row[signature]) for row in rows])
    profile /= profile.sum()
    require(len(labels) == len(set(labels)) == 96, 'Invalid source SBS96 table')
    key = hashlib.sha256(f'{settings["seed"]}\0S1_S2\0{signature}'.encode()).digest()
    rng = np.random.default_rng(int.from_bytes(key[:8], 'big'))
    base_prob = [(1-settings['gc'])/2, settings['gc']/2, settings['gc']/2, (1-settings['gc'])/2]
    observed = np.zeros(96, dtype=np.int64)
    expected = np.zeros(96)
    variance = np.zeros(96)
    started = time.perf_counter()
    for _ in range(settings['n_sequences']):
        sequence = ''.join(rng.choice(list('ACGT'), size=settings['length'], p=base_prob))
        oracle = conditional_channel_probabilities(sequence, labels, profile)
        expected += oracle
        variance += oracle * (1-oracle)
        event = core.signature_mutation(sequence, mutational_signature='custom', custom_signature_path=path, column=signature,
                                        start=1, end=len(sequence)-1, rng=rng)
        require(event.e == 0 and 1 <= event.idx_target < len(sequence)-1, 'Invalid sampled event/boundary')
        label = event.mut_type
        if core.parse_motif(label)[1] not in 'CT':
            label = core.reverse_motif(label)
        index = labels.index(label)
        require(oracle[index] > 0, 'Sampled a channel without eligible opportunity')
        observed[index] += 1
    checks, residual = channel_checks(observed,expected,variance,settings['n_sequences'])
    cosine = lambda a,b: float(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)))
    return dict(Signature=signature, n_sequences=settings['n_sequences'], oracle_check_passed=checks['passed'],
                channel_checks=checks,
                cosine_to_original_profile=cosine(observed,profile), cosine_to_conditional_expectation=cosine(observed,expected),
                elapsed_seconds=time.perf_counter()-started,
                channels=[dict(Signature=signature, channel=label, observed_count=int(observed[i]),
                    expected_count=float(expected[i]), source_fraction=float(profile[i]),
                    standardized_residual=float(residual[i])) for i,label in enumerate(labels)])


def validate(args):
    require(args.workers > 0 and args.n_sequences > 0 and args.length >= 20 and 0 < args.gc < 1
            and args.seed >= 0, 'Invalid sampler validation parameters')
    package_helpers.load_package(args.sinabro_source)
    source = package_helpers.source
    path = source / 'sinabro/data/COSMIC_v3.3.1_SBS_GRCh37.txt'
    with path.open(newline='') as stream:
        names = next(csv.reader(stream, delimiter='\t'))[1:]
    selected = args.signatures or names
    require(len(selected) == len(set(selected)) and set(selected) <= set(names), 'Unknown/duplicate signatures')
    output = args.out_dir.resolve()
    require(not output.exists(), 'Choose a new sampler output directory')
    settings = dict(sinabro_source=str(source), signature_path=str(path), seed=args.seed,
                    n_sequences=args.n_sequences, length=args.length, gc=args.gc)
    report = dict(validation_passed=False, original_paper_design=dict(sequence_length=args.length,
        gc_fraction=args.gc, sequences_per_signature=args.n_sequences, mutations_per_sequence=1),
        seed=args.seed, signature_sha256=digest(path), package_source_sha256=package_helpers.package_digest(),
        script_sha256=digest(__file__), workers=args.workers, signatures=len(selected),
        total_single_substitutions=args.n_sequences*len(selected), results=[])
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.paper-sampler-',dir=output.parent) as temporary:
        stage = Path(temporary) / 'results'
        stage.mkdir()
        tasks = [(settings,name) for name in selected]
        rows, channels = [], []
        with mp.get_context('spawn').Pool(args.workers) as pool:
            for result in pool.imap(sample_signature,tasks):
                channels.extend(result.pop('channels'))
                rows.append(result)
                print(f'{result["Signature"]}: cosine={result["cosine_to_original_profile"]:.6f}; oracle={result["oracle_check_passed"]}',flush=True)
        write_csv(stage/'signature_sampler_summary.csv',rows)
        write_csv(stage/'signature_channel_counts.csv',channels)
        fig,ax=plt.subplots(figsize=(max(8,len(selected)*.12),3))
        ax.plot(range(len(rows)),[r['cosine_to_original_profile'] for r in rows],'.',label='COSMIC profile')
        ax.plot(range(len(rows)),[r['cosine_to_conditional_expectation'] for r in rows],'.',label='Sequence-conditional oracle')
        ax.set_xticks(range(len(rows)),[r['Signature'] for r in rows],rotation=90,fontsize=5)
        ax.set(ylim=(0,1.02),ylabel='Cosine similarity')
        ax.legend()
        save_plot(fig,stage,'figureS2_sampler_cosine')
        # All 96 bins remain in channel order, separately for each signature.
        for name in selected:
            group=[r for r in channels if r['Signature']==name]
            fig,ax=plt.subplots(figsize=(12,2.8))
            ax.bar(np.arange(96),[r['source_fraction'] for r in group],alpha=.5,label='COSMIC')
            ax.bar(np.arange(96),[r['observed_count']/args.n_sequences for r in group],alpha=.5,label='Simulated')
            ax.set(title=name,xlabel='SBS96 channel',ylabel='Fraction')
            ax.legend()
            save_plot(fig,stage,'figureS1_'+name)
        require(digest(path)==report['signature_sha256'] and package_helpers.package_digest()==report['package_source_sha256'],
                'Package/signature changed during sampler validation')
        report.update(validation_passed=all(row['oracle_check_passed'] for row in rows), results=rows,
            mean_cosine_to_original=float(np.mean([r['cosine_to_original_profile'] for r in rows])),
            interpretation='Observed channel fractions depend on sequence opportunities. Original-profile cosine is descriptive; sampler correctness is compared with an independently enumerated conditional oracle. Original 0.974 is not enforced.',
            outputs=[dict(filename=p.name,sha256=digest(p)) for p in sorted(stage.iterdir())])
        (stage/'sampler_validation_summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        stage.rename(output)
    require(report['validation_passed'],'Sampler conditional-oracle validation failed; inspect report')
    print(json.dumps({k:v for k,v in report.items() if k not in ('results','outputs')},indent=2),flush=True)
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sinabro-source',type=Path,required=True)
    parser.add_argument('--out-dir',type=Path,required=True)
    parser.add_argument('--n-sequences',type=int,default=10000)
    parser.add_argument('--length',type=int,default=2200)
    parser.add_argument('--gc',type=float,default=.4)
    parser.add_argument('--seed',type=int,default=42)
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--signatures',nargs='+')
    validate(parser.parse_args(argv))


if __name__=='__main__':
    main()
