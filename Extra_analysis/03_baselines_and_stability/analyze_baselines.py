"""Exact one-substitution synonymous-neighbour baselines versus saved RUMS."""
import argparse
import csv
import itertools
import json
from pathlib import Path
import re
import tempfile
import time

from common import (np, stats, prior, require, digest, write_csv, correlation,
                    load_human, check_output, finish, add_arguments, REPO)
from Bio import SeqIO
from Bio.Data import CodonTable

BASES = 'ACGT'
CODONS = [''.join(x) for x in itertools.product(BASES, repeat=3)]
TABLE = CodonTable.unambiguous_dna_by_id[1]
AA = [TABLE.forward_table.get(c, '*') for c in CODONS]
SYN = np.zeros((64, 3, 4), dtype=bool)
for i, codon in enumerate(CODONS):
    for pos, alt in itertools.product(range(3), range(4)):
        mutant = codon[:pos] + BASES[alt] + codon[pos+1:]
        SYN[i, pos, alt] = BASES[alt] != codon[pos] and TABLE.forward_table.get(mutant, '*') == AA[i]
ENCODE = np.full(256, 255, dtype=np.uint8)
for i, base in enumerate(BASES):
    ENCODE[ord(base)] = i


def signature_table(path, names):
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        key = reader.fieldnames[0]
        rows = list(reader)
    labels = [r[key] for r in rows]
    require(len(labels) == len(set(labels)) == 96, 'Expected 96 distinct channels')
    require(all(re.fullmatch(r'[ACGT]\[[CT]>[ACGT]\][ACGT]', s) and s[2] != s[4] for s in labels),
            'Invalid SBS96 channel')
    weights = np.array([[float(r[name]) for name in names] for r in rows])
    require(np.isfinite(weights).all() and (weights >= 0).all() and (weights.sum(axis=0) > 0).all(),
            'Invalid signature table')
    return labels, weights / weights.sum(axis=0)


def channel_lookup(labels):
    index = {s: i for i, s in enumerate(labels)}
    result = np.full((64, 4), -1, dtype=int)
    for i, context in enumerate(CODONS):
        left, ref, right = [BASES.index(x) for x in context]
        for alt in range(4):
            if alt == ref:
                continue
            a, b, c, d = (left, ref, alt, right) if ref in (1, 3) else (3-right, 3-ref, 3-alt, 3-left)
            result[i, alt] = index[f'{BASES[a]}[{BASES[b]}>{BASES[c]}]{BASES[d]}']
    require((np.bincount(result[result >= 0], minlength=96) == 2).all(), 'Reverse-complement channel mapping failed')
    return result


def opportunities(sequence, lookup):
    encoded = ENCODE[np.frombuffer(sequence.upper().encode('ascii'), dtype=np.uint8)]
    require(len(encoded) >= 5 and (len(encoded)-2) % 3 == 0 and (encoded < 4).all(), 'Invalid CDS-plus DNA')
    cds = encoded[1:-1]
    triplets = cds.reshape(-1, 3)
    codon = np.repeat(triplets[:, 0]*16 + triplets[:, 1]*4 + triplets[:, 2], 3)
    frame = np.tile(np.arange(3), len(triplets))
    contexts = encoded[:-2]*16 + encoded[1:-1]*4 + encoded[2:]
    all_counts, neutral_counts = np.zeros(96, dtype=np.int64), np.zeros(96, dtype=np.int64)
    for alt in range(4):
        keep = cds != alt
        channels = lookup[contexts[keep], alt]
        neutral = SYN[codon[keep], frame[keep], alt]
        all_counts += np.bincount(channels, minlength=96)
        neutral_counts += np.bincount(channels[neutral], minlength=96)
    require(all_counts.sum() == 3*len(cds) and (neutral_counts <= all_counts).all(), 'Invalid neighbour counts')
    return all_counts, neutral_counts


def analyze(args):
    started = time.perf_counter()
    data, selected = load_human(args)
    require(digest(args.signatures) == data['manifest']['configuration']['signature_sha256'],
            'Signature source differs from the saved simulation')
    labels, weights = signature_table(args.signatures, data['signatures'])
    lookup = channel_lookup(labels)
    output = check_output(args.out_dir, [data['folder'], data['fasta'].parent, args.signatures.parent])
    chosen = {data['ids'][i]: n for n, i in enumerate(selected)}
    total = np.empty((len(selected), 96), dtype=np.int64)
    neutral = np.empty_like(total)
    seen = set()
    for record in SeqIO.parse(data['fasta'], 'fasta'):
        identifier = record.id.split('|')[0]
        if identifier not in chosen:
            continue
        require(identifier not in seen, 'Duplicate FASTA ID')
        seen.add(identifier)
        i = chosen[identifier]
        total[i], neutral[i] = opportunities(str(record.seq), lookup)
        if len(seen) % 1000 == 0:
            print(f'Exact baseline counts {len(seen)}/{len(selected)}', flush=True)
    require(seen == set(chosen), 'Missing selected CDSs')
    denominator = total @ weights
    require((denominator > 0).all(), 'No available initial event for a saved complete case')
    weighted = (neutral @ weights) / denominator
    uniform = neutral.sum(axis=1) / total.sum(axis=1)
    observed = data['matrix'][selected]
    require(np.isfinite(weighted).all() and (weighted >= 0).all() and (weighted <= 1+1e-12).all(), 'Invalid baseline probability')
    summary, discordant = [], []
    for j, name in enumerate(data['signatures']):
        r_rank = stats.rankdata(-observed[:, j], method='average')
        b_rank = stats.rankdata(-weighted[:, j], method='average')
        differences = np.abs(r_rank-b_rank)
        top = max(1, int(np.ceil(.05*len(selected))))
        # Rank cutoffs include all tied observations at the top-5% threshold.
        r_top = observed[:, j] >= np.sort(observed[:, j])[-top]
        b_top = weighted[:, j] >= np.sort(weighted[:, j])[-top]
        summary.append(dict(Signature=name, n_genes=len(selected), mean_rums=float(observed[:, j].mean()),
            mean_weighted_one_step_probability=float(weighted[:, j].mean()),
            spearman_uniform_vs_rums=correlation(uniform, observed[:, j]),
            spearman_weighted_vs_rums=correlation(weighted[:, j], observed[:, j]),
            median_absolute_rank_shift=float(np.median(differences)),
            top5_rums_count=int(r_top.sum()), top5_baseline_count=int(b_top.sum()),
            top5_intersection=int((r_top & b_top).sum()), top5_union=int((r_top | b_top).sum())))
        order = sorted(range(len(selected)), key=lambda n: (-differences[n], data['ids'][selected[n]]))[:10]
        for n in order:
            i = selected[n]
            discordant.append(dict(Signature=name, Transcript=data['ids'][i],
                gene_id=data['annotations'][i]['gene_id'], gene_name=data['annotations'][i]['gene_name'],
                cds_nt=data['features'][i]['cds_nt'], gc3=data['features'][i]['gc3'],
                uniform_one_step_probability=float(uniform[n]), weighted_one_step_probability=float(weighted[n,j]),
                rums=float(observed[n,j]), rums_rank=float(r_rank[n]), baseline_rank=float(b_rank[n]),
                absolute_rank_shift=float(differences[n])))
    with tempfile.TemporaryDirectory(prefix='.one-step-', dir=output.parent) as temp:
        stage = Path(temp) / 'results'
        stage.mkdir()
        write_csv(stage/'baseline_correlations_by_signature.csv', summary)
        write_csv(stage/'rank_discordant_examples.csv', discordant)
        write_csv(stage/'uniform_one_step_baseline.csv', [dict(Transcript=data['ids'][i],
            gene_id=data['annotations'][i]['gene_id'], cds_nt=data['features'][i]['cds_nt'],
            neutral_neighbours=int(neutral[n].sum()), total_neighbours=int(total[n].sum()),
            uniform_synonymous_fraction=float(uniform[n])) for n,i in enumerate(selected)])
        with (stage/'weighted_one_step_baseline.csv').open('w', newline='') as stream:
            writer = csv.writer(stream);writer.writerow(['Transcript', *data['signatures']])
            writer.writerows([[data['ids'][i], *weighted[n]] for n,i in enumerate(selected)])
        fig, axes = prior.paper.plt.subplots(1, 2, figsize=(9,4))
        axes[0].scatter([r['mean_weighted_one_step_probability'] for r in summary],
                        [r['mean_rums'] for r in summary], s=18)
        axes[0].set(xlabel='Mean initial weighted synonymous probability', ylabel='Mean saved RUMS')
        axes[1].plot(range(len(summary)), [r['spearman_weighted_vs_rums'] for r in summary], '.', label='Signature weighted')
        axes[1].plot(range(len(summary)), [r['spearman_uniform_vs_rums'] for r in summary], '.', label='Uniform neighbours')
        axes[1].set(xlabel='Signature index in declared table order', ylabel='Spearman rank correlation', ylim=(-1,1))
        axes[1].legend()
        prior.paper.save_plot(fig, stage, 'one_step_baseline_comparison')
        prior.paper.read_export(data['folder'])
        require(digest(args.signatures) == data['manifest']['configuration']['signature_sha256']
                and digest(data['fasta']) == data['manifest']['input_sha256'], 'Source changed during analysis')
        report = dict(analysis_passed=True, stage='one_step_baselines', fresh_simulation=False,
            bounded_subset=args.allow_subset, representative_genes=len(selected), signatures=len(summary),
            input_export_manifest_sha256=data['manifest_sha256'], fasta_sha256=digest(data['fasta']),
            signature_sha256=digest(args.signatures), script_sha256=digest(__file__),
            common_sha256=digest(Path(__file__).parent/'common.py'),
            uniform_definition='Synonymous single-base neighbours / (3 * CDS nt), equal event weights',
            weighted_definition='Initial synonymous event weight / initial total event weight, both strands equal',
            criterion='Standard genetic code protein identity; terminal/internal stop symbols retained',
            sequence_contexts_updated_in_baseline=False, baseline_is_rums_or_fixed_n_robustness=False,
            genetic_code_and_codon_usage_adjusted_causal_inference=False,
            interpretation='Exact first-step baselines, not a functional assay, adaptation null, n-robustness or multivariable causal test.',
            across_signature_mean_spearman=correlation(
                np.array([r['mean_weighted_one_step_probability'] for r in summary]),
                np.array([r['mean_rums'] for r in summary])),
            focal_profiles=[r for r in summary if r['Signature'] in ('SBS2','SBS13','SBS19')],
            elapsed_seconds=time.perf_counter()-started)
        finish(stage, output, report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    parser.add_argument('--signatures', type=Path, required=True)
    return analyze(parser.parse_args(argv))


if __name__ == '__main__':
    main()
