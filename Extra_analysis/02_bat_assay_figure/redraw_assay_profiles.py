"""Redraw the existing six bat APOBEC motif means on a common percentage axis."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import sys
import tempfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / '06_Bat_genome/experimental_apobec'))
from prepare_profiles import ENZYMES, digest, load_profiles, require

MOTIFS = [a+b+'C' for a in 'TGCA' for b in 'TGCA']


def redraw(profiles_dir, out_dir):
    profiles_dir, out_dir = profiles_dir.resolve(), out_dir.resolve()
    load_profiles(profiles_dir)  # Authenticate both source values and normalized profiles.
    require(not out_dir.exists(), 'Choose a new output directory')
    require(out_dir != profiles_dir and profiles_dir not in out_dir.parents,
            'Output must be outside the prepared source directory')
    source = profiles_dir / 'source_values.csv'
    source_hash = digest(source)
    with source.open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    indexed = {(row['profile_id'], row['motif']): row for row in rows}
    require(len(rows) == len(indexed) == 96, 'Expected 96 distinct enzyme/motif means')
    values = []
    for sample, enzyme, species, *_ in ENZYMES:
        group = [indexed[enzyme, motif] for motif in MOTIFS]
        require(all(row['sample_code'] == sample and row['enzyme_species'] == species for row in group),
                'Unexpected sample/enzymatic profile mapping')
        vector = [float(row['percentage']) for row in group]
        require(all(math.isfinite(x) and 0 <= x <= 100 for x in vector), 'Invalid percentage')
        values.append(vector)
    ymax = max(5, 5 * math.ceil(max(max(v) for v in values) / 5))
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.assay-plot-', dir=out_dir.parent) as temp:
        stage = Path(temp) / 'figure'
        stage.mkdir()
        if 'MPLCONFIGDIR' not in os.environ:
            os.environ['MPLCONFIGDIR'] = str(Path(temp) / 'mpl-cache')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 3, figsize=(11, 6.4), sharey=True)
        for index, (ax, (_, enzyme, species, *_), vector) in enumerate(zip(axes.flat, ENZYMES, values)):
            species_name = 'Molossus molossus' if species == 'molMol2' else 'Myotis myotis'
            ax.bar(range(16), vector, color='#739bb8', edgecolor='black', linewidth=.5)
            ax.set_title(f'{chr(97+index)}. {species_name}\n{enzyme.split("_", 1)[1]}', fontsize=10)
            ax.set_xticks(range(16))
            ax.set_xticklabels(MOTIFS, rotation=60, ha='right', fontsize=8)
            ax.set_ylim(0, ymax)
            ax.set_yticks(range(0, ymax+1, 5))
            ax.tick_params(axis='y', labelleft=True)
            ax.set_ylabel('Mean C-to-T rate (%)')
            ax.set_xlabel('NNC motif (terminal C is mutated)')
        fig.tight_layout()
        for suffix in ('pdf', 'png'):
            fig.savefig(stage / ('figure6_shared_axis.' + suffix), dpi=200)
        plt.close(fig)
        caption = ('Figure 6. Motif preferences of six purified bat APOBEC enzymes. '
                   'Bars show mean C-to-T mutation percentages across occurrences of each NNC motif '
                   'in the 369-bp artificial DNA substrate. All panels use the same percentage axis; '
                   'values are the supplied motif means, without normalization to a sum of 100%. '
                   'Only C-to-T changes are represented. Independent reaction replicate numbers, '
                   'replicate-level variation and negative-control background rates are not specified '
                   'in the supplied records; replicate-based error bars are therefore not shown.\n')
        (stage / 'figure6_caption_draft.txt').write_text(caption)
        require(digest(source) == source_hash, 'Source changed during plotting')
        report = dict(redraw_passed=True, new_assay_data=False, new_rums_simulation=False,
            input_source_values_sha256=source_hash, input_workbook_sha256=json.loads(
                (profiles_dir / 'preparation_manifest.json').read_text())['workbook_sha256'],
            script_sha256=digest(__file__), values_used=96, y_limits=[0, ymax],
            plotted_statistic='Provided motif-average C-to-T percentage; no profile renormalization',
            averaging_unit='Repeated occurrences of the same NNC within one substrate sequence',
            independent_reaction_replicates=None, independent_replicates_inferred_as_one=False,
            replicate_error_bars_available=False, background_rate_available=False,
            caption_status='Draft for author review; missing assay metadata must be resolved or disclosed',
            outputs=[dict(filename=p.name, sha256=digest(p)) for p in sorted(stage.iterdir())])
        (stage / 'figure6_redraw_summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        require(not out_dir.exists(), 'Output appeared during plotting')
        stage.rename(out_dir)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profiles-dir', type=Path, required=True,
                        help='Existing prepared six-profile directory with authenticated source_values.csv')
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(redraw(args.profiles_dir, args.out_dir), indent=2))


if __name__ == '__main__':
    main()
