"""Deterministic GENCODE length-stratified sample shared by validation scripts."""
BINS = [('short', 0, 1000), ('medium', 1000, 3000),
        ('long', 3000, 10000), ('very_long', 10000, 10**12)]
SIGNATURES = ['SBS1', 'SBS2', 'SBS13']


def select_transcripts(input_file):
    from Bio import SeqIO
    candidates = {label: [] for label, _, _ in BINS}
    populations = {label: 0 for label, _, _ in BINS}
    exclusions = {label: 0 for label, _, _ in BINS}
    for record in SeqIO.parse(input_file, 'fasta'):
        text = str(record.seq)
        cds = text[1:-1]
        label = next(name for name, low, high in BINS if low < len(cds) <= high)
        populations[label] += 1
        # Timing ordinary complete ORFs first; retain the full preprocessing dataset.
        internal_stop = any(cds[i:i + 3] in ('TAA', 'TAG', 'TGA')
                            for i in range(0, len(cds) - 3, 3))
        if internal_stop:
            exclusions[label] += 1
            continue
        fields = record.id.split('|')
        candidates[label].append(dict(transcript_id=fields[0], gene_id=fields[1],
                                      transcript_name=fields[4], gene=fields[5],
                                      fasta_id=record.id, cds_nt=len(cds), bin=label))
    selected = []
    used_genes = set()
    for label, _, _ in BINS:
        group = sorted(candidates[label], key=lambda row: (row['cds_nt'], row['transcript_id']))
        if len(group) < 3:
            raise ValueError(f'Not enough eligible transcripts in {label}')
        quantiles = [.1, .5, 1.] if label == 'very_long' else [.1, .5, .9]
        for quantile in quantiles:
            target = group[round(quantile * (len(group) - 1))]['cds_nt']
            available = [row for row in group if row['gene_id'] not in used_genes]
            choice = min(available, key=lambda row: (abs(row['cds_nt'] - target), row['transcript_id']))
            selected.append(dict(choice, target_length_quantile=quantile))
            used_genes.add(choice['gene_id'])
    keep = {row['fasta_id'] for row in selected}
    records = [record for record in SeqIO.parse(input_file, 'fasta') if record.id in keep]
    if len(records) != 12:
        raise ValueError('Selection must contain exactly 12 transcripts')
    return selected, records, populations, exclusions
