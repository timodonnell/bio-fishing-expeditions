#!/usr/bin/env python3
"""Summarize the deposited proteins and reviewed matches around the viral-context lead."""
import json
from pathlib import Path

from Bio import SeqIO

from acquire import ROOT, sha256, utc
from search import BLAST_FIELDS
from summarize_metavr import read_tsv, write_tsv


def main():
    output = ROOT / 'results/metavr/soil_neighbors'
    locations = json.loads((output / 'locations.json').read_text())
    labels = {r.id: r.description for r in SeqIO.parse(ROOT / 'data/cache/references/swissprot.faa', 'fasta')}
    hits = read_tsv(output / 'swissprot.tsv', BLAST_FIELDS)
    rows = []
    for row in locations:
        matches = [h for h in hits if h['qseqid'] == row['protein']]
        best = max(matches, key=lambda h: float(h['bitscore'])) if matches else {}
        rows.append({'protein': row['protein'], 'length': row['length'],
                     'locations': json.dumps(row['locations']), 'best_match': best.get('sseqid', ''),
                     'evalue': best.get('evalue', ''), 'bitscore': best.get('bitscore', ''),
                     'description': labels.get(best.get('sseqid', ''), '')})
    write_tsv(output / 'gene_review.tsv', rows)
    genome = 'IMGVR_UViG_3300079393_000044'
    candidate = next(r for r in rows if r['protein'].endswith('_14604_15686'))
    neighbors = [next(r for r in rows if r['protein'].endswith(suffix))
                 for suffix in ['_14195_14554', '_15715_16191', '_17133_17618']]
    assert all(float(r['evalue']) < 1e-8 for r in neighbors)
    inputs = [output / 'locations.json', output / 'swissprot.tsv',
              ROOT / 'data/cache/references/swissprot.faa']
    report = {'completed_utc': utc(), 'genome': genome, 'proteins': len(rows),
              'candidate': candidate, 'nearby_viral_matches': neighbors,
              'input_sha256': {str(p): sha256(p) for p in inputs},
              'script_sha256': sha256(Path(__file__)),
              'gene_review_sha256': sha256(output / 'gene_review.tsv'),
              'assessment': 'Multiple adjacent genes have reviewed iridovirus homologs, supporting viral assembly context. Ig-variable-like similarity does not establish antibody/TCR ancestry. Read linkage and phylogenetic discrimination remain pending.'}
    (output / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print('Summarized', len(rows), 'proteins; viral context supported, receptor ancestry unresolved')


if __name__ == '__main__':
    main()
