#!/usr/bin/env python3
"""Freeze a deterministic family/molecule/length/annotation-stratified pilot."""
import collections
import csv
import gzip
import hashlib
import json
import random
import time
from pathlib import Path

from Bio import SeqIO

from acquire import ROOT, sha256, utc


def fields(record):
    source = next((f for f in record.features if f.type == 'source'), None)
    qualifiers = source.qualifiers if source else {}
    lineage = record.annotations.get('taxonomy', [])
    family = next((x for x in reversed(lineage) if x.endswith('viridae')), 'unassigned')
    cds = sum(f.type == 'CDS' for f in record.features)
    n = len(record)
    length_bin = '<10kb' if n < 10000 else ('10-100kb' if n < 100000 else '>=100kb')
    molecule = record.annotations.get('molecule_type', 'unknown')
    return {
        'accession': record.id, 'length': n, 'molecule': molecule,
        'family': family, 'cds_count': cds,
        'host_annotation': '|'.join(qualifiers.get('host', []) + qualifiers.get('lab_host', [])),
        'description': record.description,
        'stratum': '|'.join([family, molecule, length_bin, 'CDS' if cds else 'no_CDS'])}


def main():
    start = time.monotonic()
    cache = ROOT / 'data/cache/refseq'
    cache.mkdir(parents=True, exist_ok=True)
    (ROOT / 'results').mkdir(exist_ok=True)
    source = ROOT / 'data/raw/refseq/viral.1.genomic.gbff.gz'
    records, groups, exclusions = [], collections.defaultdict(list), []
    with gzip.open(source, 'rt') as handle:
        for record in SeqIO.parse(handle, 'genbank'):
            row = fields(record)
            records.append(row)
            if not record.seq.defined:
                exclusions.append({'accession':record.id, 'reason':'undefined_sequence_in_GenBank_record',
                                   'length':len(record)})
                continue
            groups[row['stratum']].append(record.id)
    assert len({r['accession'] for r in records}) == len(records)
    rng = random.Random(20261112)
    for name in sorted(groups):
        groups[name].sort()
        rng.shuffle(groups[name])
    selected = []
    eligible = len(records) - len(exclusions)
    while len(selected) < min(10000, eligible):
        for name in sorted(groups):
            if groups[name]:
                selected.append(groups[name].pop())
            if len(selected) == min(10000, eligible):
                break
    selected_set = set(selected)
    (ROOT / 'manifests/pilot_accessions.txt').write_text('\n'.join(sorted(selected)) + '\n')
    with (cache / 'genomes.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(records)
    counts = collections.Counter()
    with gzip.open(source, 'rt') as handle, (cache / 'pilot.fna').open('w') as nt, \
            (cache / 'pilot.faa').open('w') as aa, (cache / 'pilot_proteins.tsv').open('w') as mapping:
        writer = csv.writer(mapping, delimiter='\t')
        writer.writerow(['id', 'protein_accession', 'genome', 'location', 'product'])
        for record in SeqIO.parse(handle, 'genbank'):
            if record.id not in selected_set:
                continue
            SeqIO.write(record, nt, 'fasta')
            counts['pilot_genomes'] += 1
            counts['pilot_bases'] += len(record)
            for index, feature in enumerate(record.features):
                if feature.type != 'CDS':
                    continue
                sequence = feature.qualifiers.get('translation', [None])[0]
                if not sequence:
                    counts['cds_without_translation'] += 1
                    continue
                protein = feature.qualifiers.get('protein_id', ['unassigned'])[0]
                identifier = f'{protein}|{record.id}|{index}'
                aa.write(f'>{identifier}\n{sequence}\n')
                writer.writerow([identifier, protein, record.id, str(feature.location),
                                 feature.qualifiers.get('product', [''])[0]])
                counts['pilot_proteins'] += 1
                counts['pilot_residues'] += len(sequence)
    with (cache / 'all_proteins.faa').open('w') as output:
        for filename in ['viral.1.protein.faa.gz', 'viral.wp_protein.1.protein.faa.gz']:
            with gzip.open(ROOT / 'data/raw/refseq' / filename, 'rt') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    SeqIO.write(record, output, 'fasta')
                    counts['downloaded_proteins'] += 1
                    counts['downloaded_residues'] += len(record)
    (ROOT / 'manifests/corpus_exclusions.json').write_text(json.dumps(exclusions,indent=2)+'\n')
    summary = dict(counts, available_genomes=len(records), eligible_genomes=eligible,
                   excluded_undefined_sequences=len(exclusions),
                   available_bases=sum(r['length'] for r in records),
                   strata=len(groups), selected_strata=len({r['stratum'] for r in records if r['accession'] in selected_set}),
                   selected_with_host_annotation=sum(bool(r['host_annotation']) for r in records if r['accession'] in selected_set),
                   selection='Round-robin across family/molecule/length/annotation strata; shuffle seed 20261112',
                   host_stratification='Not implemented: host qualifiers retained; host taxonomy not inferred',
                   refseq_release=(ROOT / 'data/raw/refseq/RELEASE_NUMBER').read_text().strip(),
                   source_sha256=sha256(source), completed_utc=utc(), seconds=time.monotonic()-start,
                   outputs={p.name:sha256(p) for p in cache.iterdir() if p.is_file()})
    (ROOT / 'results/corpus_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k:v for k,v in summary.items() if k != 'outputs'}, indent=2))


if __name__ == '__main__':
    main()
