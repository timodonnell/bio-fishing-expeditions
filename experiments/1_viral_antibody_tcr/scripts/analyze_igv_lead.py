#!/usr/bin/env python3
"""Recover and verify exact alignments for the MetaVR candidate with viral flanks."""
import csv
import json
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO

from acquire import ROOT, sha256, utc
from search import BLAST_FIELDS
from summarize_metavr import read_tsv, write_tsv


def main():
    output = ROOT / 'results/metavr/lead_followup'
    output.mkdir(exist_ok=True)
    identifier = 'IMGVR_UViG_3300079393_000044|Ga0643767_0002625|Ga0643767_0002625_14604_15686'
    neighbors = ROOT / 'results/metavr/soil_neighbors'
    query = next(r for r in SeqIO.parse(neighbors / 'proteins.faa', 'fasta') if r.id == identifier)
    SeqIO.write([query], output / 'query.faa', 'fasta')
    reference_path = ROOT / 'data/cache/references/swissprot.faa'
    previous = json.loads((neighbors / 'swissprot.run.json').read_text())
    assert sha256(reference_path) in previous['inputs'].values()
    executable = previous['search_command'][0]
    database = previous['search_command'][previous['search_command'].index('-db') + 1]
    fields = BLAST_FIELDS + ['nident', 'positive', 'gaps', 'qseq', 'sseq']
    command = [executable, '-query', str(output / 'query.faa'), '-db', database,
               '-evalue', '0.01', '-num_threads', '4', '-max_target_seqs', '1000000',
               '-seg', 'yes', '-outfmt', '6 ' + ' '.join(fields),
               '-out', str(output / 'swissprot_alignments.tsv')]
    subprocess.run(command, check=True)
    hits = read_tsv(output / 'swissprot_alignments.tsv', fields)
    hits.sort(key=lambda h: float(h['bitscore']), reverse=True)
    ids = {r['sseqid'] for r in hits}
    references = {r.id: r for r in SeqIO.parse(reference_path, 'fasta') if r.id in ids}
    best = next(r for r in hits if r['sseqid'] == 'sp|P01637|KV5A5_MOUSE')
    assert best == hits[0]
    original = next(r for r in read_tsv(neighbors / 'swissprot.tsv', BLAST_FIELDS)
                    if r['qseqid'] == identifier and r['sseqid'] == best['sseqid'])
    assert {key: best[key] for key in BLAST_FIELDS} == original
    reference = references[best['sseqid']]
    assert best['qseq'].replace('-', '') == str(query.seq)[int(best['qstart'])-1:int(best['qend'])]
    assert best['sseq'].replace('-', '') == str(reference.seq)[int(best['sstart'])-1:int(best['send'])]
    assert sum(a == b for a, b in zip(best['qseq'], best['sseq'])) == int(best['nident']) == 33
    metadata_path = ROOT / 'data/raw/candidates/P01637.json'
    metadata = json.loads(metadata_path.read_text())
    assert metadata['sequence']['value'] == str(reference.seq)
    assert not metadata.get('genes')
    regions = [f for f in metadata['features'] if f['type'] == 'Region']
    annotated = []
    qpos, spos = int(best['qstart'])-1, int(best['sstart'])-1
    for column, (q, s) in enumerate(zip(best['qseq'], best['sseq'])):
        qpos += q != '-'
        spos += s != '-'
        region = next((f['description'] for f in regions if s != '-' and
                       f['location']['start']['value'] <= spos <= f['location']['end']['value']), '')
        annotated.append({'column_1based': column+1, 'query_position': qpos if q != '-' else None,
                          'reference_position': spos if s != '-' else None,
                          'query_residue': q, 'reference_residue': s,
                          'identity': q == s, 'reference_region': region})
    assert not any(r['reference_region'] == 'Complementarity-determining-3' for r in annotated)
    write_tsv(output / 'alignment_columns.tsv', annotated)
    top = [{k: h[k] for k in BLAST_FIELDS} | {'description': references[h['sseqid']].description}
           for h in hits[:12]]
    write_tsv(output / 'top_reviewed_matches.tsv', top)
    SeqIO.write([reference, references['sp|P01638|KV5A6_MOUSE']], output / 'antibody_matches.faa', 'fasta')
    lines = ['BLASTP local alignment; coordinates refer to deposited precursor sequences.',
             'Candidate 21-107; P01637 23-108; 33/89 identities; 49/89 positives; 5 gap positions.',
             'E = 3.18e-8 against reviewed UniProt 2026_04; bit score = 55.1.',
             'The reference CDR3 (109-117) and framework 4 are outside this alignment.', '']
    q, s = best['qseq'], best['sseq']
    for offset in range(0, len(q), 45):
        a, b = q[offset:offset+45], s[offset:offset+45]
        lines += ['Candidate  ' + a, '           ' + ''.join('|' if x == y else ' ' for x, y in zip(a,b)),
                  'Mouse T1   ' + b, '']
    (output / 'alignment.txt').write_text('\n'.join(lines) + '\n')
    inputs = [neighbors / 'proteins.faa', neighbors / 'swissprot.tsv', neighbors / 'locations.json',
              neighbors / 'gene_review.tsv', metadata_path, ROOT / 'data/raw/candidates/P01638.json',
              ROOT / 'results/metavr/hits.tsv', reference_path]
    figure_data = {'genome': identifier.split('|')[0], 'genome_length': 17620, 'candidate_id': identifier,
                   'candidate_sequence': str(query.seq), 'best_alignment': best, 'alignment_columns': annotated,
                   'reference': {'accession': 'P01637', 'name': metadata['proteinDescription']['recommendedName']['fullName']['value'],
                                 'sequence': str(reference.seq), 'entry_audit': metadata['entryAudit'],
                                 'gene_or_allele_in_record': None, 'features': metadata['features']},
                   'profile': {'name': 'IGK_V', 'query_start': 25, 'query_end': 109, 'global_domain_E': 0.00039545778330796823},
                   'neighbors': read_tsv(neighbors / 'gene_review.tsv'),
                   'query_coverage_percent': 100*(int(best['qend'])-int(best['qstart'])+1)/len(query),
                   'allele_assignment': 'Not established; short distant local similarity and near-tied alternative matches',
                   'novelty': 'Unresolved; viral Ig-like domains have published precedents; no natural antibody/TCR capture established'}
    (output / 'figure_data.json').write_text(json.dumps(figure_data, indent=2) + '\n')
    products = ['query.faa', 'antibody_matches.faa', 'swissprot_alignments.tsv', 'top_reviewed_matches.tsv',
                'alignment_columns.tsv', 'alignment.txt', 'figure_data.json']
    report = {'status': 'complete', 'completed_utc': utc(), 'command': [sys.executable, *sys.argv],
              'script_sha256': sha256(Path(__file__)), 'search_command': command,
              'binary_sha256': sha256(Path(executable)), 'version': subprocess.check_output([executable,'-version'],text=True),
              'database_release': 'UniProt 2026_04, 576073 reviewed proteins',
              'input_sha256': {str(p): sha256(p) for p in inputs},
              'output_sha256': {name: sha256(output/name) for name in products},
              'checks': ['Original BLAST score/coordinates reproduced', 'Aligned sequence slices match FASTA',
                         '33 identities recomputed', 'Current UniProt sequence equals frozen sequence',
                         'No gene/allele field in P01637', 'Reference CDR3 not aligned']}
    (output / 'analysis_manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    print('Verified exact alignment and reference annotation;', len(hits), 'reviewed comparison rows')


if __name__ == '__main__':
    main()
