#!/usr/bin/env python3
"""Validate completed candidate comparisons, assembly mappings, and review artifacts."""
import collections
import json
from pathlib import Path

from Bio import SeqIO

from acquire import ROOT, sha256, utc
from metavr_context import locate
from summarize_metavr import read_tsv


def main():
    output = ROOT / 'results/metavr'
    memo, checks = {}, []

    def check(path, expected):
        path = Path(path).resolve()
        if path not in memo:
            memo[path] = sha256(path)
        assert memo[path] == expected, ('Checksum mismatch', str(path))
        checks.append(str(path.relative_to(ROOT)))

    bulk_path = ROOT / 'results/scale_integrity_validation.json'
    bulk = json.loads(bulk_path.read_text())
    assert bulk['status'] == 'passed'
    scan = json.loads((output / 'scan.json').read_text())
    check(ROOT / 'scripts/scan_bgzf.py', scan['identity']['script_sha256'])
    for name, digest in scan['summary']['outputs'].items():
        check(output / name, digest)
    runs = list(output.rglob('*.run.json'))
    for path in sorted(runs):
        run = json.loads(path.read_text())
        assert run['status'] == 'complete', path
        check(ROOT / 'scripts/search.py', run['script_sha256'])
        for filename, digest in run['inputs'].items():
            check(filename, digest)
        check(path.with_name(path.name.replace('.run.json', '.tsv')), run['output_sha256'])

    proteins = {r.id: str(r.seq).rstrip('*') for r in SeqIO.parse(output / 'candidates.faa', 'fasta')}
    context = output / 'context'
    ranges = json.loads((context / 'ranges.json').read_text())
    assert ranges['status'] == 'complete'
    check(ROOT / 'scripts/metavr_context.py', ranges['script_sha256'])
    check(output / 'candidates.faa', ranges['proteins_sha256'])
    check(context / 'loci.json', ranges['loci_sha256'])
    mapped = set()
    loci = json.loads((context / 'loci.json').read_text())
    for row in loci:
        genome_path = context / (row['genome'] + '.fna')
        metadata_path = context / (row['genome'] + '.json')
        check(genome_path, row['genome_sha256'])
        check(metadata_path, row['metadata_sha256'])
        genome = str(SeqIO.read(genome_path, 'fasta').seq)
        metadata = json.loads(metadata_path.read_text())
        assert len(genome) == int(metadata['uvig_length'])
        for protein in row['proteins']:
            matches = locate(proteins[protein['id']], genome, int(metadata['genetic_code']))
            assert matches and matches == protein['locations']
            mapped.add(protein['id'])
    assert mapped == set(proteins)
    for filename in ['mouse_identity.json', 'human_identity.json']:
        evidence = json.loads((output / 'host_comparison' / filename).read_text())
        for path, digest in evidence['input_sha256'].items():
            check(path, digest)
    neighbor_dir = output / 'soil_neighbors'
    neighbor_ranges = json.loads((neighbor_dir / 'ranges.json').read_text())
    assert neighbor_ranges['status'] == 'complete'
    for manifest in [ranges, neighbor_ranges]:
        for item in manifest['ranges'].values():
            check(item['path'], item['sha256'])
            assert Path(item['path']).stat().st_size == item['bytes']
    check(ROOT / 'scripts/extract_metavr_neighbors.py', neighbor_ranges['script_sha256'])
    check(neighbor_dir / 'proteins.faa', neighbor_ranges['protein_fasta_sha256'])
    check(neighbor_dir / 'locations.json', neighbor_ranges['locations_sha256'])
    genome = str(SeqIO.read(context / (neighbor_ranges['genome'] + '.fna'), 'fasta').seq)
    check(context / (neighbor_ranges['genome'] + '.fna'), neighbor_ranges['genome_sha256'])
    neighbors = {r.id: str(r.seq).rstrip('*') for r in SeqIO.parse(neighbor_dir / 'proteins.faa', 'fasta')}
    locations = json.loads((neighbor_dir / 'locations.json').read_text())
    for row in locations:
        assert row['locations'] == locate(neighbors[row['protein']], genome, 11)
        assert row['locations']
    assert len(locations) == len(neighbors) == 27
    neighbor_summary = json.loads((neighbor_dir / 'summary.json').read_text())
    check(ROOT / 'scripts/summarize_metavr_neighbors.py', neighbor_summary['script_sha256'])
    check(neighbor_dir / 'gene_review.tsv', neighbor_summary['gene_review_sha256'])
    for path, digest in neighbor_summary['input_sha256'].items():
        check(path, digest)
    summary = json.loads((output / 'summary.json').read_text())
    check(output / 'candidate_review.tsv', summary['candidate_table_sha256'])
    check(ROOT / 'data/raw/metavr/IMG_full_metadata.tsv.gz', summary['source_metadata_sha256'])
    assessment = json.loads((output / 'assessment_summary.json').read_text())
    check(output / 'candidate_assessments.tsv', assessment['candidate_table_sha256'])
    for path, digest in assessment['prior_corpus_inputs'].items():
        check(path, digest)
    reviews = read_tsv(output / 'candidate_assessments.tsv')
    assert {r['protein'] for r in reviews} == set(proteins)
    assert dict(collections.Counter(r['assessment'] for r in reviews)) == assessment['assessments']
    unresolved = [r for r in reviews if r['assessment'] == 'unresolved_environmental_IgV_like']
    assert len(unresolved) == assessment['unresolved_environmental_records'] == 9
    assert len({proteins[r['protein']] for r in unresolved}) == assessment['unresolved_exact_distinct_sequences'] == 7
    report = {'completed_utc': utc(), 'status': 'passed', 'hash_checks': len(checks),
              'distinct_files_hashed': len(memo), 'completed_search_runs': len(runs),
              'candidate_proteins_mapped': len(mapped), 'candidate_genomes': len(loci),
              'neighbor_proteins_mapped': len(locations), 'assessment_counts_verified': True,
              'previous_bulk_validation_sha256': sha256(bulk_path),
              'previous_bulk_hash_checks': bulk['hash_checks'],
              'checks': checks,
              'scope': 'Post-scan artifact integrity and translated mappings. Bulk archive checks reused from the separate completed validation; no inference of antibody/TCR ancestry.'}
    (ROOT / 'results/scale_review_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print('passed', len(checks), 'hash checks;', len(runs), 'completed search runs;', len(mapped), 'candidate mappings')


if __name__ == '__main__':
    main()
