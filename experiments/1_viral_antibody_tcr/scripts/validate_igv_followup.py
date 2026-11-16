#!/usr/bin/env python3
"""Check the candidate follow-up's hashes and plotted measurements against source evidence."""
import csv
import json
from pathlib import Path

from acquire import ROOT, sha256, utc


def main():
    output = ROOT / 'results/metavr/lead_followup'
    checks = []

    def check(path, digest):
        path = Path(path)
        assert sha256(path) == digest, path
        checks.append(str(path.relative_to(ROOT)))

    for name, script in [('analysis_manifest.json', 'analyze_igv_lead.py'),
                         ('figure_manifest.json', 'plot_igv_lead.py')]:
        manifest = json.loads((output / name).read_text())
        assert manifest['status'] == 'complete'
        check(ROOT / 'scripts' / script, manifest['script_sha256'])
        for path, digest in manifest['input_sha256'].items():
            check(path, digest)
        for name, digest in manifest.get('output_sha256', {}).items():
            check(output / name, digest)
        for name, digest in manifest.get('outputs', {}).items():
            check(ROOT / name, digest)
    data = json.loads((output / 'figure_data.json').read_text())
    with (ROOT / 'results/metavr/hits.tsv').open() as handle:
        hit = next(r for r in csv.DictReader(handle, delimiter='\t')
                   if r['target'].startswith(data['candidate_id'] + '|scan_record='))
    assert hit['query'] == data['profile']['name']
    assert int(hit['ali_from']) == data['profile']['query_start']
    assert int(hit['ali_to']) == data['profile']['query_end']
    assert float(hit['domain_iE']) == data['profile']['global_domain_E']
    assert data['reference']['gene_or_allele_in_record'] is None
    columns = data['alignment_columns']
    assert len(columns) == 89 and sum(c['identity'] for c in columns) == 33
    assert sum(c['query_position'] is not None for c in columns) == 87
    assert [(c['query_position'], c['reference_position']) for c in columns
            if c['query_residue'] == c['reference_residue'] == 'C'] == [(41,43),(107,108)]
    assert all(c['reference_region'] != 'Complementarity-determining-3' for c in columns)
    nr = json.loads((output / 'nr/remote_blast.json').read_text())
    check(output / 'query.faa', nr['query_sha256'])
    report = {'completed_utc': utc(), 'status': 'passed', 'hash_checks': len(checks),
              'checks': checks, 'profile_coordinates_match_original_scan': True,
              'identity_count_and_conserved_cysteines_verified': True,
              'reference_CDR3_outside_alignment': True,
              'remote_nr_status_at_validation': nr['status'],
              'scope': 'Artifact and illustration accuracy, not a donor-allele assignment or novelty proof'}
    (output / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print('passed', len(checks), 'hash checks; plotted measurements agree with source evidence')


if __name__ == '__main__':
    main()
