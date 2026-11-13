#!/usr/bin/env python3
"""Report retrieval, including near-identity exclusions; never call ancestry from scores."""
import collections
import csv
import json

from acquire import ROOT, sha256, utc
from search import BLAST_FIELDS


def main():
    results = ROOT / 'results'
    manifest = {r['id']: r for r in csv.DictReader((ROOT / 'manifests/benchmark.tsv').open(), delimiter='\t')}
    domains = {r['id']: r for r in json.loads((ROOT / 'data/cache/references/domains.json').read_text())}
    references = {r['Entry']: r for r in csv.DictReader((ROOT / 'manifests/references.tsv').open(), delimiter='\t')}
    effective_hmm_sequences = (json.loads((results/'corpus_summary.json').read_text())['downloaded_proteins']
                               + json.loads((results/'inphared_protein_summary.json').read_text())['records'])
    hmm_scale = effective_hmm_sequences / 1000000
    # Sensitive, small-database BLAST results are used only to identify leakage.
    near = set()
    for row in csv.DictReader((results / 'benchmark_blastp.tsv').open(), delimiter='\t', fieldnames=BLAST_FIELDS):
        coverage = (abs(int(row['send']) - int(row['sstart'])) + 1) / int(row['slen'])
        if float(row['pident']) >= 90 and coverage >= 0.8:
            near.add(row['sseqid'])
    hmm, blast = set(), set()
    for row in csv.DictReader((results / 'benchmark_hmm_global.tsv').open(), delimiter='\t'):
        if float(row['domain_iE']) * hmm_scale <= 1e-3 and int(row['ali_to']) - int(row['ali_from']) + 1 >= 50:
            hmm.add(row['target'])
    for row in csv.DictReader((results / 'benchmark_blastp_global.tsv').open(), delimiter='\t', fieldnames=BLAST_FIELDS):
        if float(row['evalue']) <= 1e-3 and int(row['length']) >= 50:
            blast.add(row['sseqid'])
    table = []
    for name, row in manifest.items():
        source, category = row['source_id'], row['category']
        accession = domains.get(source, {}).get('accession', source)
        family = domains.get(source, {}).get('family', references.get(accession, {}).get('family', ''))
        # Exclude fragments if their source domain leaks, even if the fragment alignment does not.
        leakage = name in near or ('heldout_domain__' + source) in near
        table.append(dict(row, family=family, near_identity_excluded=leakage,
                          hmm_retrieved=name in hmm, blastp_retrieved=name in blast,
                          union_retrieved=name in hmm | blast))
    with (results / 'calibration_targets.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(table)
    summary = []
    for category in sorted({r['category'] for r in table}):
        all_rows = [r for r in table if r['category'] == category]
        rows = [r for r in all_rows if not r['near_identity_excluded']] if category.startswith(('heldout', 'fragment')) else all_rows
        summary.append({'category':category, 'total':len(all_rows), 'evaluated':len(rows),
                        'excluded':len(all_rows)-len(rows),
                        **{method:sum(r[method+'_retrieved'] for r in rows) for method in ['hmm','blastp','union']}})
    by_family = []
    for family in sorted({r['family'] for r in table if r['category']=='heldout_domain'}):
        rows = [r for r in table if r['category']=='heldout_domain' and r['family']==family and not r['near_identity_excluded']]
        by_family.append({'family':family,'evaluated':len(rows),'hmm':sum(r['hmm_retrieved'] for r in rows),
                          'blastp':sum(r['blastp_retrieved'] for r in rows)})
    result = {'completed_utc':utc(),'main_threshold':{'E_or_independent_domain_E':0.001,'aligned_aa':50},
              'effective_HMM_sequences':effective_hmm_sequences,'effective_BLAST_residues':1000000000,
              'HMM_rescaling':'Independent domain E-values from Z=1,000,000 multiplied by combined RefSeq+INPHARED protein record count / 1,000,000; overlapping records retained conservatively',
              'near_identity_rule':'BLASTP >=90% identity over >=80% target span; fragments also inherit source-domain exclusion',
              'categories':summary,'heldout_domains_by_family':by_family,
              'interpretation':'Retrieval only. Hard-negative matches remain decoys, not discoveries. Broad negative-result sensitivity is not established.',
              'uncompleted_controls':['verified engineered carrier genomes','camelid VHH','shark IgNAR',
                                      'systematic frameshift/splice cases','divergence series','family-holdout calibration'],
              'input_sha256':{p.name:sha256(p) for p in [results/'benchmark_blastp.tsv',results/'benchmark_hmm_global.tsv',results/'benchmark_blastp_global.tsv']}}
    (results/'calibration_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
