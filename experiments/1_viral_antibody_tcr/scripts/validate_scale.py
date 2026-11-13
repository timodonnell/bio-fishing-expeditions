#!/usr/bin/env python3
"""Check full MetaVR coverage, artifact integrity, global scoring, and controls."""
import csv
import hashlib
import json
import math
from pathlib import Path

from Bio import SeqIO

from acquire import ROOT, sha256, utc


def main():
    output = ROOT/'results/metavr'
    scan = json.loads((output/'scan.json').read_text())
    download = json.loads((ROOT/'manifests/metavr_download.json').read_text())
    checks=[]
    memo={}
    def check(path, expected):
        path=Path(path).resolve()
        if path not in memo:
            memo[path]=sha256(path)
        if memo[path]!=expected:
            raise ValueError('Checksum mismatch: '+str(path))
        checks.append(str(path.relative_to(ROOT)))
    assert scan['status']==download['status']=='complete'
    assert scan['identity']['target_sha256']==download['sha256']
    check(ROOT/download['path'],download['sha256'])
    check(ROOT/'scripts/scan_bgzf.py',scan['identity']['script_sha256'])
    check(ROOT/'references/receptor_profiles.hmm',scan['identity']['profiles_sha256'])
    assert download['bytes']==scan['identity']['target_bytes']
    assert sum(r['bytes'] for r in download['completed_ranges'].values())==download['bytes']
    chunk=download['config']['range_bytes']
    assert set(map(int,download['completed_ranges']))==set(range((download['bytes']+chunk-1)//chunk))
    checkpoints=[]
    raw={}
    for path in sorted((output/'shards').glob('*.json')):
        row=json.loads(path.read_text())
        assert row['shard']==len(checkpoints)
        assert row['start_virtual_offset']==(checkpoints[-1]['end_virtual_offset'] if checkpoints else 0)
        assert row['first_record_ordinal']==sum(c['records'] for c in checkpoints)
        for filename,expected in row['outputs'].items():
            check(filename,expected)
        with path.with_suffix('.tsv').open() as handle:
            for hit in csv.DictReader(handle,delimiter='\t'):
                key=(row['shard'],hit['query'],hit['target'],hit['ali_from'],hit['ali_to'],hit['hmm_from'],hit['hmm_to'])
                assert key not in raw
                raw[key]=hit
        checkpoints.append(row)
    summary=scan['summary']
    assert len(checkpoints)==summary['shards']
    assert checkpoints[-1]['end_virtual_offset'] >> 16 == download['bytes']
    assert sum(r['records'] for r in checkpoints)==summary['protein_records']
    assert sum(r['residues'] for r in checkpoints)==summary['residues']
    assert sum(r['headers_without_pipe'] for r in checkpoints)==summary['headers_without_pipe']==0
    assert summary['effective_target_count']==summary['protein_records']
    assert len(raw)==summary['raw_domain_rows']
    z=summary['effective_target_count']/scan['identity']['parameters']['reference_Z']
    retained={};main_ids=set()
    with (output/'hits.tsv').open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            key=(int(row['shard']),row['query'],row['target'],row['ali_from'],row['ali_to'],row['hmm_from'],row['hmm_to'])
            original=raw[key]
            for field in ['sequence_E','domain_iE']:
                assert math.isclose(float(row[field]),float(original[field])*z,rel_tol=1e-12)
            assert float(row['domain_iE'])<=0.01
            expected='main' if float(row['domain_iE'])<=.001 and int(row['ali_to'])-int(row['ali_from'])+1>=50 else 'rescue'
            assert row['tier']==expected
            if expected=='main':main_ids.add(row['target'])
            retained[key]=row
    expected={key for key,row in raw.items() if float(row['domain_iE'])*z<=.01}
    assert set(retained)==expected
    proteins={r.id:str(r.seq) for r in SeqIO.parse(output/'candidates.faa','fasta')}
    assert set(proteins)=={r['target'] for r in retained.values()}
    assert len(proteins)==summary['retained_candidate_proteins']
    assert len(main_ids)==summary['main_candidate_proteins']
    for name,digest in summary['outputs'].items():check(output/name,digest)
    for report_name,script in [('scale_pipeline_validation.json','scan_bgzf.py'),
                               ('scale_download_validation.json','acquire_large.py'),
                               ('scale_context_validation.json','metavr_context.py')]:
        control=json.loads((ROOT/'results'/report_name).read_text())
        assert control['status']=='passed'
        check(ROOT/'scripts'/script,control['script_sha256'])
    report={'status':'passed','completed_utc':utc(),'hash_checks':len(checks),
            'distinct_files_hashed':len(memo),'protein_records':summary['protein_records'],
            'shards':len(checkpoints),'contiguous_BGZF_coverage_to_EOF':True,
            'global_E_value_rescaling_verified':True,'all_retained_sequences_present':True,
            'pipeline_controls_passed':True,'checks':checks,
            'scope':'Artifact and coverage integrity; not a biological ancestry verdict or complete sensitivity assessment'}
    (ROOT/'results/scale_integrity_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report['status'],report['hash_checks'],'hash checks;',report['protein_records'],'proteins')


if __name__=='__main__':main()
