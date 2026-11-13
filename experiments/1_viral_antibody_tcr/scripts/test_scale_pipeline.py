#!/usr/bin/env python3
"""Integration checks for chunk boundaries, resume integrity, and global HMM E-values."""
import csv
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import bgzf

from acquire import ROOT, sha256, utc
from scan_bgzf import read_shard


def rows(path):
    result = {}
    for row in csv.DictReader(path.open(), delimiter='\t'):
        key = tuple(row[k] for k in ['query', 'target', 'ali_from', 'ali_to', 'hmm_from', 'hmm_to'])
        result[key] = row
    return result


def compare(a, b):
    assert set(a) == set(b), (len(a), len(b), set(a)^set(b))
    for key in a:
        for field in ['sequence_E', 'domain_iE', 'sequence_score', 'domain_score']:
            assert math.isclose(float(a[key][field]), float(b[key][field]), rel_tol=1e-6, abs_tol=1e-100), (key,field)


def main():
    work = Path(tempfile.mkdtemp(prefix='scale_pipeline_test_', dir=ROOT/'work'))
    with bgzf.BgzfReader(str(ROOT/'data/cache/metavr/sample.bgz'), 'rb') as source:
        sample = read_shard(source, 4*2**20)
    control = (ROOT/'data/cache/benchmark/targets.faa').read_bytes()
    # A record crossing multiple BGZF blocks catches boundary handling errors.
    data = control + sample + b'>boundary_long_control\n' + b'A'*90000 + b'\n'
    target = work/'targets.bgz'
    with bgzf.BgzfWriter(str(target), 'wb') as handle:
        handle.write(data)
    (work/'targets.faa').write_bytes(data)
    script = ROOT/'scripts/scan_bgzf.py'
    common = ['--cpus','2','--workers','2','--effective-target-count','1000000000']
    commands = []
    def run(output, extra=(), success=True):
        command = [sys.executable, str(script), str(target), str(work/output), *common, *extra]
        commands.append(command)
        result = subprocess.run(command, capture_output=True, text=True)
        if success and result.returncode:
            raise RuntimeError(result.stdout+result.stderr)
        assert success == (result.returncode == 0), result.stderr
        return result
    run('whole', ['--shard-mib','64'])
    run('sharded', ['--shard-mib','1','--stop-after-shards','1'])
    state = json.loads((work/'sharded/scan.json').read_text())
    assert state['status']=='paused_for_restart_test' and state['completed_shards']==1
    first = work/'sharded/shards/000000.tsv'
    first_digest = sha256(first)
    saved = first.read_bytes()
    first.write_bytes(saved+b'corrupt\n')
    error = run('sharded', ['--shard-mib','1'], success=False)
    assert 'Checkpoint output changed' in error.stderr
    first.write_bytes(saved)
    run('sharded', ['--shard-mib','1'])
    assert sha256(first)==first_digest
    whole, sharded = rows(work/'whole/hits.tsv'), rows(work/'sharded/hits.tsv')
    compare(whole, sharded)
    independent = [sys.executable, str(ROOT/'scripts/search.py'), 'hmm',
                   str(ROOT/'references/receptor_profiles.hmm'), str(work/'targets.faa'),
                   str(work/'independent.tsv'), '--cpus','4','--effective-target-count','1000000000']
    commands.append(independent)
    subprocess.run(independent, check=True)
    compare(whole, rows(work/'independent.tsv'))
    summaries = [json.loads((work/name/'scan.json').read_text())['summary'] for name in ['whole','sharded']]
    for key in ['protein_records','residues','distinct_source_genome_identifiers','retained_domain_rows']:
        assert summaries[0][key]==summaries[1][key], key
    positive = [r for r in whole.values() if r['target']=='known_natural__P11364' and r['tier']=='main']
    assert positive, 'Known natural TCR positive not recovered'
    report = {'status':'passed','completed_utc':utc(),'commands':commands,
              'checks':['multiple BGZF blocks within a record','complete shard-boundary coverage',
                        'restart preserves completed shard','corrupt checkpoint rejected',
                        'sharded/unsharded hit and score equivalence',
                        'global E-values match independent original scanner at Z=1 billion',
                        'known natural FeLV TCR recovered at Z=1 billion'],
              'target_sha256':sha256(target),'script_sha256':sha256(script),
              'test_script_sha256':sha256(Path(__file__)), 'protein_records':summaries[0]['protein_records'],
              'shards_tested':summaries[1]['shards'],'domain_hits_compared':len(whole)}
    (ROOT/'results/scale_pipeline_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
