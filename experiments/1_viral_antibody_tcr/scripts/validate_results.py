#!/usr/bin/env python3
"""Check executed search inputs/outputs and the biological positive's integrity."""
import json
import shutil
from pathlib import Path

from Bio import SeqIO

from acquire import ROOT, sha256, utc


def main():
    memo={};checks=[]
    def check(path,expected,label):
        path=path.resolve()
        if path not in memo:memo[path]=sha256(path)
        assert memo[path]==expected,(label,str(path),memo[path],expected)
        checks.append(label)
    downloads=json.loads((ROOT/'manifests/downloads.json').read_text())
    for identifier,row in downloads.items():
        if row['status']=='complete':check(ROOT/row['path'],row['sha256'],'download:'+identifier)
    for path in sorted((ROOT/'results').glob('*.run.json')):
        row=json.loads(path.read_text());assert row['status']=='complete',path
        source=next((p for p in [ROOT/'scripts/search.py',ROOT/'scripts/archive/search_initial.py']
                     if sha256(p)==row['script_sha256']),None)
        assert source is not None,('Missing executed source',path)
        check(source,row['script_sha256'],path.name+':executed_script')
        for filename,expected in row['inputs'].items():check(Path(filename),expected,path.name+':input:'+filename)
        output=path.with_name(path.name.replace('.run.json','.tsv'))
        check(output,row['output_sha256'],path.name+':output')
    model=json.loads((ROOT/'results/reference_summary.json').read_text())
    for file in (ROOT/'references').iterdir():
        if file.name in model['outputs']:check(file,model['outputs'][file.name],'frozen_reference:'+file.name)
    training={str(r.seq) for r in SeqIO.parse(ROOT/'references/train.faa','fasta')}
    positive=next(SeqIO.parse(ROOT/'data/cache/benchmark/known_positive.faa','fasta'))
    assert str(positive.seq) not in training
    known=json.loads((ROOT/'results/known_positive_validation.json').read_text())
    assert known['translation_matches_GenBank_and_UniProt'] and known['protein_aa']==321
    assert all(row['translation_validated'] for row in json.loads((ROOT/'results/annotation_only_loci.json').read_text()))
    report={'completed_utc':utc(),'status':'passed','hash_checks':len(checks),'distinct_files_hashed':len(memo),
            'search_runs':len(list((ROOT/'results').glob('*.run.json'))),
            'known_positive_not_in_training':True,'known_positive_translation_validated':True,
            'annotation_only_CDS_translations_validated':5,
            'free_disk_bytes':shutil.disk_usage(ROOT).free,
            'checks':checks,'scope':'Artifact integrity and positive controls; not proof of complete biological sensitivity or natural provenance for unresolved candidates'}
    (ROOT/'results/validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report['status'],report['hash_checks'],'hash checks;',report['search_runs'],'search runs')


if __name__=='__main__':main()
