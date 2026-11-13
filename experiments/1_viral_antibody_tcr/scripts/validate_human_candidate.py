#!/usr/bin/env python3
"""Measure high-identity human-reference coverage of the heavy-variable candidate."""
import csv
import json

from Bio import SeqIO

from acquire import ROOT, sha256, utc
from search import BLAST_FIELDS


def main():
    output=ROOT/'results/metavr/host_comparison'
    genome='IMGVR_UViG_3300032031_000646'
    query_path=ROOT/'results/metavr/context'/f'{genome}.fna'
    raw_path=ROOT/'data/raw/candidates/hg38_chr14_105500000_107000000.json'
    target_path=ROOT/'data/cache/metavr/human_IGH_region.fna'
    source=json.loads(raw_path.read_text())
    assert source['genome']=='hg38' and source['chrom']=='chr14'
    assert source['start']==105500000 and source['end']==107000000 and len(source['dna'])==1500000
    expected='>hg38_chr14_105500001_107000000\n'+source['dna']+'\n'
    assert target_path.read_text()==expected
    query=SeqIO.read(query_path,'fasta')
    rows=list(csv.DictReader((output/'human_igh_blastn.tsv').open(),delimiter='\t',fieldnames=BLAST_FIELDS))
    eligible=[r for r in rows if float(r['pident'])>=99 and int(r['length'])>=500]
    selected=[]
    end=0
    while end<len(query):
        choices=[r for r in eligible if int(r['qstart'])<=end+1 and int(r['qend'])>end]
        assert choices,'Not all query bases covered'
        row=max(choices,key=lambda r:(int(r['qend']),float(r['bitscore'])))
        selected.append(row)
        end=int(row['qend'])
    report={'status':'host_locus_supported','verified_utc':utc(),'genome':genome,
        'query_length':len(query),'query_coverage_percent':100.0,
        'minimum_selected_alignment_identity_percent':min(float(r['pident']) for r in selected),
        'reference_assembly':'hg38 (GRCh38)','chromosome':'chr14',
        'reference_window_start_0based':source['start'],'reference_window_end_exclusive':source['end'],
        'selected_alignments':[dict(r,chromosome_start_1based=source['start']+min(int(r['sstart']),int(r['send'])),
            chromosome_end_inclusive=source['start']+max(int(r['sstart']),int(r['send']))) for r in selected],
        'input_sha256':{str(p):sha256(p) for p in [query_path,raw_path,target_path,output/'human_igh_blastn.tsv']},
        'assessment':'All assembly bases are covered by >=99% identity alignments to the human IGH reference region. The best alignments involve duplicated blocks, so this is not a single exact contiguous-reference match. Host-locus sequence/assembly is strongly supported; no independent viral flanks or natural capture established.'}
    (output/'human_identity.json').write_text(json.dumps(report,indent=2)+'\n')
    print(len(query),'bases fully covered by',len(selected),'human-reference alignments;',report['minimum_selected_alignment_identity_percent'],'percent minimum identity')


if __name__=='__main__':main()
