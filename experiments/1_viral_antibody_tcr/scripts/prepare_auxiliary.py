#!/usr/bin/env python3
"""Reproduce the small extraction steps originally executed during exploration."""
import argparse
import csv
import gzip
import json

from Bio import SeqIO

from acquire import ROOT, sha256


def main():
    parser=argparse.ArgumentParser();parser.add_argument('kind',choices=['swissprot','annotation','phage_candidates']);args=parser.parse_args()
    if args.kind=='swissprot':
        target=ROOT/'data/cache/references/swissprot.faa';n=length=0
        with gzip.open(ROOT/'data/raw/references/uniprot_reviewed.fasta.gz','rt') as inp,target.open('w') as out:
            for record in SeqIO.parse(inp,'fasta'):
                SeqIO.write(record,out,'fasta');n+=1;length+=len(record)
        print(n,length,sha256(target))
    elif args.kind=='annotation':
        cache=ROOT/'data/cache/genbank';cache.mkdir(parents=True,exist_ok=True)
        records=[]
        for path in sorted((ROOT/'data/raw/genbank').glob('records_*.gb')):records.extend(SeqIO.parse(path,'genbank'))
        SeqIO.write(records,cache/'annotation_selected.faa','fasta')
        with (ROOT/'manifests/genbank_annotation_records.tsv').open('w') as out:
            writer=csv.writer(out,delimiter='\t',lineterminator='\n');writer.writerow(['accession','length','organism','description','dbxrefs'])
            for r in records:writer.writerow([r.id,len(r),r.annotations.get('organism'),r.description,';'.join(r.dbxrefs)])
        profile=ROOT/'results/genbank_annotation_hmm.tsv'
        if profile.exists():
            ids={r['target'] for r in csv.DictReader(profile.open(),delimiter='\t')}
            ids.update(['DAJ54870.1','DAW92981.1','DBA58143.1','XOK31829.1','XYW61038.1'])
            SeqIO.write([r for r in records if r.id in ids],cache/'candidate_proteins.faa','fasta')
        print(len(records),sha256(cache/'annotation_selected.faa'))
    else:
        ids={r[1] for r in csv.reader((ROOT/'results/inphared_blastp.tsv').open(),delimiter='\t')}
        cache=ROOT/'data/cache/inphared'
        count=SeqIO.write((r for r in SeqIO.parse(cache/'vConTACT2_proteins.faa','fasta') if r.id in ids),cache/'candidate_proteins.faa','fasta')
        print(count,sha256(cache/'candidate_proteins.faa'))


if __name__=='__main__':main()
