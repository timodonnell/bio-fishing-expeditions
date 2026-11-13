#!/usr/bin/env python3
"""Recover every retained translated hit, including loci absent from protein candidates."""
import csv
import gzip
import json

from Bio import SeqIO

from acquire import ROOT, sha256, utc
from search import BLAST_FIELDS


def main():
    rows=list(csv.DictReader((ROOT/'results/refseq_all_tblastn.tsv').open(),delimiter='\t',fieldnames=BLAST_FIELDS))
    hits={}
    for row in rows:
        lo,hi=sorted([int(row['sstart']),int(row['send'])]);strand=1 if int(row['send'])>=int(row['sstart']) else -1
        key=(row['sseqid'],lo,hi,strand)
        hits.setdefault(key,[]).append(row)
    records={}
    wanted={key[0] for key in hits}
    with gzip.open(ROOT/'data/raw/refseq/viral.1.genomic.gbff.gz','rt') as handle:
        for record in SeqIO.parse(handle,'genbank'):
            if record.id in wanted:records[record.id]=record
    output=ROOT/'data/cache/refseq/translated_candidates.faa'
    contexts=[]
    with output.open('w') as handle:
        for (accession,lo,hi,strand),alignments in sorted(hits.items()):
            record=records[accession];sequence=record.seq[lo-1:hi]
            if strand==-1:sequence=sequence.reverse_complement()
            assert len(sequence)%3==0
            protein=str(sequence.translate(to_stop=False)).replace('*','X')
            name=f'{accession}_{lo}_{hi}_{strand}'
            handle.write(f'>{name}\n{protein}\n')
            features=[]
            for feature in record.features:
                if feature.type not in ['CDS','gene','misc_feature']:continue
                if feature.location.strand not in [strand,None]:continue
                if min(hi,int(feature.location.end))<=max(lo-1,int(feature.location.start)):continue
                features.append({'type':feature.type,'location':str(feature.location),'qualifiers':feature.qualifiers})
            contexts.append({'id':name,'genome':accession,'organism':record.annotations.get('organism'),
                             'start_1based':lo,'end_inclusive':hi,'strand':strand,
                             'best_E':min(float(r['evalue']) for r in alignments),
                             'main_queue':any(float(r['evalue'])<=.001 and int(r['length'])>=50 for r in alignments),
                             'alignments':alignments,'overlapping_features':features})
    (ROOT/'results/refseq_translated_context.json').write_text(json.dumps(contexts,indent=2)+'\n')
    print(len(contexts),'distinct translated intervals',sha256(output),utc())


if __name__=='__main__':main()
