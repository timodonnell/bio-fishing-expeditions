#!/usr/bin/env python3
"""Count/decompress available INPHARED files and measure exact RefSeq overlap."""
import argparse
import gzip
import hashlib
import json
import shutil

from Bio import SeqIO

from acquire import ROOT, sha256, utc


def main():
    parser=argparse.ArgumentParser();parser.add_argument('kind',choices=['protein','genome']);args=parser.parse_args()
    cache=ROOT/'data/cache/inphared';cache.mkdir(parents=True,exist_ok=True)
    name='vConTACT2_proteins.faa' if args.kind=='protein' else 'genomes_excluding_refseq.fa'
    reference=ROOT/('data/cache/refseq/all_proteins.faa' if args.kind=='protein' else 'data/cache/refseq/all_genomes.fna')
    source=ROOT/'data/raw/inphared'/(name+'.gz');target=cache/name
    digest=lambda sequence:hashlib.sha256(str(sequence).upper().encode()).hexdigest()
    known={digest(r.seq) for r in SeqIO.parse(reference,'fasta')}
    seen,identifiers=set(),set();n=total=overlap=0
    with gzip.open(source,'rt') as inp,target.open('w') as out:
        for record in SeqIO.parse(inp,'fasta'):
            assert record.id not in identifiers,record.id
            identifiers.add(record.id);SeqIO.write(record,out,'fasta')
            checksum=digest(record.seq);overlap+=checksum in known;seen.add(checksum)
            n+=1;total+=len(record)
    if args.kind=='genome':
        (ROOT/'manifests/inphared_genome_accessions.txt').write_text('\n'.join(sorted(identifiers))+'\n')
        for name in ['data_excluding_refseq.tsv','vConTACT2_gene_to_genome.csv']:
            with gzip.open(ROOT/'data/raw/inphared'/(name+'.gz'),'rb') as inp,(cache/name).open('wb') as out:
                shutil.copyfileobj(inp,out)
    summary={'completed_utc':utc(),'release':'14Apr2025','kind':args.kind,'records':n,'total_length':total,
             'unique_exact_sequences':len(seen),'records_identical_to_RefSeq':overlap,
             'unique_sequences_not_in_RefSeq':len(seen-known),
             'overlap_rule':'SHA-256 of uppercase sequence; genome comparisons do not normalize rotations or reverse complements',
             'source_sha256':sha256(source),'output_sha256':sha256(target)}
    (ROOT/f'results/inphared_{args.kind}_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
