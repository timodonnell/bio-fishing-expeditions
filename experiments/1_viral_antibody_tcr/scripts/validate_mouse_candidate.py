#!/usr/bin/env python3
"""Verify whole-contig host identity independently of the BLAST alignment."""
import csv
import json

from Bio import SeqIO

from acquire import ROOT, sha256, utc
from search import BLAST_FIELDS


def main():
    output=ROOT/'results/metavr/host_comparison'
    genome='IMGVR_UViG_3300051405_000263'
    query_path=ROOT/'results/metavr/context'/f'{genome}.fna'
    reference_path=ROOT/'data/raw/candidates/GRCm39_chr16_18850000_19050000.fna'
    annotation_path=ROOT/'data/raw/candidates/mouse_candidate_overlap.json'
    query=SeqIO.read(query_path,'fasta')
    reference=SeqIO.read(reference_path,'fasta')
    assert reference.id=='chromosome:GRCm39:16:18850000:19050000:1'
    rows=list(csv.DictReader((output/'mouse_igl_blastn.tsv').open(),delimiter='\t',fieldnames=BLAST_FIELDS))
    best=max(rows,key=lambda r:float(r['bitscore']))
    assert float(best['pident'])==100 and int(best['length'])==len(query)==18461
    assert int(best['qstart'])==1 and int(best['qend'])==len(query)
    start,end=sorted([int(best['sstart']),int(best['send'])])
    matched=reference.seq[start-1:end].reverse_complement()
    assert str(query.seq).upper()==str(matched).upper()
    locus=next(r for r in json.loads((ROOT/'results/metavr/context/loci.json').read_text()) if r['genome']==genome)
    location=locus['proteins'][0]['locations'][0]
    chrom_end=18850000+end-1
    protein_start=chrom_end-location['end_inclusive']+1
    protein_end=chrom_end-location['start_1based']+1
    genes=json.loads(annotation_path.read_text())
    overlapping=[g for g in genes if g['start']<=protein_end and g['end']>=protein_start]
    assert any(g['external_name']=='Iglc4' for g in overlapping)
    report={'status':'verified_host_locus','verified_utc':utc(),'genome':genome,
        'query_length':len(query),'whole_contig_identity_percent':100.0,'coverage_percent':100.0,
        'reference_assembly':'GRCm39','chromosome':'16','reference_start_1based':18850000+start-1,
        'reference_end_inclusive':chrom_end,'alignment_strand':-1,
        'protein_reference_start_1based':protein_start,'protein_reference_end_inclusive':protein_end,
        'overlapping_host_genes':overlapping,'direct_sequence_equality_verified':True,
        'input_sha256':{str(p):sha256(p) for p in [query_path,reference_path,annotation_path,output/'mouse_igl_blastn.tsv']},
        'assessment':'The complete catalog assembly is exact mouse immunoglobulin-locus DNA, including Iglc4. No evidence here establishes viral genomic incorporation; treat as host sequence misclassified in the viral catalog.'}
    (output/'mouse_identity.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report['assessment'],report['reference_start_1based'],report['reference_end_inclusive'])


if __name__=='__main__':main()
