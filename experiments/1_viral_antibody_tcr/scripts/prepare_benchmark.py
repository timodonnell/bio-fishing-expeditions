#!/usr/bin/env python3
"""Prepare separate calibration targets and verify the deposited FeLV translation."""
import csv
import json
import random
from Bio import SeqIO
from acquire import ROOT, sha256, utc


def main():
    cache = ROOT / 'data/cache/benchmark'
    cache.mkdir(parents=True,exist_ok=True)
    rows = []
    records = {}
    refs = ROOT / 'data/cache/references'
    for filename, category in [('heldout_species.faa','heldout_protein'),('heldout_domains.faa','heldout_domain'),('decoys.faa','hard_negative')]:
        for record in SeqIO.parse(refs/filename,'fasta'):
            identifier=category+'__'+record.id
            sequence=str(record.seq)
            records[identifier]=sequence
            rows.append({'id':identifier,'source_id':record.id,'category':category,'length':len(sequence)})
    domains=json.loads((refs/'domains.json').read_text())
    rng=random.Random(20261112)
    for row in domains:
        if row['split'] != 'heldout_species':
            continue
        sequence=row['sequence']
        for length in [30,60,90]:
            if len(sequence)<length:
                continue
            start=(len(sequence)-length)//2
            identifier=f"fragment{length}__{row['id']}"
            records[identifier]=sequence[start:start+length]
            rows.append({'id':identifier,'source_id':row['id'],'category':f'fragment{length}','length':length})
        shuffled=list(sequence); rng.shuffle(shuffled)
        identifier='shuffled__'+row['id']
        records[identifier]=''.join(shuffled)
        rows.append({'id':identifier,'source_id':row['id'],'category':'shuffled','length':len(sequence)})
    source=ROOT/'data/raw/candidates/X05155.gb'
    record=SeqIO.read(source,'genbank')
    feature=next(f for f in record.features if f.type=='CDS')
    translated=str(feature.extract(record.seq).translate(table=1,cds=True))
    uniprot=json.loads((ROOT/'data/raw/candidates/P11364.json').read_text())
    assert translated==feature.qualifiers['translation'][0]==uniprot['sequence']['value']
    identifier='known_natural__P11364'
    records[identifier]=translated
    rows.append({'id':identifier,'source_id':'P11364','category':'known_natural_positive','length':len(translated)})
    (cache/'targets.faa').write_text(''.join(f'>{k}\n{v}\n' for k,v in records.items()))
    (cache/'known_positive.faa').write_text(f'>P11364\n{translated}\n')
    with (cache/'known_positive.fna').open('w') as handle:
        SeqIO.write(record,handle,'fasta')
        reverse=record.reverse_complement(id='X05155.1_reverse_control',description='Software orientation control, not a second biological record')
        SeqIO.write(reverse,handle,'fasta')
    with (ROOT/'manifests/benchmark.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=['id','source_id','category','length'],delimiter='\t')
        writer.writeheader();writer.writerows(rows)
    evidence={'accession':record.id,'genome_record_bp':len(record),'protein_accession':feature.qualifiers['protein_id'][0],
              'uniprot':'P11364','protein_aa':len(translated),'cds_start_1based':int(feature.location.start)+1,
              'cds_end_1based_inclusive':int(feature.location.end),'strand':feature.location.strand,
              'translation_matches_GenBank_and_UniProt':True,'genbank_sha256':sha256(source),
              'sequence_sha256':__import__('hashlib').sha256(translated.encode()).hexdigest(),
              'features':[{'type':f.type,'location':str(f.location),'notes':f.qualifiers.get('note',[])} for f in record.features],
              'validated_utc':utc(),'claim_scope':'Natural TCR beta chain capture supported by original literature; partial viral nucleotide record, not a whole viral genome; no new raw-read validation'}
    (ROOT/'results/known_positive_validation.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print('Prepared',len(records),'benchmark targets; verified 321-aa FeLV TCR translation')


if __name__=='__main__':
    main()
