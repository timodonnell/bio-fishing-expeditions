#!/usr/bin/env python3
"""Prepare the full downloaded RefSeq nucleotide corpus and candidate context."""
import csv
import gzip
import json
import time

from Bio import SeqIO

from acquire import ROOT, sha256, utc
from search import BLAST_FIELDS


def main():
    start = time.monotonic()
    cache = ROOT / 'data/cache/refseq'
    proteins = set()
    for stem in ['pilot', 'refseq_all']:
        for row in csv.DictReader((ROOT/f'results/{stem}_hmm.tsv').open(),delimiter='\t'):
            proteins.add(row['target'].split('|')[0])
        for row in csv.DictReader((ROOT/f'results/{stem}_blastp.tsv').open(),delimiter='\t',fieldnames=BLAST_FIELDS):
            proteins.add(row['sseqid'].split('|')[0])
    hits = list(csv.DictReader((ROOT/'results/pilot_tblastn.tsv').open(),delimiter='\t',fieldnames=BLAST_FIELDS))
    nucleotide_hits = {}
    for row in hits:
        nucleotide_hits.setdefault(row['sseqid'],[]).append(row)
    contexts, translated = [], []
    n, bases, missing = 0, 0, []
    accession_list = []
    with gzip.open(ROOT/'data/raw/refseq/viral.1.genomic.gbff.gz','rt') as inp, (cache/'all_genomes.fna').open('w') as out:
        for record in SeqIO.parse(inp,'genbank'):
            if not record.seq.defined:
                missing.append(record.id)
                continue
            SeqIO.write(record,out,'fasta')
            n += 1; bases += len(record); accession_list.append(record.id)
            for index,feature in enumerate(record.features):
                if feature.type != 'CDS': continue
                protein = feature.qualifiers.get('protein_id',[''])[0]
                row = {'protein':protein,'genome':record.id,'organism':record.annotations.get('organism'),
                       'taxonomy':record.annotations.get('taxonomy',[]),'genome_description':record.description,
                       'start_1based':int(feature.location.start)+1,'end_inclusive':int(feature.location.end),
                       'strand':feature.location.strand,'location':str(feature.location),
                       'qualifiers':feature.qualifiers}
                if protein in proteins:
                    contexts.append(row)
                for hit in nucleotide_hits.get(record.id,[]):
                    begin,end=sorted([int(hit['sstart'])-1,int(hit['send'])-1])
                    strand=1 if int(hit['send'])>=int(hit['sstart']) else -1
                    if feature.location.strand==strand and min(end+1,int(feature.location.end))>max(begin,int(feature.location.start)):
                        translated.append(dict(hit,overlapping_protein=protein,product=feature.qualifiers.get('product',[''])[0]))
    (ROOT/'manifests/refseq_genome_accessions.txt').write_text('\n'.join(accession_list)+'\n')
    with (cache/'candidate_proteins.faa').open('w') as out:
        selected = [r for r in SeqIO.parse(cache/'all_proteins.faa','fasta') if r.id in proteins]
        SeqIO.write(selected,out,'fasta')
    (ROOT/'results/refseq_candidate_context.json').write_text(json.dumps(contexts,indent=2)+'\n')
    with (ROOT/'results/pilot_translated_context.tsv').open('w') as out:
        writer=csv.DictWriter(out,fieldnames=BLAST_FIELDS+['overlapping_protein','product'],delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(translated)
    summary={'completed_utc':utc(),'genomes':n,'bases':bases,'excluded_undefined':missing,
             'candidate_proteins':len(selected),'candidate_contexts':len(contexts),
             'known_positive_X05155_in_genome_corpus':'X05155.1' in accession_list,
             'known_positive_P11364_exact_sequence_in_protein_corpus':False,
             'seconds':time.monotonic()-start,
             'outputs':{p.name:sha256(p) for p in [cache/'all_genomes.fna',cache/'candidate_proteins.faa']}}
    known=str(next(SeqIO.parse(ROOT/'data/cache/benchmark/known_positive.faa','fasta')).seq)
    matches=[r.id for r in SeqIO.parse(cache/'all_proteins.faa','fasta') if str(r.seq)==known]
    summary['known_positive_P11364_exact_sequence_in_protein_corpus']=bool(matches)
    summary['known_positive_exact_matches']=matches
    (ROOT/'results/refseq_expansion_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
