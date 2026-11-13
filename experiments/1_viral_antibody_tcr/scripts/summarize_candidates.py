#!/usr/bin/env python3
"""Build transparent candidate tables; weak/annotation-only matches stay unresolved."""
import collections
import csv
import json
import re

from Bio import SeqIO

from acquire import ROOT, sha256, utc
from search import BLAST_FIELDS


def blast(path):
    return list(csv.DictReader(path.open(),delimiter='\t',fieldnames=BLAST_FIELDS))


def table(path, rows):
    if not rows:return
    with path.open('w') as out:
        writer=csv.DictWriter(out,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def main():
    results=ROOT/'results'
    labels={r.id:r.description for r in SeqIO.parse(ROOT/'data/cache/references/swissprot.faa','fasta')}
    def competing(filename):
        best={}
        for row in sorted(blast(results/filename),key=lambda r:-float(r['bitscore'])):best.setdefault(row['qseqid'],row)
        return best
    def describe(best):
        if best is None:return {'comparison_accession':'','comparison_E':'','comparison_bitscore':'','comparison_description':''}
        return {'comparison_accession':best['sseqid'],'comparison_E':best['evalue'],
                'comparison_bitscore':best['bitscore'],'comparison_description':labels[best['sseqid']]}
    def interpretation(description,evalue,annotation):
        text=description.lower()
        strong=bool(evalue!='' and float(evalue)<=1e-5)
        if strong and any(x in text for x in ['histocompatibility','hemochromatosis','zinc-alpha-2']):
            return 'competing_MHC_family_supported'
        if strong and any(x in text for x in ['ox-2 membrane','interleukin-1-binding','ig-like','outer capsid','protein a56','protein opg185','envelope glycoprotein c']):
            return 'other_Ig_like_or_viral_surface_protein_supported'
        if strong and any(x in text for x in ['polymerase','exonuclease','ubiquitin','capsid polyprotein']):
            return 'other_protein_family_supported'
        if any(x in annotation.lower() for x in ['mhc','ox-2','ig-like','ig domain','immunoglobulin domain']):
            return 'other_Ig_like_annotation_not_receptor_specific'
        return 'unresolved_weak_or_nonspecific_match'
    main_ids=set();route_counts=[]
    for corpus in ['pilot','refseq_all','inphared']:
        for method in ['hmm','blastp','tblastn']:
            path=results/f'{corpus}_{method}.tsv'
            if method=='hmm':
                rows=list(csv.DictReader(path.open(),delimiter='\t'))
                retained={r['target'] for r in rows}
                selected={r['target'] for r in rows if float(r['domain_iE'])<=.001 and int(r['ali_to'])-int(r['ali_from'])+1>=50}
            else:
                rows=blast(path);retained={r['sseqid'] for r in rows}
                selected={r['sseqid'] for r in rows if float(r['evalue'])<=.001 and int(r['length'])>=50}
            route_counts.append({'corpus':corpus,'method':method,'alignments_retained':len(rows),
                                 'targets_retained':len(retained),'main_queue_targets':len(selected)})
            if method!='tblastn' and corpus!='inphared':main_ids.update(x.split('|')[0] for x in selected)
    best=competing('refseq_competing_blastp.tsv');protein_rows=[]
    for row in json.loads((results/'refseq_candidate_context.json').read_text()):
        desc=describe(best.get(row['protein']));annotation='; '.join(row['qualifiers'].get('product',[])+row['qualifiers'].get('note',[]))
        protein_rows.append({'protein':row['protein'],'genome':row['genome'],'organism':row['organism'],
                             'start_1based':row['start_1based'],'end_inclusive':row['end_inclusive'],'strand':row['strand'],
                             'main_queue_in_pilot_or_full':row['protein'] in main_ids,'annotation':annotation,
                             **desc,'assessment':interpretation(desc['comparison_description'],desc['comparison_E'],annotation)})
    table(results/'refseq_protein_candidates.tsv',protein_rows)
    best=competing('refseq_translated_competing_blastp.tsv');translated_rows=[]
    for row in json.loads((results/'refseq_translated_context.json').read_text()):
        desc=describe(best.get(row['id']));annotation='; '.join(x for f in row['overlapping_features'] for k in ['product','note'] for x in f['qualifiers'].get(k,[]))
        translated_rows.append({'id':row['id'],'genome':row['genome'],'organism':row['organism'],
                                'start_1based':row['start_1based'],'end_inclusive':row['end_inclusive'],'strand':row['strand'],
                                'main_queue':row['main_queue'],'best_E':row['best_E'],
                                'overlapping_CDS':';'.join(x for f in row['overlapping_features'] if f['type']=='CDS' for x in f['qualifiers'].get('protein_id',[])),
                                **desc,'assessment':interpretation(desc['comparison_description'],desc['comparison_E'],annotation)})
    table(results/'refseq_translated_candidates.tsv',translated_rows)
    # Annotation-only GenBank leads are retained even when both profile and reviewed-protein searches fail.
    ids={r.id for r in SeqIO.parse(ROOT/'data/cache/genbank/candidate_proteins.faa','fasta')}
    known={'CAA28801.1','P11364.1','1304304A'}
    best=competing('genbank_competing_blastp.tsv');annotation_rows=[];loci=[]
    for path in sorted((ROOT/'data/raw/genbank').glob('records_*.gb')):
        for record in SeqIO.parse(path,'genbank'):
            if record.id not in ids:continue
            source=record.annotations.get('db_source','');coded_by=';'.join(x for f in record.features if f.type=='CDS' for x in f.qualifiers.get('coded_by',[]))
            if record.id in known:assessment='known_natural_FeLV_TCR_single_event'
            elif source.startswith('pdb:'):assessment='structure_component_no_viral_genomic_evidence'
            elif source.startswith('prf:'):assessment='legacy_protein_record_unresolved_taxonomic_provenance'
            else:assessment='annotation_only_unresolved'
            annotation_rows.append({'accession':record.id,'length':len(record),'description':record.description,
                                    'organism':record.annotations.get('organism',''),'record_source':source,'coded_by':coded_by,
                                    **describe(best.get(record.id)),'assessment':assessment})
            if assessment=='annotation_only_unresolved':
                accession=re.search(r'([A-Z]+\d+\.\d+):',coded_by)[1]
                genome=SeqIO.read(ROOT/f'data/raw/candidates/{accession}.gb','genbank')
                cds=next(f for f in genome.features if record.id in f.qualifiers.get('protein_id',[]))
                translated=str(cds.extract(genome.seq).translate(table=int(cds.qualifiers.get('transl_table',['1'])[0]),to_stop=True))
                assert translated==str(record.seq),(record.id,translated,str(record.seq))
                loci.append({'protein':record.id,'genome':accession,'genome_bp':len(genome),'protein_aa':len(record),
                             'location':str(cds.location),'translation_validated':True,'qualifiers':cds.qualifiers,
                             'genome_source':next(f.qualifiers for f in genome.features if f.type=='source'),
                             'genbank_sha256':sha256(ROOT/f'data/raw/candidates/{accession}.gb'),
                             'assessment':'Annotation-only unresolved; no convincing receptor-specific match; no raw-read validation'})
    table(results/'genbank_annotation_candidates.tsv',annotation_rows)
    (results/'annotation_only_loci.json').write_text(json.dumps(loci,indent=2)+'\n')
    table(results/'search_route_counts.tsv',route_counts)
    summary={'completed_utc':utc(),'refseq_protein_candidates':len(protein_rows),
             'refseq_protein_assessments':dict(collections.Counter(r['assessment'] for r in protein_rows)),
             'refseq_translated_intervals':len(translated_rows),'refseq_translated_main_intervals':sum(r['main_queue'] for r in translated_rows),
             'refseq_translated_assessments':dict(collections.Counter(r['assessment'] for r in translated_rows)),
             'genbank_annotation_candidates':len(annotation_rows),'genbank_annotation_assessments':dict(collections.Counter(r['assessment'] for r in annotation_rows)),
             'caution':'Assessments are triage, not exhaustive evolutionary proof. Unresolved entries remain open. Multiple accession representations of FeLV are one biological precedent.'}
    (results/'candidate_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
