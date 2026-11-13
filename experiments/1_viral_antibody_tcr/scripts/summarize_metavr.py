#!/usr/bin/env python3
"""Summarize executed MetaVR calibration and candidate comparisons without calling ancestry."""
import collections
import csv
import gzip
import json
from pathlib import Path

from Bio import SeqIO

from acquire import ROOT, sha256, utc
from search import BLAST_FIELDS


def read_tsv(path, fields=None):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t', fieldnames=fields))


def write_tsv(path, rows):
    if not rows:
        path.write_text('')
        return
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    output=ROOT/'results/metavr'
    scan=json.loads((output/'scan.json').read_text())
    assert scan['status']=='complete'
    n=scan['summary']['protein_records']
    control_run=json.loads((output/'benchmark_hmm.run.json').read_text())
    assert control_run['status']=='complete' and control_run['parameters']['Z']==n
    recovered={r['target'] for r in read_tsv(output/'benchmark_hmm.tsv')
               if float(r['domain_iE'])<=.001 and int(r['ali_to'])-int(r['ali_from'])+1>=50}
    controls=read_tsv(ROOT/'results/calibration_targets.tsv')
    table=[]
    for row in controls:
        table.append({key:row[key] for key in ['id','source_id','category','length','family','near_identity_excluded']} |
                     {'hmm_retrieved_at_metavr_size':row['id'] in recovered})
    additional=json.loads((ROOT/'results/scale_antibody_controls.json').read_text())
    for row in additional['controls']:
        table.append({'id':row['id'],'source_id':row['pdb']+'_entity'+str(row['entity']),
            'category':row['category'],'length':row['length'],'family':row['category'],
            'near_identity_excluded':'False','hmm_retrieved_at_metavr_size':row['id'] in recovered})
    write_tsv(output/'calibration_targets.tsv',table)
    categories=[]
    for category in sorted({r['category'] for r in table}):
        all_rows=[r for r in table if r['category']==category]
        eligible=[r for r in all_rows if not (category.startswith(('heldout','fragment')) and r['near_identity_excluded']=='True')]
        categories.append({'category':category,'total':len(all_rows),'evaluated':len(eligible),
                           'retrieved':sum(r['hmm_retrieved_at_metavr_size'] for r in eligible)})
    assert 'known_natural__P11364' in recovered
    calibration={'completed_utc':utc(),'effective_targets':n,'categories':categories,
                 'method':'HMM only, direct run at observed full-corpus Z; original near-identity exclusions retained',
                 'additional_control_scope':'One published VHH and one VNAR construct; exact training duplicates excluded; these are spot checks, not complete family/divergence panels',
                 'uncompleted_controls':json.loads((ROOT/'results/calibration_summary.json').read_text())['uncompleted_controls'],
                 'inputs':{str(p):sha256(p) for p in [output/'benchmark_hmm.tsv',ROOT/'results/calibration_targets.tsv']}}
    (output/'calibration_summary.json').write_text(json.dumps(calibration,indent=2)+'\n')
    hits=read_tsv(output/'hits.tsv')
    by_protein=collections.defaultdict(list)
    for row in hits:
        by_protein[row['target']].append(row)
    labels={r.id:r.description for r in SeqIO.parse(ROOT/'data/cache/references/swissprot.faa','fasta')}
    comparisons=collections.defaultdict(list)
    if by_protein:
        for row in read_tsv(output/'candidate_swissprot.tsv',BLAST_FIELDS):
            comparisons[row['qseqid']].append(row)
    context_path=output/'context/loci.json'
    context={r['genome']:r for r in json.loads(context_path.read_text())} if context_path.exists() else {}
    taxon_ids={r.split('_')[2] for r in context if r.startswith('IMGVR_UViG_')}
    samples={}
    sample_source=ROOT/'data/raw/metavr/IMG_full_metadata.tsv.gz'
    with gzip.open(sample_source,'rt') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            if row.get('IMG Genome ID') in taxon_ids:
                samples[row['IMG Genome ID']]={k:v for k,v in row.items() if k}
    write_tsv(output/'source_samples.tsv',list(samples.values()))
    review=[]
    for protein,domains in by_protein.items():
        best=min(domains,key=lambda r:float(r['domain_iE']))
        lo,hi=int(best['ali_from']),int(best['ali_to'])
        overlapping=[]
        for row in comparisons[protein]:
            overlap=max(0,min(hi,int(row['qend']))-max(lo,int(row['qstart']))+1)/(hi-lo+1)
            if overlap>=.5:
                overlapping.append((row,overlap))
        match,overlap=max(overlapping,key=lambda pair:float(pair[0]['bitscore'])) if overlapping else ({},0)
        genome=protein.split('|')[0]
        locus=context.get(genome,{})
        metadata_path=output/'context'/f'{genome}.json'
        metadata=json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
        sample=samples.get(str(metadata.get('taxon_oid','')), {})
        locations=next((p['locations'] for p in locus.get('proteins',[]) if p['id']==protein),[])
        review.append({'protein':protein,'genome':genome,'length':best['target_length'],
            'source_protein_identifier':protein.split('|scan_record=')[0],
            'tier':'main' if any(d['tier']=='main' for d in domains) else 'rescue',
            'best_profile':best['query'],'domain_iE':best['domain_iE'],'ali_from':lo,'ali_to':hi,
            'all_profiles':';'.join(sorted({d['query'] for d in domains})),
            'overlapping_comparison':match.get('sseqid',''),'comparison_E':match.get('evalue',''),
            'comparison_bitscore':match.get('bitscore',''),'comparison_pident':match.get('pident',''),
            'comparison_aligned_aa':match.get('length',''),'domain_overlap_fraction':overlap,
            'comparison_description':labels.get(match.get('sseqid',''),''),
            'genome_bp':metadata.get('uvig_length',''),'quality':metadata.get('quality',''),
            'viral_confidence':metadata.get('viral_confidence',''),'viral_hallmarks':metadata.get('n_hallmarks',''),
            'genomad_score':metadata.get('genomad_score',''),'checkv_contamination':metadata.get('checkv_contamination',''),
            'source':metadata.get('source',''),'host_taxonomy':metadata.get('host_taxonomy',''),
            'sample_name':sample.get('Genome Name / Sample Name',''),
            'study':sample.get('Study Name',''),'ecosystem':sample.get('Ecosystem',''),
            'ictv_taxonomy':metadata.get('ictv_taxonomy',''),'exact_translation_locations':json.dumps(locations),
            'assessment':'Requires manual ancestry/provenance review; no automated natural-capture call'})
    write_tsv(output/'candidate_review.tsv',review)
    report={'completed_utc':utc(),'coverage':scan['summary'],'calibration':categories,
            'candidates_with_context':sum(bool(r['genome_bp']) for r in review),
            'candidates_mapping_exactly_to_genome':sum(r['exact_translation_locations']!='[]' for r in review),
            'candidate_table_sha256':sha256(output/'candidate_review.tsv'),
            'source_metadata_sha256':sha256(sample_source),
            'interpretation':'Similarity screen and candidate-specific assembly checks, not evidence by itself of natural acquisition'}
    (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
