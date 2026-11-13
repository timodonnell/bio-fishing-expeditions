#!/usr/bin/env python3
"""Record evidence-based triage and exact candidate overlap with the prior corpora."""
import collections
import csv
import hashlib
import json

from Bio import SeqIO

from acquire import ROOT, sha256, utc
from summarize_metavr import read_tsv, write_tsv


def main():
    output=ROOT/'results/metavr'
    rows=read_tsv(output/'candidate_review.tsv')
    mouse=json.loads((output/'host_comparison/mouse_identity.json').read_text())
    human=json.loads((output/'host_comparison/human_identity.json').read_text())
    assert mouse['direct_sequence_equality_verified'] and human['query_coverage_percent']==100
    for row in rows:
        genome=row['genome']
        if genome==mouse['genome']:
            assessment='host_mouse_Iglc4_locus'
            evidence='Entire 18461-bp assembly exactly equals the reverse complement of GRCm39 chr16:19008356-19026816; candidate overlaps Iglc4.'
        elif genome==human['genome']:
            assessment='host_human_IGH_locus_supported'
            evidence='All 10225 query bases covered by two >=99.413%-identity hg38 IGH-region alignments; duplicated blocks prevent a single exact-copy claim.'
        elif 'histocompatibility' in row['comparison_description'].lower():
            assessment='competing_MHC_family'
            evidence='Strong reviewed-protein homology overlaps the HMM-matched domain; MHC ancestry is distinct from antibody/TCR ancestry.'
        elif 'CMRF35-like molecule 6' in row['comparison_description']:
            assessment='competing_CD300_family'
            evidence='99.138% identity over 116 aa to mouse CD300c/CMRF35-like molecule 6; an Ig-superfamily receptor, not an antibody or TCR chain.'
        else:
            assert row['source']=='Metagenome'
            assessment='unresolved_environmental_IgV_like'
            evidence='Ig-variable-like domain similarity; receptor-specific ancestry and natural viral incorporation are unproven. Viral provenance, broader homology, and phylogeny require follow-up.'
        row['assessment']=assessment
        row['evidence']=evidence
    sequences={r.id:str(r.seq) for r in SeqIO.parse(output/'candidates.faa','fasta')}
    hashes=collections.defaultdict(list)
    for identifier,sequence in sequences.items():
        hashes[hashlib.sha256(sequence.encode()).hexdigest()].append(identifier)
    overlaps=collections.defaultdict(list)
    corpus_paths=[ROOT/'data/cache/refseq/all_proteins.faa',ROOT/'data/cache/inphared/vConTACT2_proteins.faa']
    for path in corpus_paths:
        for record in SeqIO.parse(path,'fasta'):
            digest=hashlib.sha256(str(record.seq).encode()).hexdigest()
            for candidate in hashes.get(digest,[]):
                overlaps[candidate].append(str(path.relative_to(ROOT))+':'+record.id)
    for row in rows:
        row['exact_prior_corpus_matches']=';'.join(overlaps[row['protein']])
    write_tsv(output/'candidate_assessments.tsv',rows)
    unresolved=[r for r in rows if r['assessment']=='unresolved_environmental_IgV_like']
    with (output/'environmental_candidates.faa').open('w') as handle:
        for row in unresolved:
            handle.write(f'>{row["protein"]}\n{sequences[row["protein"]]}\n')
    report={'completed_utc':utc(),'assessments':dict(collections.Counter(r['assessment'] for r in rows)),
        'by_tier':{tier:dict(collections.Counter(r['assessment'] for r in rows if r['tier']==tier)) for tier in ['main','rescue']},
        'unresolved_environmental_records':len(unresolved),
        'unresolved_exact_distinct_sequences':len({sequences[r['protein']] for r in unresolved}),
        'unresolved_with_high_viral_confidence':sum(r['viral_confidence']=='High' for r in unresolved),
        'unresolved_with_viral_hallmark':sum(int(r['viral_hallmarks'])>0 for r in unresolved),
        'candidate_records_with_exact_prior_corpus_matches':sum(bool(overlaps[r['protein']]) for r in rows),
        'new_supported_natural_antibody_captures':0,'new_supported_natural_TCR_captures':0,
        'candidate_table_sha256':sha256(output/'candidate_assessments.tsv'),
        'prior_corpus_inputs':{str(p):sha256(p) for p in corpus_paths},
        'scope':'Triage only. Exact overlap is measured for the 26 retained candidates, not across all MetaVR records. Unresolved candidates remain open.'}
    (output/'assessment_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
