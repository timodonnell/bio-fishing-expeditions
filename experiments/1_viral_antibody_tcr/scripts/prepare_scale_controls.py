#!/usr/bin/env python3
"""Add two published single-domain antibody constructs as held-out spot checks."""
import json
from pathlib import Path

from Bio import SeqIO

from acquire import ROOT, sha256, utc


def main():
    source=ROOT/'data/raw/references'
    training={str(r.seq) for r in SeqIO.parse(ROOT/'references/train.faa','fasta')}
    controls=[]
    fasta=[]
    for pdb,entity,category,taxon in [('1MVF',1,'camelid_VHH_spot_check',9838),
                                      ('1T6V',2,'shark_VNAR_spot_check',7801)]:
        path=source/f'{pdb}_entity{entity}.json'
        entry_path=source/f'{pdb}_entry.json'
        data=json.loads(path.read_text())
        entry=json.loads(entry_path.read_text())
        assert taxon in {r['ncbi_taxonomy_id'] for r in data['rcsb_entity_source_organism']}
        sequence=data['entity_poly']['pdbx_seq_one_letter_code_can'].replace('\n','')
        assert sequence not in training
        identifier=f'{category}__{pdb}_entity{entity}'
        fasta.append(f'>{identifier}\n{sequence}\n')
        controls.append({'id':identifier,'category':category,'pdb':pdb,'entity':entity,
            'length':len(sequence),'source_taxon':taxon,
            'description':data['rcsb_polymer_entity']['pdbx_description'],
            'author_chains':data['rcsb_polymer_entity_container_identifiers']['auth_asym_ids'],
            'citation':entry['rcsb_primary_citation'],
            'source_url':f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{entity}',
            'entity_sha256':sha256(path),'entry_sha256':sha256(entry_path),
            'exact_sequence_in_training':False,
            'scope':'Published recombinant structural construct, full deposited polymer including tags; retrieval spot check, not a natural viral capture or broad family panel'})
    output=ROOT/'references/scale_antibody_controls.faa'
    output.write_text(''.join(fasta))
    benchmark=ROOT/'data/cache/benchmark/scale_targets.faa'
    benchmark.write_bytes((ROOT/'data/cache/benchmark/targets.faa').read_bytes()+output.read_bytes())
    report={'prepared_utc':utc(),'controls':controls,
            'output_sha256':sha256(output),'combined_benchmark_sha256':sha256(benchmark),
            'original_benchmark_sha256':sha256(ROOT/'data/cache/benchmark/targets.faa'),
            'training_unchanged':True,'profiles_unchanged':True}
    (ROOT/'results/scale_antibody_controls.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Prepared',len(controls),'additional held-out structural constructs; training unchanged')


if __name__=='__main__':main()
