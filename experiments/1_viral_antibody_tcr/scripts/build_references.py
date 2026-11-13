#!/usr/bin/env python3
"""Build host-only sequence references and pilot HMMs; keep FeLV out of training."""
import collections
import csv
import json
import re
import subprocess

import pyhmmer
from Bio import SeqIO

from acquire import ROOT, sha256, utc

MUSCLE = '/home/bizon/git/ANARCI/bin/muscle'


def classify(name):
    name = name.lower()
    if 'lambda-like' in name or 'alternate reading frame' in name:
        return None
    match = re.match(r'(?:probable non-functional )?t[- ]cell receptor (alpha|beta|gamma|delta)', name)
    if match:
        family = {'alpha':'TRA', 'beta':'TRB', 'gamma':'TRG', 'delta':'TRD'}[match[1]]
    else:
        match = re.match(r'(?:immunoglobulin|ig) (heavy|kappa|lambda|gamma|mu|epsilon|alpha|delta)', name)
        if not match:
            return None
        family = {'kappa':'IGK','lambda':'IGL'}.get(match[1], 'IGH')
    if re.search(r'variable|\bv(?:-[a-z0-9]+)? region', name):
        return family + '_V'
    if re.search(r'constant|\bc region', name):
        return family + '_C'
    return family + '_chain'


def main():
    cache = ROOT / 'data/cache/references'
    cache.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader((ROOT / 'data/raw/references/uniprot_expanded.tsv').open(), delimiter='\t'))
    references, domains, decoys, exclusions = [], [], [], []
    for row in sorted(rows, key=lambda r:r['Entry']):
        family = classify(row['Protein names'])
        if not family:
            if any(term in row['Protein names'].lower() for term in
                   ['immunoglobulin-like', 'immunoglobulin superfamily', 'fc receptor',
                    'fc region receptor', 'immunoglobulin g-binding', 'lambda-like']):
                decoys.append(row)
            exclusions.append({'accession':row['Entry'], 'reason':'not_unambiguously_IG_TR',
                               'name':row['Protein names']})
            continue
        if row['Organism (ID)'] == '11768' or row['Entry'] == 'P11364':
            raise AssertionError('Viral sequence entered host reference panel')
        split = 'train' if row['Organism (ID)'] in ['9606', '10090'] else 'heldout_species'
        row.update(family=family, split=split)
        references.append(row)
        annotations = re.findall(r'DOMAIN ([<>]?\d+)\.\.([<>]?\d+); /note="([^"]+)"', row['Domain [FT]'])
        selected = [(a,b,n) for a,b,n in annotations if 'Ig-like' in n]
        for index, (a,b,note) in enumerate(selected):
            a, b = int(a.strip('<>')), int(b.strip('<>'))
            if family.endswith('_chain'):
                domain_family = family.replace('_chain', '_V' if index == 0 else '_C')
            else:
                domain_family = family
            sequence = row['Sequence'][a-1:b]
            if not 50 <= len(sequence) <= 180:
                continue
            domains.append({'id':f"{row['Entry']}_{a}_{b}", 'accession':row['Entry'],
                            'family':domain_family, 'taxid':row['Organism (ID)'],
                            'split':split, 'sequence':sequence, 'start':a, 'end':b,
                            'boundary_source':'UniProt_DOMAIN',
                            'fuzzy_boundary':bool(re.search(r'[<>]', annotations[index][0]+annotations[index][1])),
                            'fragment':'Fragment' in row['Protein names']})
        if not selected and (family.endswith('_V') or family.startswith('TR')) and 80 <= len(row['Sequence']) <= 200:
            domains.append({'id':row['Entry']+'_whole', 'accession':row['Entry'],
                            'family':family, 'taxid':row['Organism (ID)'], 'split':split,
                            'sequence':row['Sequence'], 'start':1, 'end':len(row['Sequence']),
                            'boundary_source':'whole_named_region', 'fuzzy_boundary':True,
                            'fragment':'Fragment' in row['Protein names']})
    for split in ['train', 'heldout_species']:
        with (cache / f'{split}.faa').open('w') as handle:
            for row in references:
                if row['split'] == split:
                    handle.write(f">{row['Entry']}\n{row['Sequence']}\n")
    with (cache / 'decoys.faa').open('w') as handle:
        for row in decoys:
            handle.write(f">{row['Entry']}\n{row['Sequence']}\n")
    with (cache / 'heldout_domains.faa').open('w') as handle:
        for row in domains:
            if row['split'] == 'heldout_species':
                handle.write(f">{row['id']}\n{row['sequence']}\n")
    # Grouped, nonredundant training sequences only; held-out species never enter these profiles.
    groups = collections.defaultdict(dict)
    for row in domains:
        if row['split'] == 'train' and not row['fragment']:
            groups[row['family']].setdefault(row['sequence'], row['id'])
    alphabet = pyhmmer.easel.Alphabet.amino()
    background = pyhmmer.plan7.Background(alphabet)
    model_info = []
    with (cache / 'receptor_profiles.hmm').open('wb') as output:
        for family, sequences in sorted(groups.items()):
            fasta, aligned = cache / f'{family}.faa', cache / f'{family}.afa'
            fasta.write_text(''.join(f'>{name}\n{sequence}\n' for sequence,name in sequences.items()))
            command = [MUSCLE, '-in', str(fasta), '-out', str(aligned), '-maxiters', '2', '-quiet']
            subprocess.run(command, check=True)
            with pyhmmer.easel.MSAFile(str(aligned), digital=True, alphabet=alphabet, format='afa') as handle:
                msa = handle.read()
            msa.name = family.encode()
            builder = pyhmmer.plan7.Builder(alphabet, seed=20261112)
            hmm, _, _ = builder.build_msa(msa, background)
            hmm.write(output)
            model_info.append({'family':family, 'training_sequences':len(sequences), 'length':hmm.M,
                               'alignment_command':command, 'alignment_sha256':sha256(aligned)})
    with (ROOT / 'manifests/references.tsv').open('w') as handle:
        names = ['Entry','Entry Name','Protein names','Organism (ID)','Organism','family','split','Length','EMBL']
        writer = csv.DictWriter(handle, fieldnames=names, delimiter='\t', extrasaction='ignore')
        writer.writeheader(); writer.writerows(references)
    (cache / 'domains.json').write_text(json.dumps(domains, indent=2)+'\n')
    (cache / 'decoys.json').write_text(json.dumps(decoys, indent=2)+'\n')
    summary = {'completed_utc':utc(), 'source_sha256':sha256(ROOT / 'data/raw/references/uniprot_expanded.tsv'),
               'references':len(references), 'species':len({r['Organism (ID)'] for r in references}),
               'families':dict(collections.Counter(r['family'] for r in references)),
               'splits':dict(collections.Counter(r['split'] for r in references)), 'decoys':len(decoys),
               'domain_splits':dict(collections.Counter(r['split'] for r in domains)), 'profiles':model_info,
               'muscle_path':MUSCLE, 'muscle_sha256':sha256(__import__('pathlib').Path(MUSCLE)),
               'pyhmmer':pyhmmer.__version__,
               'limitations':['No camelid VHH or IgNAR-specific panel yet',
                              'TRD and other rare families have limited training diversity',
                              'Whole named regions used where domain bounds unavailable',
                              'Holdout near-identity filtering is performed during calibration'],
               'outputs':{p.name:sha256(p) for p in cache.iterdir() if p.is_file()}}
    (ROOT / 'results/reference_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['outputs','profiles']},indent=2))


if __name__ == '__main__':
    main()
