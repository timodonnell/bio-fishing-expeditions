#!/usr/bin/env python3
"""Recover deposited neighboring proteins for one retained candidate assembly."""
import argparse
import json
import sys
from pathlib import Path

from Bio import SeqIO, bgzf

from acquire import ROOT, sha256, utc
from acquire_large import atomic_json
from metavr_context import Ffx, Ranges, SOURCES, BASE, locate


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('genome')
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    manifest_path=args.output/'ranges.json'
    manifest=json.loads(manifest_path.read_text()) if manifest_path.exists() else {'sources':{},'ranges':{}}
    SOURCES['protein_index']={'url':BASE+'METAVR.faa.bgz.ffx','bytes':4042159286,
                             'etag':'"f0ee74b6-659fa44e7c540"'}
    index=Ffx(Ranges('protein_index',manifest))
    entries=index.find(args.genome)
    if not entries:
        raise ValueError('No indexed proteins for genome')
    context=ROOT/'results/metavr/context'
    genome=SeqIO.read(context/f'{args.genome}.fna','fasta')
    metadata=json.loads((context/f'{args.genome}.json').read_text())
    archive=ROOT/'data/raw/metavr/METAVR.faa.bgz'
    rows=[]
    with bgzf.BgzfReader(str(archive),'rb') as handle, (args.output/'proteins.faa').open('w') as out:
        for entry in entries:
            n,bases,width=entry['sequence_length'],entry['line_bases'],entry['line_width']
            handle.seek(entry['virtual_offset'])
            sequence=handle.read(n+((n-1)//bases)*(width-bases)).replace(b'\n',b'').replace(b'\r',b'').decode()
            assert len(sequence)==n and set(sequence.upper())<=set('ACDEFGHIKLMNPQRSTVWYBXZJUO*')
            matches=locate(sequence.rstrip('*'),str(genome.seq),int(metadata['genetic_code']))
            if not matches:
                raise ValueError('Neighbor protein does not match genome: '+entry['id'])
            out.write(f'>{entry["header"]}\n{sequence}\n')
            rows.append({'protein':entry['id'],'length':n,'index_entry':entry,'locations':matches})
    (args.output/'locations.json').write_text(json.dumps(rows,indent=2)+'\n')
    manifest.update(status='complete',completed_utc=utc(),command=[sys.executable,*sys.argv],
                    script_sha256=sha256(Path(__file__)),genome=args.genome,proteins=len(rows),
                    protein_index_record_count=index.record_count,
                    archive_sha256=json.loads((ROOT/'manifests/metavr_download.json').read_text())['sha256'],
                    genome_sha256=sha256(context/f'{args.genome}.fna'),
                    protein_fasta_sha256=sha256(args.output/'proteins.faa'),
                    locations_sha256=sha256(args.output/'locations.json'))
    atomic_json(manifest_path,manifest)
    print('Retrieved and mapped',len(rows),'neighbor proteins for',args.genome)


if __name__=='__main__':main()
