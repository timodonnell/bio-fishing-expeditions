#!/usr/bin/env python3
"""Checkpointed HMM screen of every protein in a BGZF FASTA; no ancestry calls."""
import argparse
import csv
import concurrent.futures
import fcntl
import gzip
import hashlib
import io
import json
import multiprocessing
import os
import re
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path

import Bio
from Bio import bgzf
import pyhmmer

from acquire import ROOT, sha256, utc
from acquire_large import atomic_json

BASE_Z = 1_000_000
FIELDS = ['query', 'target', 'sequence_E', 'sequence_score', 'domain_iE',
          'domain_score', 'bias', 'ali_from', 'ali_to', 'hmm_from', 'hmm_to',
          'target_length', 'model_length', 'hmm_alignment', 'target_alignment']


def read_shard(reader, minimum_bytes):
    """Read to a FASTA boundary and return bytes with a reproducible BGZF offset.

    Only the final 64 KiB is reread to position the checkpoint exactly before
    the next header. This also handles sequences larger than one BGZF block.
    """
    parts, size = [], 0
    while True:
        position = reader.tell()
        block = reader.read(65536)
        if not block:
            break
        split = block.rfind(b'\n>') + 1
        if size + len(block) >= minimum_bytes and split:
            parts.append(block[:split])
            reader.seek(position)
            reader.read(split)
            break
        parts.append(block)
        size += len(block)
    data = b''.join(parts)
    if data and not data.startswith(b'>'):
        raise ValueError('Shard does not start at a FASTA header')
    return data


def search_block(hmms, sequences, cpus):
    rows = []
    text = lambda value: value.decode() if isinstance(value, bytes) else value
    for hits in pyhmmer.hmmer.hmmsearch(hmms, sequences, cpus=cpus,
                                      E=100, domE=100, Z=BASE_Z, seed=42):
        for hit in hits:
            for domain in hit.domains:
                if domain.i_evalue > 0.01:
                    continue
                alignment = domain.alignment
                rows.append(dict(zip(FIELDS, [text(hits.query.name), text(hit.name),
                    hit.evalue, hit.score, domain.i_evalue, domain.score, domain.bias,
                    alignment.target_from, alignment.target_to, alignment.hmm_from,
                    alignment.hmm_to, hit.length, hits.query.M,
                    alignment.hmm_sequence, alignment.target_sequence])))
    return rows


def write_tsv(path, rows, fields=FIELDS):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def process_shard(data, number, shard_start, shard_end, output, coverage, profiles, cpus):
    started = time.monotonic()
    with pyhmmer.plan7.HMMFile(str(profiles)) as handle:
        hmms = list(handle)
    raw_hash = hashlib.sha256(data).hexdigest()
    ids = re.findall(rb'^>(\S+)', data, re.MULTILINE)
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate protein identifiers within shard')
    genomes = {name.split(b'|', 1)[0] for name in ids}
    with pyhmmer.easel.SequenceFile(io.BytesIO(data), digital=True,
            alphabet=pyhmmer.easel.Alphabet.amino()) as handle:
        block = handle.read_block()
    if len(block) != len(ids):
        raise ValueError('FASTA header/parser count mismatch')
    rows = search_block(hmms, block, cpus)
    keep = {row['target'] for row in rows}
    prefix = output/'shards'/f'{number:06d}'
    write_tsv(prefix.with_suffix('.tsv'), rows)
    with prefix.with_suffix('.faa').open('wb') as handle:
        for sequence in block:
            name = sequence.name.decode() if isinstance(sequence.name, bytes) else sequence.name
            if name in keep:
                sequence.textize().write(handle)
    genome_path = coverage/f'{number:06d}.txt.gz'
    with gzip.open(genome_path, 'wb') as handle:
        handle.write(b'\n'.join(sorted(genomes))+b'\n')
    return {'shard': number, 'completed_utc': utc(), 'records': len(block),
            'residues': sum(len(s) for s in block), 'uncompressed_bytes': len(data),
            'fasta_sha256': raw_hash, 'start_virtual_offset': shard_start,
            'end_virtual_offset': shard_end, 'source_genome_ids_in_shard': len(genomes),
            'headers_without_pipe': sum(b'|' not in name for name in ids),
            'domain_rows': len(rows), 'candidate_proteins': len(keep),
            'worker_seconds': time.monotonic()-started,
            'worker_maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'outputs': {str(p): sha256(p) for p in [prefix.with_suffix('.tsv'),
                       prefix.with_suffix('.faa'), genome_path]}}


def finalize(output, coverage, checkpoints, effective_count=None):
    count = sum(row['records'] for row in checkpoints)
    z = effective_count or count
    if z < BASE_Z:
        raise ValueError('Final search space below retention reference; use a >=1M calibration override')
    genomes = set()
    rows = []
    candidates = set()
    for checkpoint in checkpoints:
        number = checkpoint['shard']
        with gzip.open(coverage / f'{number:06d}.txt.gz', 'rt') as handle:
            genomes.update(line.rstrip('\n') for line in handle)
        with (output/'shards'/f'{number:06d}.tsv').open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                row['sequence_E'] = float(row['sequence_E']) * z / BASE_Z
                row['domain_iE'] = float(row['domain_iE']) * z / BASE_Z
                if row['domain_iE'] <= 0.01:
                    span = int(row['ali_to']) - int(row['ali_from']) + 1
                    row['tier'] = 'main' if row['domain_iE'] <= 0.001 and span >= 50 else 'rescue'
                    row['shard'] = number
                    rows.append(row)
                    candidates.add(row['target'])
    write_tsv(output/'hits.tsv', rows, FIELDS+['tier', 'shard'])
    seen = {}
    from Bio import SeqIO
    with (output/'candidates.faa').open('w') as out:
        for checkpoint in checkpoints:
            path = output/'shards'/f"{checkpoint['shard']:06d}.faa"
            for record in SeqIO.parse(path, 'fasta'):
                if record.id not in candidates:
                    continue
                sequence = str(record.seq)
                if record.id in seen and seen[record.id] != sequence:
                    raise ValueError('Candidate identifier has conflicting sequences')
                if record.id not in seen:
                    out.write(f'>{record.id}\n{sequence}\n')
                seen[record.id] = sequence
    if set(seen) != candidates:
        raise ValueError('Missing retained candidate sequence')
    return {'protein_records': count, 'residues': sum(row['residues'] for row in checkpoints),
            'uncompressed_fasta_bytes': sum(row['uncompressed_bytes'] for row in checkpoints),
            'distinct_source_genome_identifiers': len(genomes), 'shards': len(checkpoints),
            'headers_without_pipe': sum(row['headers_without_pipe'] for row in checkpoints),
            'effective_target_count': z, 'raw_domain_rows': sum(row['domain_rows'] for row in checkpoints),
            'retained_domain_rows': len(rows), 'retained_candidate_proteins': len(candidates),
            'main_candidate_proteins': len({r['target'] for r in rows if r['tier']=='main'}),
            'outputs': {name: sha256(output/name) for name in ['hits.tsv', 'candidates.faa']},
            'genome_count_scope': 'Distinct header token before first |; proteins only, not genomes searched by DNA',
            'overlap_with_prior_corpora': 'Not deduplicated across catalogs'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('target', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--profiles', type=Path, default=ROOT/'references/receptor_profiles.hmm')
    parser.add_argument('--cpus', type=int, default=14)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--shard-mib', type=int, default=256)
    parser.add_argument('--effective-target-count', type=int)
    parser.add_argument('--stop-after-shards', type=int, help='Intentional restart test; not complete coverage')
    parser.add_argument('--minimum-free-gib', type=int, default=50)
    args = parser.parse_args()
    if not 1 <= args.cpus * args.workers <= 16 or min(args.cpus, args.workers, args.shard_mib) < 1:
        parser.error('Use positive worker/shard sizes and at most 16 total search threads')
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'shards').mkdir(exist_ok=True)
    coverage_key = args.output.name + '-' + hashlib.sha256(str(args.output.resolve()).encode()).hexdigest()[:12]
    coverage = ROOT/'data/cache/metavr/coverage'/coverage_key
    coverage.mkdir(parents=True, exist_ok=True)
    lock = (args.output/'scan.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    manifest_path = args.output/'scan.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {
        'started_utc': utc(), 'status': 'initializing', 'attempts': []}
    identity = {'target_sha256': sha256(args.target), 'target_bytes': args.target.stat().st_size,
                'profiles_sha256': sha256(args.profiles), 'script_sha256': sha256(Path(__file__)),
                'pyhmmer': pyhmmer.__version__, 'biopython': Bio.__version__,
                'parameters': {'cpus_per_worker': args.cpus, 'workers': args.workers,
                    'shard_bytes': args.shard_mib*2**20,
                    'reference_Z': BASE_Z, 'E': 100, 'domE': 100, 'seed': 42,
                    'retain_reference_independent_domain_E': 0.01,
                    'final_Z_override': args.effective_target_count}}
    if 'identity' in manifest and manifest['identity'] != identity:
        raise ValueError('Input, source, environment, or parameters changed since checkpoint')
    manifest.update(identity=identity, target=str(args.target), profiles=str(args.profiles),
                    git_commit=subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip())
    manifest['attempts'].append({'started_utc': utc(), 'command': [sys.executable, *sys.argv]})
    checkpoints = []
    for path in sorted((args.output/'shards').glob('*.json')):
        row = json.loads(path.read_text())
        if row['shard'] != len(checkpoints) or row['start_virtual_offset'] != (
                checkpoints[-1]['end_virtual_offset'] if checkpoints else 0):
            raise ValueError('Noncontiguous checkpoint coverage')
        for filename, expected in row['outputs'].items():
            if sha256(Path(filename)) != expected:
                raise ValueError(f'Checkpoint output changed: {filename}')
        checkpoints.append(row)
    manifest.update(status='running', completed_shards=len(checkpoints))
    atomic_json(manifest_path, manifest)
    started = time.monotonic()
    try:
        with bgzf.BgzfReader(str(args.target), 'rb') as reader, \
             concurrent.futures.ProcessPoolExecutor(max_workers=args.workers,
                 mp_context=multiprocessing.get_context('spawn')) as pool:
            if checkpoints:
                reader.seek(checkpoints[-1]['end_virtual_offset'])
            pending = []
            finished_reading = False
            while pending or not finished_reading:
                while len(pending) < args.workers and not finished_reading:
                    number = len(checkpoints) + len(pending)
                    if args.stop_after_shards and number >= args.stop_after_shards:
                        finished_reading = True
                        break
                    shard_start = reader.tell()
                    data = read_shard(reader, args.shard_mib*2**20)
                    if not data:
                        finished_reading = True
                        break
                    if shutil.disk_usage(ROOT).free < args.minimum_free_gib*2**30:
                        raise RuntimeError('Free-space reserve reached')
                    pending.append(pool.submit(process_shard, data, number, shard_start,
                        reader.tell(), args.output, coverage, args.profiles, args.cpus))
                    del data
                if not pending:
                    break
                checkpoint = pending.pop(0).result()
                number = checkpoint['shard']
                prefix = args.output/'shards'/f'{number:06d}'
                atomic_json(prefix.with_suffix('.json'), checkpoint)
                checkpoints.append(checkpoint)
                manifest.update(completed_shards=len(checkpoints), updated_utc=utc(),
                                searched_records=sum(c['records'] for c in checkpoints),
                                compressed_position=checkpoint['end_virtual_offset'] >> 16)
                atomic_json(manifest_path, manifest)
                print(f"shard {number}: {checkpoint['records']:,} proteins, "
                      f"{checkpoint['candidate_proteins']} candidates; "
                      f"total {manifest['searched_records']:,}; "
                      f"{100*manifest['compressed_position']/identity['target_bytes']:.1f}% of archive; "
                      f"worker {checkpoint['worker_seconds']:.1f}s", flush=True)
                if args.stop_after_shards and len(checkpoints) >= args.stop_after_shards:
                    manifest.update(status='paused_for_restart_test', updated_utc=utc())
                    atomic_json(manifest_path, manifest)
                    return
        summary = finalize(args.output, coverage, checkpoints, args.effective_target_count)
        manifest.update(status='complete', finished_utc=utc(), summary=summary,
                        elapsed_this_attempt_seconds=time.monotonic()-started,
                        maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        atomic_json(manifest_path, manifest)
        print(json.dumps(summary, indent=2), flush=True)
    except BaseException as error:
        manifest.update(status='interrupted', error=repr(error), updated_utc=utc())
        atomic_json(manifest_path, manifest)
        raise


if __name__ == '__main__':
    main()
