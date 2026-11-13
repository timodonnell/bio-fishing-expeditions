#!/usr/bin/env python3
"""Download public inputs atomically and retain a checksum/provenance manifest."""
import argparse
import fcntl
import hashlib
import json
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def validate_content(path, destination):
    with path.open('rb') as handle:
        prefix = handle.read(200)
    if destination.suffix == '.pdf' and not prefix.startswith(b'%PDF-'):
        raise RuntimeError('Expected PDF; server returned a different content format')
    if destination.suffix == '.gz' and not prefix.startswith(b'\x1f\x8b'):
        raise RuntimeError('Expected gzip data; server returned a different content format')
    if destination.suffix == '.gb' and not prefix.startswith(b'LOCUS'):
        raise RuntimeError('Expected GenBank record')
    if destination.suffix == '.json':
        json.loads(path.read_text())


def acquire(spec, manifest):
    destination = ROOT / spec['path']
    destination.parent.mkdir(parents=True, exist_ok=True)
    old = manifest.get(spec['id'])
    if (old and old.get('status') == 'complete' and old['url'] == spec['url']
            and destination.exists() and sha256(destination) == old['sha256']
            and (not spec.get('expected_sha256') or old['sha256'] == spec['expected_sha256'])):
        try:
            validate_content(destination, destination)
            print('cached', spec['id'], old['bytes'], flush=True)
            return old
        except (RuntimeError, ValueError):
            pass
    entry = dict(spec, started_utc=utc(), status='running')
    partial = destination.with_name(destination.name + '.part')
    for attempt in range(3):
        try:
            with requests.get(spec['url'], stream=True, timeout=(30, 120),
                              headers={'User-Agent': 'bio-fishing-expeditions/1.0'}) as response:
                entry.update(http_status=response.status_code, final_url=response.url,
                             headers={key: response.headers.get(key) for key in
                                      ['Content-Length', 'Content-Type', 'Last-Modified', 'ETag',
                                       'X-UniProt-Release', 'X-UniProt-Release-Date']})
                response.raise_for_status()
                total = 0
                digest = hashlib.sha256()
                with partial.open('wb') as handle:
                    for block in response.iter_content(1024 * 1024):
                        if shutil.disk_usage(ROOT).free < 50 * 2**30:
                            raise RuntimeError('Less than 50 GiB filesystem headroom')
                        total += len(block)
                        if total > spec.get('max_bytes', 2 * 2**30):
                            raise RuntimeError('Download exceeds declared maximum size')
                        handle.write(block)
                        digest.update(block)
                validate_content(partial, destination)
                if spec.get('expected_sha256') and digest.hexdigest() != spec['expected_sha256']:
                    raise RuntimeError('Checksum does not match pinned input')
                partial.replace(destination)
                entry.update(status='complete', bytes=total, sha256=digest.hexdigest(),
                             finished_utc=utc())
                print('downloaded', spec['id'], total, flush=True)
                return entry
        except (requests.RequestException, RuntimeError, ValueError) as error:
            entry.update(status='failed', error=str(error), finished_utc=utc())
            if attempt < 2:
                time.sleep(2**attempt)
    partial.unlink(missing_ok=True)
    print('failed', spec['id'], entry['error'], flush=True)
    return entry


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('spec', type=Path)
    args = parser.parse_args()
    manifest_path = ROOT / 'manifests/downloads.json'
    manifest_path.parent.mkdir(exist_ok=True)
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    failed = False
    for spec in json.loads(args.spec.read_text()):
        entry = acquire(spec, manifest)
        failed |= entry['status'] != 'complete'
        # Merge under a lock: independent downloads must not overwrite each other's provenance.
        with manifest_path.with_suffix('.lock').open('w') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
            manifest[spec['id']] = entry
            temporary = manifest_path.with_suffix('.tmp')
            temporary.write_text(json.dumps(manifest, indent=2) + '\n')
            temporary.replace(manifest_path)
    raise SystemExit(int(failed))


if __name__ == '__main__':
    main()
