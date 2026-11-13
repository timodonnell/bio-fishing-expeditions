#!/usr/bin/env python3
"""Resume a pinned public archive in place, checking every completed HTTP range."""
import argparse
import concurrent.futures
import fcntl
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import requests

from acquire import ROOT, sha256, utc


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    os.replace(temporary, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    dest = ROOT / config['destination']
    dest.parent.mkdir(parents=True, exist_ok=True)
    state_path = ROOT / 'manifests/metavr_download.json'
    lock = state_path.with_suffix('.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    config_hash = sha256(args.config)
    state = json.loads(state_path.read_text()) if state_path.exists() else {
        'config_sha256': config_hash, 'config': config, 'started_utc': utc(),
        'command': [sys.executable, *sys.argv], 'script_sha256': sha256(Path(__file__)),
        'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'requests_version': requests.__version__, 'completed_ranges': {}, 'attempts': []}
    if state['config_sha256'] != config_hash:
        raise ValueError('Configuration changed; use a separately named acquisition')
    if dest.exists():
        if state.get('status') != 'complete' or sha256(dest) != state['sha256']:
            raise ValueError('Existing final archive is not validated')
        print('Already complete and checksum verified', flush=True)
        return
    response = requests.head(config['url'], timeout=60)
    response.raise_for_status()
    if (int(response.headers['Content-Length']) != config['bytes'] or
            response.headers.get('ETag') != config['etag'] or
            response.headers.get('Last-Modified') != config['last_modified']):
        raise ValueError('Source identity changed')
    state['response_headers'] = dict(response.headers)
    state['attempts'].append({'started_utc': utc(), 'script_sha256': sha256(Path(__file__))})
    partial = dest.with_suffix(dest.suffix + '.part')
    size, chunk = config['bytes'], config['range_bytes']
    if size > config['maximum_archive_bytes']:
        raise ValueError('Archive exceeds configured budget')
    allocated = partial.stat().st_blocks * 512 if partial.exists() else 0
    if shutil.disk_usage(dest.parent).free - max(0, size-allocated) < config['minimum_free_bytes']:
        raise RuntimeError('Insufficient space to finish while preserving free-space reserve')
    fd = os.open(partial, os.O_RDWR | os.O_CREAT, 0o644)
    os.ftruncate(fd, size)
    try:
        # A crash can lose writes after a checkpoint; hashes, not file size, establish completion.
        invalid = []
        for key, row in state['completed_ranges'].items():
            data = os.pread(fd, row['bytes'], int(key)*chunk)
            if len(data) != row['bytes'] or hashlib.sha256(data).hexdigest() != row['sha256']:
                invalid.append(key)
        for key in invalid:
            del state['completed_ranges'][key]
        state['status'] = 'running'
        atomic_json(state_path, state)

        def fetch(index):
            start, end = index*chunk, min(size, (index+1)*chunk)-1
            for attempt in range(5):
                try:
                    if shutil.disk_usage(dest.parent).free < config['minimum_free_bytes'] + chunk*config['workers']:
                        raise RuntimeError('Free-space reserve reached')
                    result = requests.get(config['url'], headers={
                        'Range': f'bytes={start}-{end}', 'If-Match': config['etag'],
                        'Accept-Encoding': 'identity'}, timeout=(30, 120))
                    result.raise_for_status()
                    if (result.status_code != 206 or
                            result.headers.get('Content-Range') != f'bytes {start}-{end}/{size}' or
                            result.headers.get('ETag') != config['etag']):
                        raise ValueError('Untrusted range response')
                    data = result.content
                    if len(data) != end-start+1:
                        raise ValueError('Truncated range')
                    digest = hashlib.sha256(data).hexdigest()
                    written = 0
                    while written < len(data):
                        n = os.pwrite(fd, memoryview(data)[written:], start+written)
                        if n <= 0:
                            raise OSError('Short write')
                        written += n
                    return index, {'bytes': len(data), 'sha256': digest}
                except (requests.RequestException, ValueError):
                    if attempt == 4:
                        raise
                    time.sleep(2**attempt)

        pending_indices = [i for i in range((size+chunk-1)//chunk)
                           if str(i) not in state['completed_ranges']]
        started = time.monotonic()
        initial = sum(row['bytes'] for row in state['completed_ranges'].values())
        with concurrent.futures.ThreadPoolExecutor(max_workers=config['workers']) as pool:
            futures = {pool.submit(fetch, i): i for i in pending_indices}
            for future in concurrent.futures.as_completed(futures):
                index, row = future.result()
                state['completed_ranges'][str(index)] = row
                state['downloaded_bytes'] = sum(r['bytes'] for r in state['completed_ranges'].values())
                state['updated_utc'] = utc()
                atomic_json(state_path, state)
                if len(state['completed_ranges']) % 16 == 0:
                    elapsed = time.monotonic()-started
                    print(f"{state['downloaded_bytes']/size:.1%} "
                          f"{state['downloaded_bytes']/2**30:.2f} GiB "
                          f"{(state['downloaded_bytes']-initial)/elapsed/2**20:.1f} MiB/s", flush=True)
        os.fsync(fd)
        state['sha256'] = sha256(partial)
        if state['downloaded_bytes'] != size:
            raise ValueError('Incomplete range coverage')
        os.replace(partial, dest)
        state.update(status='complete', finished_utc=utc(), bytes=size,
                     path=str(dest.relative_to(ROOT)))
        atomic_json(state_path, state)
        print('Complete', size, state['sha256'], flush=True)
    except BaseException as error:
        state.update(status='interrupted', error=repr(error), updated_utc=utc())
        atomic_json(state_path, state)
        raise
    finally:
        os.close(fd)


if __name__ == '__main__':
    main()
