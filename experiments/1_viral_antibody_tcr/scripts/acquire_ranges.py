#!/usr/bin/env python3
"""Bounded parallel HTTP ranges for large public files; verifies ranges and ETag."""
import argparse
import concurrent.futures
import fcntl
import json
import shutil
import time
from pathlib import Path

import requests

from acquire import ROOT, acquire, sha256, utc, validate_content


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('spec',type=Path)
    args=parser.parse_args()
    manifest_path=ROOT/'manifests/downloads.json'
    for spec in json.loads(args.spec.read_text()):
        manifest=json.loads(manifest_path.read_text())
        dest=ROOT/spec['path']
        if manifest.get(spec['id'],{}).get('status')=='complete' or not spec['path'].endswith('.gz'):
            entry=acquire(spec,manifest)
        else:
            started=utc()
            head=requests.head(spec['url'],timeout=(30,60));head.raise_for_status()
            size=int(head.headers['Content-Length']);etag=head.headers.get('ETag')
            assert size <= spec.get('max_bytes',2*2**30)
            assert shutil.disk_usage(ROOT).free-size > 50*2**30
            dest.parent.mkdir(parents=True,exist_ok=True)
            partial=dest.with_name(dest.name+'.part')
            prefix=partial.stat().st_size if partial.exists() else 0
            # A partial prefix from the earlier single-stream transfer is retained only when
            # the preflight HEAD ETag matches today's HEAD. Otherwise download from byte zero.
            old_heads=json.loads((ROOT/'manifests/remote_sizes.json').read_text())['requests']
            prior=next((r for r in old_heads if r['url']==spec['url']),{})
            if not etag or prior.get('headers',{}).get('ETag')!=etag:prefix=0
            prefix_hash=sha256(partial) if prefix else None
            chunks=ROOT/'data/cache/download_ranges'/spec['id'];chunks.mkdir(parents=True,exist_ok=True)
            ranges=[(i,min(i+4*2**20,size)-1) for i in range(prefix,size,4*2**20)]
            def fetch(bounds):
                begin,end=bounds;out=chunks/f'{begin}-{end}'
                for attempt in range(3):
                    try:
                        headers={'Range':f'bytes={begin}-{end}','If-Match':etag,'Accept-Encoding':'identity'}
                        response=requests.get(spec['url'],headers=headers,timeout=(30,120))
                        assert response.status_code==206,(response.status_code,response.text[:100])
                        assert response.headers['Content-Range']==f'bytes {begin}-{end}/{size}'
                        assert response.headers.get('ETag')==etag
                        assert len(response.content)==end-begin+1
                        out.write_bytes(response.content)
                        return out
                    except (requests.RequestException,AssertionError):
                        if attempt==2:raise
                        time.sleep(2**attempt)
            with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
                paths=list(pool.map(fetch,ranges))
            with partial.open('r+b' if prefix else 'wb') as output:
                if prefix:output.seek(prefix);output.truncate()
                for path in paths:
                    output.write(path.read_bytes());path.unlink()
            assert partial.stat().st_size==size
            validate_content(partial,dest)
            checksum=sha256(partial)
            if spec.get('expected_sha256'):assert checksum==spec['expected_sha256']
            partial.replace(dest)
            entry=dict(spec,status='complete',started_utc=started,finished_utc=utc(),bytes=size,
                       sha256=checksum,headers=dict(head.headers),http_status=206,final_url=head.url,
                       acquisition='8 parallel HTTP range workers; 4 MiB chunks; If-Match and Content-Range verified',
                       resumed_prefix_bytes=prefix,resumed_prefix_sha256=prefix_hash)
        with manifest_path.with_suffix('.lock').open('w') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            manifest=json.loads(manifest_path.read_text());manifest[spec['id']]=entry
            tmp=manifest_path.with_suffix('.tmp');tmp.write_text(json.dumps(manifest,indent=2)+'\n');tmp.replace(manifest_path)
        print(spec['id'],entry['status'],entry.get('bytes'),flush=True)


if __name__=='__main__':main()
