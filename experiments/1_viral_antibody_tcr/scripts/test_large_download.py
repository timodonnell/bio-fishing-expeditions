#!/usr/bin/env python3
"""Fault-injection test of range recovery without downloading a real archive."""
import hashlib
import gc
import json
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import requests

import acquire_large
from acquire import ROOT, sha256, utc


def main():
    directory = Path(tempfile.mkdtemp(prefix='range_test_', dir=ROOT/'work'))
    (directory/'manifests').mkdir()
    payload = bytes(range(256))*40
    config = {'url':'https://test.invalid/archive.bgz','bytes':len(payload),'etag':'"test"',
              'last_modified':'test-date','destination':'data/archive.bgz','range_bytes':1024,
              'workers':1,'minimum_free_bytes':0,'maximum_archive_bytes':len(payload)}
    config_path = directory/'config.json'
    config_path.write_text(json.dumps(config))
    calls = []
    failure = [True]
    headers = {'Content-Length':str(len(payload)),'ETag':'"test"','Last-Modified':'test-date'}
    class Response:
        status_code = 200
        content = b''
        def __init__(self):
            self.headers = dict(headers)
        def raise_for_status(self):
            pass
    def get(url, headers, timeout):
        start, end = map(int, headers['Range'].split('=')[1].split('-'))
        calls.append(start//1024)
        if failure[0] and start == 5*1024:
            raise requests.ConnectionError('Injected interruption')
        result = Response()
        result.status_code = 206
        result.headers['Content-Range'] = f'bytes {start}-{end}/{len(payload)}'
        result.content = payload[start:end+1]
        return result
    with patch.object(acquire_large,'ROOT',directory), \
         patch.object(sys,'argv',[str(ROOT/'scripts/acquire_large.py'),str(config_path)]), \
         patch.object(acquire_large.requests,'head',return_value=Response()), \
         patch.object(acquire_large.requests,'get',side_effect=get), \
         patch.object(acquire_large.time,'sleep',return_value=None):
        try:
            acquire_large.main()
            raise AssertionError('Expected failure')
        except requests.ConnectionError:
            pass
        # A real failed CLI process releases its lock on exit. Collect the
        # executor/exception cycle when exercising multiple invocations here.
        gc.collect()
        state = json.loads((directory/'manifests/metavr_download.json').read_text())
        assert state['status']=='interrupted' and state['completed_ranges']
        damaged = min(map(int,state['completed_ranges']))
        partial = directory/'data/archive.bgz.part'
        with partial.open('r+b') as handle:
            handle.seek(damaged*1024)
            handle.write(b'BAD')
        failure[0] = False
        calls.clear()
        acquire_large.main()
        assert damaged in calls and 5 in calls
        assert (directory/'data/archive.bgz').read_bytes()==payload
        complete = json.loads((directory/'manifests/metavr_download.json').read_text())
        assert complete['sha256']==hashlib.sha256(payload).hexdigest()
        calls.clear()
        acquire_large.main()
        assert not calls
        # A changed source is rejected before an incomplete archive is resumed.
        (directory/'data/archive.bgz').unlink()
        changed = Response()
        changed.headers['ETag']='"changed"'
        with patch.object(acquire_large.requests,'head',return_value=changed):
            try:
                acquire_large.main()
                raise AssertionError('Changed source accepted')
            except ValueError as error:
                assert 'Source identity changed' in str(error)
    report = {'status':'passed','completed_utc':utc(),
              'checks':['interrupted range retrieval','corrupt completed range redownloaded',
                        'resumed bytes identical to source','completed download reuse verified',
                        'changed source ETag rejected'],
              'script_sha256':sha256(ROOT/'scripts/acquire_large.py'),
              'test_script_sha256':sha256(Path(__file__)),
              'method':'Mock HTTP responses, real file writes and checksums; no external test server'}
    (ROOT/'results/scale_download_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
