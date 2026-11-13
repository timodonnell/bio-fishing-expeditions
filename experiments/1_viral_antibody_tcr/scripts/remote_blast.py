#!/usr/bin/env python3
"""Submit or poll a candidate-specific public NCBI nucleotide comparison."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import re
import sys
from pathlib import Path

import requests

from acquire import sha256, utc
from acquire_large import atomic_json

URL='https://blast.ncbi.nlm.nih.gov/Blast.cgi'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['submit','poll'])
    parser.add_argument('output',type=Path)
    parser.add_argument('--query',type=Path)
    parser.add_argument('--program',choices=['blastn','blastp'],default='blastn')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    path=args.output/'remote_blast.json'
    if args.action=='submit':
        if path.exists():
            raise ValueError('Existing request; poll its RID instead of submitting again')
        if not args.query:
            parser.error('--query required for submit')
        params={'CMD':'Put','PROGRAM':'blastn','MEGABLAST':'on','DATABASE':'core_nt',
                'EXPECT':'1e-20','HITLIST_SIZE':'20','FILTER':'mL','WORD_SIZE':'28',
                'tool':'bio-fishing-expeditions'}
        if args.program=='blastp':
            params.update(PROGRAM='blastp',DATABASE='nr',EXPECT='1e-5',WORD_SIZE='3')
            del params['MEGABLAST']
        response=requests.post(URL,data=dict(params,QUERY=args.query.read_text()),timeout=120)
        response.raise_for_status()
        (args.output/'submission.html').write_bytes(response.content)
        rid=re.search(r'RID\s*=\s*(\S+)',response.text)
        rtoe=re.search(r'RTOE\s*=\s*(\S+)',response.text)
        record={'request_parameters':params,'query_file':str(args.query),'query_sha256':sha256(args.query),
                'response_sha256':hashlib.sha256(response.content).hexdigest(),
                'RID':rid[1] if rid else None,'RTOE':rtoe[1] if rtoe else None,
                'status':'submitted' if rid else 'submission_failed','submitted_utc':utc(),
                'command':[sys.executable,*sys.argv],'script_sha256':sha256(Path(__file__))}
        atomic_json(path,record)
    else:
        record=json.loads(path.read_text())
        if record['status']=='complete':
            print('Already complete',record['RID'])
            return
        previous=record.get('last_polled_utc',record['submitted_utc'])
        elapsed=(datetime.now(timezone.utc)-datetime.fromisoformat(previous)).total_seconds()
        if elapsed<60:
            print('Respecting NCBI polling interval; retry after',round(60-elapsed),'seconds')
            return
        response=requests.get(URL,params={'CMD':'Get','RID':record['RID'],'FORMAT_TYPE':'JSON2_S'},timeout=120)
        response.raise_for_status()
        record['last_polled_utc']=utc()
        record['poll_command']=[sys.executable,*sys.argv]
        record['poll_script_sha256']=sha256(Path(__file__))
        try:
            data=response.json()
        except ValueError:
            data={}
        if 'BlastOutput2' in data:
            destination=args.output/'results.json'
            destination.write_bytes(response.content)
            record.update(status='complete',result_sha256=sha256(destination),finished_utc=utc())
        else:
            (args.output/'last_poll.html').write_bytes(response.content)
            status=re.search(r'Status=(\w+)',response.text)
            record['status']=status[1].lower() if status else 'response_requires_inspection'
        atomic_json(path,record)
    print(record['RID'],record['status'])


if __name__=='__main__':main()
