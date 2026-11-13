#!/usr/bin/env python3
"""Retrieve only candidate genomes through MetaVR's FASTARS2 index and HTTP ranges.

Index layout reference: Matt-BF/fastars src/ffx/{format,block}.rs, Git blobs
3dcc2d09d9e8072bde70ad463f16b1018605149f and 78f4e10e199855d422bda67590943ab8aec1a9ad.
This independent reader accepts format v2 only and validates its directory.
"""
import argparse
import bisect
import hashlib
import io
import json
import struct
import sys
import time
from pathlib import Path

import requests
import zstandard
from Bio import SeqIO, bgzf
from Bio.Seq import Seq

from acquire import ROOT, sha256, utc
from acquire_large import atomic_json

BASE = 'https://portal.nersc.gov/cfs/m342/METAVR/'
SOURCES = {
    'index': {'url':BASE+'METAVR.fna.bgz.ffx', 'bytes':378283829,
              'etag':'"168c2735-659fa441225c0"'},
    'fasta': {'url':BASE+'METAVR.fna.bgz', 'bytes':71552362108,
              'etag':'"10a8da667c-659fa54a41740"'}}


class Ranges:
    def __init__(self, name, manifest):
        self.name, self.spec, self.manifest = name, SOURCES[name], manifest
        self.cache = ROOT/'data/cache/metavr/remote_ranges'/name
        self.cache.mkdir(parents=True, exist_ok=True)
        response = requests.head(self.spec['url'], timeout=60)
        response.raise_for_status()
        if (response.headers.get('ETag') != self.spec['etag'] or
                int(response.headers['Content-Length']) != self.spec['bytes']):
            raise ValueError('Remote source changed: '+name)
        manifest['sources'][name] = dict(self.spec, response_headers=dict(response.headers))

    def get(self, start, size):
        if size <= 0 or start < 0 or start+size > self.spec['bytes']:
            raise ValueError('Out-of-bounds range')
        path = self.cache/f'{start}_{size}.bin'
        key = self.name+':'+path.name
        previous = self.manifest['ranges'].get(key)
        if previous and path.exists() and sha256(path)==previous['sha256']:
            return path.read_bytes()
        for attempt in range(4):
            try:
                response = requests.get(self.spec['url'], headers={
                    'Range':f'bytes={start}-{start+size-1}', 'If-Match':self.spec['etag'],
                    'Accept-Encoding':'identity'}, timeout=(30,90))
                response.raise_for_status()
                if (response.status_code != 206 or response.headers.get('ETag') != self.spec['etag'] or
                        response.headers.get('Content-Range') != f'bytes {start}-{start+size-1}/{self.spec["bytes"]}' or
                        len(response.content) != size):
                    raise ValueError('Invalid range response')
                data = response.content
                path.write_bytes(data)
                self.manifest['ranges'][key] = {'offset':start,'bytes':size,
                    'sha256':hashlib.sha256(data).hexdigest(),'path':str(path),'retrieved_utc':utc()}
                return data
            except requests.RequestException:
                if attempt==3:
                    raise
                time.sleep(2**attempt)


class RemoteFile:
    def __init__(self, ranges):
        self.ranges, self.position = ranges, 0
        self.page_start, self.page = -1, b''
    def tell(self):
        return self.position
    def seek(self, position, whence=0):
        if whence != 0 or position < 0:
            raise ValueError('Only absolute seeks are supported')
        self.position = position
        return position
    def read(self, size=-1):
        if size < 0:
            raise ValueError('Unbounded reads forbidden')
        pieces = []
        while size and self.position < self.ranges.spec['bytes']:
            start = self.position // 2**20 * 2**20
            if start != self.page_start:
                self.page = self.ranges.get(start, min(2**20, self.ranges.spec['bytes']-start))
                self.page_start = start
            offset = self.position-start
            piece = self.page[offset:offset+size]
            pieces.append(piece)
            self.position += len(piece)
            size -= len(piece)
        return b''.join(pieces)
    def close(self):
        pass


class Ffx:
    def __init__(self, ranges):
        self.ranges = ranges
        header = ranges.get(0,80)
        magic, version, length, count, blocks, offset, size, directory, directory_size, target, reserved, flags = struct.unpack('<8sII6QIIQ', header)
        if (magic!=b'FASTARS2' or version!=2 or length!=80 or offset!=80 or
                offset+size!=directory or directory+directory_size!=ranges.spec['bytes'] or reserved or flags):
            raise ValueError('Unsupported/corrupt FFX header')
        data = ranges.get(directory,directory_size)
        self.entries, self.last_ids = [], []
        cursor, next_offset, next_record = 0,80,0
        for _ in range(blocks):
            block_offset, stored, raw, first, number, id_length, flags = struct.unpack_from('<7Q',data,cursor)
            cursor += 56
            last_id = data[cursor:cursor+id_length].decode()
            cursor += id_length
            if (block_offset!=next_offset or first!=next_record or flags & ~15 or
                    not min(stored,raw,number) or (self.last_ids and last_id<self.last_ids[-1])):
                raise ValueError('Invalid FFX directory')
            self.entries.append((block_offset,stored,raw,first,number,flags,last_id))
            self.last_ids.append(last_id)
            next_offset += stored
            next_record += number
        if cursor!=len(data) or next_offset!=directory or next_record!=count:
            raise ValueError('FFX directory totals disagree')
        self.record_count = count

    def block(self, entry):
        offset, stored, raw_size, first, count, flags, last_id = entry
        payload = self.ranges.get(offset,stored)
        data = zstandard.ZstdDecompressor().decompress(payload,max_output_size=raw_size) if flags & 1 else payload
        if len(data)!=raw_size:
            raise ValueError('Wrong uncompressed index length')
        cursor, previous_id, previous_offset = 0,b'',0
        def take(n):
            nonlocal cursor
            value = data[cursor:cursor+n]
            cursor += n
            if len(value)!=n:
                raise ValueError('Truncated index field')
            return value
        def integer():
            value = 0
            for shift in range(0,64,7):
                byte = take(1)[0]
                if shift==63 and byte>1:
                    raise ValueError('Overflowing index integer')
                value |= (byte & 127) << shift
                if not byte & 128:
                    return value
            raise ValueError('Invalid index integer')
        result=[]
        for i in range(count):
            if i and flags & 2:
                prefix, suffix = integer(), integer()
                if prefix>len(previous_id):
                    raise ValueError('Invalid front-coded identifier')
                identifier = previous_id[:prefix]+take(suffix)
            else:
                identifier = take(integer())
            header = identifier + (take(integer()) if flags & 8 else b'')
            encoded = integer()
            virtual = previous_offset + ((encoded >> 1) ^ -(encoded & 1)) if i and flags & 4 else encoded
            sequence_length, bases, width = integer(), integer(), integer()
            if not bases or width<bases or virtual<0 or identifier<previous_id:
                raise ValueError('Invalid index record')
            result.append({'id':identifier.decode(),'header':header.decode(),'virtual_offset':virtual,
                           'sequence_length':sequence_length,'line_bases':bases,'line_width':width})
            previous_id, previous_offset = identifier,virtual
        if cursor!=len(data) or result[-1]['id']!=last_id:
            raise ValueError('Index block disagrees with directory')
        return result

    def find(self, prefix):
        matches=[]
        for entry in self.entries[bisect.bisect_left(self.last_ids,prefix):]:
            for row in self.block(entry):
                if row['id']==prefix or row['id'].startswith(prefix+'|'):
                    matches.append(row)
                elif row['id']>prefix:
                    return matches
        return matches


def locate(protein, genome, code):
    # Allow an initiating Met translated from an alternative start codon only
    # at the first residue. The remainder must match exactly, without gaps.
    locations=[]
    for strand, nucleotides in [(1,Seq(genome)),(-1,Seq(genome).reverse_complement())]:
        for frame in range(3):
            translated=str(nucleotides[frame:len(nucleotides)-(len(nucleotides)-frame)%3].translate(table=code))
            probe=protein[1:] if protein.startswith('M') else protein
            cursor=0
            while True:
                found=translated.find(probe,cursor)
                if found<0:
                    break
                cursor=found+1
                start=found-1 if protein.startswith('M') else found
                if start<0:
                    continue
                lo=frame+start*3
                hi=lo+len(protein)*3
                coding=nucleotides[lo:hi]
                if protein.startswith('M'):
                    from Bio.Data import CodonTable
                    if str(coding[:3]) not in CodonTable.unambiguous_dna_by_id[code].start_codons:
                        continue
                if strand==1:
                    begin,end=lo+1,hi
                else:
                    begin,end=len(genome)-hi+1,len(genome)-lo
                locations.append({'start_1based':begin,'end_inclusive':end,'strand':strand,
                                  'genetic_code':code,'match':'exact protein; initiator Met rule allowed'})
    return locations


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('proteins',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    manifest_path=args.output/'ranges.json'
    manifest=json.loads(manifest_path.read_text()) if manifest_path.exists() else {'sources':{},'ranges':{}}
    manifest.update(command=[sys.executable,*sys.argv],script_sha256=sha256(Path(__file__)),
                    proteins_sha256=sha256(args.proteins),zstandard=zstandard.__version__,status='running')
    index=Ffx(Ranges('index',manifest))
    remote=RemoteFile(Ranges('fasta',manifest))
    reader=bgzf.BgzfReader(fileobj=remote,mode='rb')
    proteins=list(SeqIO.parse(args.proteins,'fasta'))
    results=[]
    try:
        for identifier in sorted({p.id.split('|')[0] for p in proteins}):
            entries=index.find(identifier)
            if len(entries)!=1:
                raise ValueError(f'{identifier}: expected one genome, got {len(entries)}')
            entry=entries[0]
            n,bases,width=entry['sequence_length'],entry['line_bases'],entry['line_width']
            if n>100_000_000:
                raise ValueError('Genome exceeds bounded retrieval limit')
            reader.seek(entry['virtual_offset'])
            sequence=reader.read(n+((n-1)//bases)*(width-bases)).replace(b'\n',b'').replace(b'\r',b'').decode()
            if len(sequence)!=n or set(sequence.upper())-set('ACGTRYSWKMBDHVN'):
                raise ValueError('Retrieved genome length/alphabet invalid')
            genome_path=args.output/f'{identifier}.fna'
            genome_path.write_text(f'>{entry["header"]}\n{sequence}\n')
            api_path=args.output/f'{identifier}.json'
            if not api_path.exists():
                response=requests.get('https://meta-virome.org/api/uvigs/'+identifier,timeout=60)
                response.raise_for_status()
                api_path.write_bytes(response.content)
            metadata=json.loads(api_path.read_text())
            if int(metadata['uvig_length'])!=n or metadata['uvig']!=identifier:
                raise ValueError('Metadata and sequence index disagree')
            code=int(metadata['genetic_code'])
            matching=[p for p in proteins if p.id.split('|')[0]==identifier]
            results.append({'genome':identifier,'index_entry':entry,'genome_sha256':sha256(genome_path),
                'metadata_sha256':sha256(api_path),'metadata_url':'https://meta-virome.org/api/uvigs/'+identifier,
                'retrieved_utc':utc(),'proteins':[{'id':p.id,'length':len(p),
                    'locations':locate(str(p.seq).rstrip('*'),sequence,code)} for p in matching]})
            atomic_json(args.output/'loci.json',results)
            atomic_json(manifest_path,manifest)
            print(identifier,n,'bp',len(matching),'proteins',flush=True)
        manifest.update(status='complete',completed_utc=utc(),genomes=len(results),
                        bytes_in_retained_ranges=sum(row['bytes'] for row in manifest['ranges'].values()),
                        loci_sha256=sha256(args.output/'loci.json'))
    except BaseException as error:
        manifest.update(status='interrupted',error=repr(error))
        raise
    finally:
        atomic_json(manifest_path,manifest)


if __name__=='__main__':
    main()
