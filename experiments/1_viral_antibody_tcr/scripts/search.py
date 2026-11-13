#!/usr/bin/env python3
"""Logged protein/profile/translated searches; no automatic natural-capture calls."""
import argparse
import csv
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import pyhmmer
from Bio import SeqIO

from acquire import ROOT, sha256, utc

BLAST_DIR = Path('/home/bizon/anaconda3/envs/imprint/bin')
BLAST_FIELDS = ['qseqid','sseqid','pident','length','qstart','qend','sstart','send','evalue','bitscore','qlen','slen','qcovhsp']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('method', choices=['hmm','blastp','tblastn','blastn'])
    parser.add_argument('query', type=Path)
    parser.add_argument('target', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--cpus',type=int,default=16)
    parser.add_argument('--effective-target-count',type=int)
    parser.add_argument('--effective-db-residues',type=int)
    args = parser.parse_args()
    assert 1 <= args.cpus <= 16
    args.output.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    run = {'started_utc':utc(), 'command':[sys.executable,*sys.argv],
           'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
           'script_sha256':sha256(Path(__file__)),
           'inputs':{str(p):sha256(p) for p in [args.query,args.target]},
           'method':args.method, 'threads':args.cpus, 'status':'running'}
    manifest = args.output.with_suffix('.run.json')
    manifest.write_text(json.dumps(run,indent=2)+'\n')
    try:
        if args.method == 'hmm':
            with pyhmmer.plan7.HMMFile(str(args.query)) as handle:
                hmms = list(handle)
            with pyhmmer.easel.SequenceFile(str(args.target),digital=True,
                                            alphabet=pyhmmer.easel.Alphabet.amino()) as handle:
                sequences = handle.read_block()
            count = len(sequences)
            z = args.effective_target_count or count
            run.update(pyhmmer=pyhmmer.__version__, targets=count, profiles=len(hmms),
                       parameters={'E':100,'domE':100,'Z':z,'retain_independent_domain_E':0.01})
            with args.output.open('w') as handle:
                writer=csv.writer(handle,delimiter='\t')
                writer.writerow(['query','target','sequence_E','sequence_score','domain_iE','domain_score',
                                 'bias','ali_from','ali_to','hmm_from','hmm_to','target_length','model_length',
                                 'hmm_alignment','target_alignment'])
                for hits in pyhmmer.hmmer.hmmsearch(hmms,sequences,cpus=args.cpus,E=100,domE=100,Z=z):
                    for hit in hits:
                        for domain in hit.domains:
                            if domain.i_evalue > 0.01:
                                continue
                            a = domain.alignment
                            name=lambda x:x.decode() if isinstance(x,bytes) else x
                            writer.writerow([name(hits.query.name),name(hit.name),hit.evalue,hit.score,
                                             domain.i_evalue,domain.score,domain.bias,a.target_from,a.target_to,
                                             a.hmm_from,a.hmm_to,hit.length,hits.query.M,
                                             a.hmm_sequence,a.target_sequence])
        else:
            nucleotide = args.method in ['tblastn','blastn']
            db = ROOT / 'data/cache/blast' / (sha256(args.target)[:16] + ('_nt' if nucleotide else '_aa'))
            db.parent.mkdir(parents=True,exist_ok=True)
            make = [str(BLAST_DIR/'makeblastdb'),'-in',str(args.target),'-dbtype','nucl' if nucleotide else 'prot',
                    '-out',str(db)]
            if not db.with_suffix('.nin' if nucleotide else '.pin').exists():
                subprocess.run(make,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            command = [str(BLAST_DIR/args.method),'-query',str(args.query),'-db',str(db),'-evalue','0.01',
                       '-num_threads',str(args.cpus),'-max_target_seqs','1000000',
                       '-outfmt','6 '+' '.join(BLAST_FIELDS),'-out',str(args.output)]
            if args.method in ['blastp','tblastn']:
                command += ['-seg','yes']
            if args.effective_db_residues:
                command += ['-dbsize',str(args.effective_db_residues)]
            run.update(search_command=command,index_command=make,
                       version=subprocess.check_output([str(BLAST_DIR/args.method),'-version'],text=True),
                       binary_sha256=sha256(BLAST_DIR/args.method))
            with args.output.with_suffix('.stderr.txt').open('w') as errors:
                subprocess.run(command,check=True,stdout=subprocess.PIPE,stderr=errors)
            run['columns']=BLAST_FIELDS
        run.update(status='complete', exit_status=0, output_sha256=sha256(args.output))
    except Exception as error:
        run.update(status='failed',error=repr(error),exit_status=1)
        raise
    finally:
        run.update(finished_utc=utc(),seconds=time.monotonic()-started,
                   process_maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   children_maxrss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
        manifest.write_text(json.dumps(run,indent=2)+'\n')
    print(args.output,run['status'],round(run['seconds'],2),'seconds',flush=True)


if __name__ == '__main__':
    main()
