#!/usr/bin/env python3
"""Runs the long-read host removal rules end to end through Snakemake.

A random host genome, and Nanopore-like reads cut from it and from random
sequence, so each read's host fraction is known by construction. Runs
host_minimap2_index and host_filter_ont with the tools on PATH (no conda):
needs snakemake, minimap2, samtools, seqkit, pandas.

Run from the repository root:

    python3 resources/test/test_long_read_host.py
"""

import gzip
import os
import random
import shutil
import subprocess
import sys
import tempfile

import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
FAILURES = []


def check(name, condition, detail=''):
    print(("PASS  " if condition else "FAIL  ") + name)
    if not condition:
        FAILURES.append(name)
        if detail:
            print("      " + str(detail)[-1500:])


def fastq_names(path):
    with gzip.open(path, 'rt') as fh:
        return [l[1:].split()[0] for i, l in enumerate(fh) if i % 4 == 0]


def main():
    missing = [t for t in ('snakemake', 'minimap2', 'samtools', 'seqkit') if not shutil.which(t)]
    if missing:
        sys.exit("Needs on PATH: " + ", ".join(missing))
    rng = random.Random(11)
    seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    tmp = tempfile.mkdtemp(prefix='magmaker_lr_host_',
                           dir=os.environ.get('TMPDIR'))
    os.symlink(os.path.join(REPO_ROOT, 'resources'), os.path.join(tmp, 'resources'))
    host = seq(200000)
    with open(os.path.join(tmp, 'host.fna'), 'w') as fh:
        fh.write(">chr1\n" + host + "\n")
    reads = {
        "r_host": host[10000:16000],                       # 100%
        "r_host60": host[50000:53000] + seq(2000),         # 60%
        "r_part40": seq(3000) + host[90000:92000],         # 40%: kept, partial
        "r_clean1": seq(5000),
        "r_clean2": seq(4000),
    }
    ont = os.path.join(tmp, 'S.ont.fastq.gz')
    with gzip.open(ont, 'wt') as fh:
        for k, v in reads.items():
            fh.write(f"@{k} runid=x\n{v}\n+\n{'5' * len(v)}\n")
    meta = os.path.join(tmp, 'metadata.txt')
    with open(meta, 'w') as fh:
        fh.write("Sample\tSequencing_Run\tPlatform\tR1_fp\tR2_fp\n")
        fh.write("S\tRun_1\tillumina\tresources/test/test_reads/John_R1.fastq.gz\t"
                 "resources/test/test_reads/John_R2.fastq.gz\n")
        fh.write(f"S\tONT1\tont\t{ont}\tNA\n")
    with open(os.path.join(tmp, 'over.yaml'), 'w') as fh:
        yaml.safe_dump({'metadata': meta,
                        'host_filter': {'boardercollie': False,
                                        'genome': os.path.join(tmp, 'host.fna'),
                                        'db_dir': tmp}}, fh)
    target = 'output/qc/host_filter_ont/S.counts.tsv'
    p = subprocess.run(['snakemake', '-s', os.path.join(REPO_ROOT, 'Snakefile'),
                        '--configfile', os.path.join(REPO_ROOT, 'config.yaml'),
                        os.path.join(tmp, 'over.yaml'), '--cores', '2', '-d', tmp, target],
                       capture_output=True, text=True)
    check("host_minimap2_index and host_filter_ont run", p.returncode == 0, p.stdout + p.stderr)
    if p.returncode != 0:
        return
    out = os.path.join(tmp, 'output', 'qc', 'host_filter_ont')
    kept = sorted(fastq_names(os.path.join(out, 'nonhost', 'S.ont.fastq.gz')))
    check("nonhost reads are exactly the reads under 50% host",
          kept == ['r_clean1', 'r_clean2', 'r_part40'], kept)
    ids = sorted(gzip.open(os.path.join(out, 'host', 'S.host_ids.txt.gz'), 'rt').read().split())
    check("host read IDs", ids == ['r_host', 'r_host60'], ids)
    bam_reads = subprocess.run(['samtools', 'view', os.path.join(out, 'host', 'S.bam')],
                               capture_output=True, text=True).stdout.split('\n')
    check("host BAM holds the alignments, not unaligned reads",
          {l.split('\t')[0] for l in bam_reads if l} == {'r_host', 'r_host60', 'r_part40'})
    c = dict(l.rstrip('\n').split('\t', 1) for l in open(os.path.join(out, 'S.counts.tsv')))
    check("counts", c['reads_in'] == '5' and c['reads_host'] == '2' and c['reads_kept'] == '3'
          and c['reads_kept_with_partial_host_alignment'] == '1', c)
    check("index is one part (no split index)",
          'multi-part' not in open(os.path.join(tmp, 'output', 'logs', 'qc', 'host_filter_ont',
                                                'host_minimap2_index.log')).read())

    # Assembly is planned from the nonhost reads, once per long-read assembler.
    drafts = ['output/assemble/flye/S.draft.fasta', 'output/assemble/myloasm/S.draft.fasta']
    p = subprocess.run(['snakemake', '-s', os.path.join(REPO_ROOT, 'Snakefile'),
                        '--configfile', os.path.join(REPO_ROOT, 'config.yaml'),
                        os.path.join(tmp, 'over.yaml'), '--cores', '2', '-d', tmp, '-n', *drafts],
                       capture_output=True, text=True)
    out = p.stdout + p.stderr
    stats = dict(l.split()[:2] for l in out.split('Job stats:')[-1].splitlines()[3:] if len(l.split()) == 2)
    check("flye and myloasm planned once each, nothing else rerun",
          p.returncode == 0 and stats.get('flye') == '1' and stats.get('myloasm') == '1'
          and stats.get('total') == '2', out[-1500:])
    shutil.rmtree(tmp)


if __name__ == '__main__':
    main()
    print("\n%d failure(s)" % len(FAILURES))
    sys.exit(1 if FAILURES else 0)
