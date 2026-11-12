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


def boardercollie_mode():
    """With BoarderCollie on, the long-read index is BoarderCollie's shared
    one: built once per host key under index_dir and linked into the run
    folder; a second run folder with the same host builds no index."""
    bc_dir = os.environ.get('BOARDERCOLLIE_DIR')
    if not bc_dir:
        check("BOARDERCOLLIE_DIR set for the BoarderCollie-mode test", False)
        return
    rng = random.Random(5)
    seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    base = tempfile.mkdtemp(prefix='magmaker_lr_bc_', dir=os.environ.get('TMPDIR'))
    host = seq(100000)
    host_fa = os.path.join(base, 'host1.fna')
    with open(host_fa, 'w') as fh:
        fh.write(">chr1\n" + host + "\n")
    ont = os.path.join(base, 'S.ont.fastq.gz')
    with gzip.open(ont, 'wt') as fh:
        for k, v in {"r_host": host[1000:7000], "r_clean": seq(5000)}.items():
            fh.write(f"@{k}\n{v}\n+\n{'5' * len(v)}\n")
    shared = {'library_dir': os.path.join(base, 'lib'), 'index_dir': os.path.join(base, 'idx'),
              'work_dir': os.path.join(base, 'work'),
              'hosts': [{'name': 'H', 'genomes': [{'fasta': host_fa, 'tax_id': 9606}]}],
              'settings': {'download_only': True}}
    logs = []
    for run in ('runA', 'runB'):
        tmp = os.path.join(base, run)
        os.makedirs(tmp)
        os.symlink(os.path.join(REPO_ROOT, 'resources'), os.path.join(tmp, 'resources'))
        meta = os.path.join(tmp, 'metadata.txt')
        with open(meta, 'w') as fh:
            fh.write("Sample\tSequencing_Run\tPlatform\tR1_fp\tR2_fp\n")
            fh.write("S\tRun_1\tillumina\tresources/test/test_reads/John_R1.fastq.gz\t"
                     "resources/test/test_reads/John_R2.fastq.gz\n")
            fh.write(f"S\tONT1\tont\t{ont}\tNA\n")
        with open(os.path.join(tmp, 'over.yaml'), 'w') as fh:
            yaml.safe_dump({'metadata': meta,
                            'host_filter': {'boardercollie': True, 'boardercollie_dir': bc_dir},
                            'boardercollie': shared}, fh)
        p = subprocess.run(['snakemake', '-s', os.path.join(REPO_ROOT, 'Snakefile'),
                            '--configfile', os.path.join(REPO_ROOT, 'config.yaml'),
                            os.path.join(tmp, 'over.yaml'), '--cores', '2', '-d', tmp,
                            # no small-file checksums, as for a real (gigabyte) host index
                            '--max-checksum-file-size', '0',
                            'output/qc/host_filter_ont/S.counts.tsv'],
                           env={**os.environ, 'BOARDERCOLLIE_RUN_ID': 'test'},
                           capture_output=True, text=True)
        logs.append(p.stdout + p.stderr)
        check(f"{run}: long-read host removal runs with BoarderCollie on", p.returncode == 0,
              logs[-1])
        if p.returncode != 0:
            return
        kept = fastq_names(os.path.join(tmp, 'output', 'qc', 'host_filter_ont', 'nonhost',
                                        'S.ont.fastq.gz'))
        check(f"{run}: host read removed, clean read kept", kept == ['r_clean'], kept)
        shared_mmi = [os.path.join(r, f) for r, _, fs in os.walk(os.path.join(base, 'idx', 'minimap2'))
                      for f in fs if f.endswith('.mmi')]
        linked = os.path.join(tmp, 'output', 'qc', 'host_reference', 'H.map-ont.mmi')
        check(f"{run}: one shared index, linked into the run folder",
              len(shared_mmi) == 1 and os.path.exists(linked)
              and os.stat(linked).st_ino == os.stat(shared_mmi[0]).st_ino, shared_mmi)
    check("the first run builds the shared index, not MAGmaker's own",
          'rule boardercollie_minimap2_index:' in logs[0]
          and 'rule host_minimap2_index:' not in logs[0])
    first = os.path.join(base, 'runA')
    again = subprocess.run(['snakemake', '-s', os.path.join(REPO_ROOT, 'Snakefile'),
                            '--configfile', os.path.join(REPO_ROOT, 'config.yaml'),
                            os.path.join(first, 'over.yaml'), '--cores', '2', '-d', first,
                            '--max-checksum-file-size', '0', '-n',
                            'output/qc/host_filter_ont/S.counts.tsv'],
                           env={**os.environ, 'BOARDERCOLLIE_RUN_ID': 'test'},
                           capture_output=True, text=True)
    check("the first run folder has nothing to redo after the second ran",
          again.returncode == 0 and 'Nothing to be done' in again.stdout + again.stderr,
          (again.stdout + again.stderr)[-1500:])
    check("the second run folder builds no index",
          'rule boardercollie_minimap2_index:' not in logs[1]
          and 'rule host_minimap2_index:' not in logs[1]
          and 'rule boardercollie_minimap2_out_H:' in logs[1], logs[1][-1500:])
    shutil.rmtree(base)


if __name__ == '__main__':
    main()
    boardercollie_mode()
    print("\n%d failure(s)" % len(FAILURES))
    sys.exit(1 if FAILURES else 0)
