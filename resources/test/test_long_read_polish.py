#!/usr/bin/env python3
"""Runs short-read polishing and DNA Control Strand removal end to end.

A synthetic bacterial contig is given to MAGmaker as a Flye draft with three
planted base errors, beside a contig of a synthetic lambda. Paired short reads
are simulated, error-free, from the true contig and go through MAGmaker's own
fastp and host filter. Polishing must restore the true sequence (3 changes);
lambda removal must take out the lambda contig. Flye itself is not run: the
draft is written where rule flye would write it, after the long-read host
filter has run, so Snakemake treats it as up to date.

Needs on PATH: snakemake, minimap2, samtools, seqkit, bowtie2, fastp, pigz,
bwa, polypolish, pypolca. Run from the repository root:

    python3 resources/test/test_long_read_polish.py
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
TOOLS = ('snakemake', 'minimap2', 'samtools', 'seqkit', 'bowtie2', 'fastp', 'pigz',
         'bwa', 'polypolish', 'pypolca')
FAILURES = []


def check(name, condition, detail=''):
    print(("PASS  " if condition else "FAIL  ") + name)
    if not condition:
        FAILURES.append(name)
        if detail:
            print("      " + str(detail)[-2000:])


def snakemake(tmp, *targets):
    return subprocess.run(['snakemake', '-s', os.path.join(REPO_ROOT, 'Snakefile'),
                           '--configfile', os.path.join(REPO_ROOT, 'config.yaml'),
                           os.path.join(tmp, 'over.yaml'), '--cores', '4', '-d', tmp, *targets],
                          capture_output=True, text=True)


def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split()[0]
            out[name] = []
        elif line:
            out[name].append(line)
    return {k: ''.join(v) for k, v in out.items()}


def main():
    missing = [t for t in TOOLS if not shutil.which(t)]
    if missing:
        sys.exit("Needs on PATH: " + ", ".join(missing))
    rng = random.Random(23)
    seq = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    rc = lambda s: s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    tmp = tempfile.mkdtemp(prefix='magmaker_lr_polish_', dir=os.environ.get('TMPDIR'))
    os.symlink(os.path.join(REPO_ROOT, 'resources'), os.path.join(tmp, 'resources'))

    host, truth, lam = seq(100000), seq(30000), seq(48502)
    open(os.path.join(tmp, 'host.fna'), 'w').write('>chr1\n' + host + '\n')
    open(os.path.join(tmp, 'lambda.fa'), 'w').write('>NC_001416.1\n' + lam + '\n')

    # Paired 150 bp reads, 400 bp fragments, about 60x, error-free, both strands.
    with gzip.open(os.path.join(tmp, 'S_R1.fastq.gz'), 'wt') as r1, \
            gzip.open(os.path.join(tmp, 'S_R2.fastq.gz'), 'wt') as r2:
        for i in range(6000):
            s = rng.randrange(0, len(truth) - 400)
            frag = truth[s:s + 400]
            if rng.random() < 0.5:
                frag = rc(frag)
            r1.write(f'@p{i}/1\n{frag[:150]}\n+\n{"I" * 150}\n')
            r2.write(f'@p{i}/2\n{rc(frag)[:150]}\n+\n{"I" * 150}\n')
    with gzip.open(os.path.join(tmp, 'S.ont.fastq.gz'), 'wt') as fh:
        for i in range(5):
            s = rng.randrange(0, len(truth) - 8000)
            fh.write(f'@ont{i}\n{truth[s:s + 8000]}\n+\n{"5" * 8000}\n')

    with open(os.path.join(tmp, 'metadata.txt'), 'w') as fh:
        fh.write('Sample\tSequencing_Run\tPlatform\tR1_fp\tR2_fp\n')
        fh.write(f'S\tRun_1\tillumina\t{tmp}/S_R1.fastq.gz\t{tmp}/S_R2.fastq.gz\n')
        fh.write(f'S\tONT1\tont\t{tmp}/S.ont.fastq.gz\tNA\n')
    with open(os.path.join(tmp, 'over.yaml'), 'w') as fh:
        yaml.safe_dump({'metadata': os.path.join(tmp, 'metadata.txt'),
                        'long_read_assemblers': ['flye'],
                        'host_filter': {'boardercollie': False,
                                        'genome': os.path.join(tmp, 'host.fna'), 'db_dir': tmp},
                        'params': {'long_reads': {'control_strand_fasta': os.path.join(tmp, 'lambda.fa')}}},
                       fh)

    p = snakemake(tmp, 'output/qc/host_filter_ont/S.counts.tsv')
    check("long-read host filter runs", p.returncode == 0, p.stdout + p.stderr)

    # The Flye draft: the true contig with three substitutions, and a lambda contig.
    draft = list(truth)
    planted = [5000, 15000, 25000]
    for i in planted:
        draft[i] = 'A' if draft[i] != 'A' else 'C'
    asm = os.path.join(tmp, 'output', 'assemble', 'flye')
    os.makedirs(os.path.join(asm, 'S_flye'), exist_ok=True)
    with open(os.path.join(asm, 'S.draft.fasta'), 'w') as fh:
        fh.write('>S__flye__contig_1\n' + ''.join(draft) + '\n>S__flye__contig_2\n' + lam[1000:21000] + '\n')
    with open(os.path.join(asm, 'S.contig_info.tsv'), 'w') as fh:
        fh.write('contig\toriginal_name\tlength\tcircular\tdepth\n'
                 f'S__flye__contig_1\tcontig_1\t{len(truth)}\tyes\t30\n'
                 'S__flye__contig_2\tcontig_2\t20000\tno\t5\n')
    open(os.path.join(asm, 'S_flye', 'assembly_info.txt'), 'w').write('#seq_name\n')

    p = snakemake(tmp, 'output/assemble/flye/S.contigs.fasta')
    out = p.stdout + p.stderr
    check("fastp, host filter, polishing and lambda removal run", p.returncode == 0, out)
    if p.returncode == 0:
        check("Flye itself was not rerun", 'rule flye:' not in out)
        final = read_fasta(os.path.join(asm, 'S.contigs.fasta'))
        check("lambda contig removed, bacterial contig kept", list(final) == ['S__flye__contig_1'], list(final))
        check("polishing restored the true sequence",
              final.get('S__flye__contig_1') == truth,
              sum(a != b for a, b in zip(final.get('S__flye__contig_1', ''), truth)))
        rows = {l.split('\t')[0]: l.rstrip('\n').split('\t')
                for l in open(os.path.join(asm, 'S.polish_changes.tsv')) if not l.startswith(('#', 'contig\t'))}
        check("3 changes recorded for the bacterial contig", rows['S__flye__contig_1'][4] == '3',
              rows['S__flye__contig_1'])
        check("draft kept beside the polished contigs", os.path.exists(os.path.join(asm, 'S.draft.fasta')))
        removed = read_fasta(os.path.join(asm, 'S.control_strand.fasta'))
        check("lambda contig kept aside", list(removed) == ['S__flye__contig_2'])
    shutil.rmtree(tmp)


if __name__ == '__main__':
    main()
    print("\n%d failure(s)" % len(FAILURES))
    sys.exit(1 if FAILURES else 0)
