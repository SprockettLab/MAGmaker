#!/usr/bin/env python3
"""long_read_assembler_exclude: a long-read sample that skips flye and myloasm
but keeps its long-read host removal and its hybrid assembly.

Two samples, John and Paul, both with short and long reads, flye and hybrid
metaSPAdes on.
  - Stage 1 dry run (`all`) with Paul excluded: flye is planned for John only;
    long-read host removal and hybrid assembly for both.
  - The same with nothing excluded: flye for both (the control, so the first
    check cannot pass by flye being absent altogether).
  - generate_binning_config, run for real with Paul excluded: John__flye,
    John__metaspades_hybrid and Paul__metaspades_hybrid rows, no Paul__flye.
  - A name with no long reads, or a bare string, stops the run at load with a
    message naming the problem.

Needs snakemake (with scikit-bio) on PATH. Run from the repository root:

    python3 resources/test/test_long_read_assembler_exclude.py
"""

import csv
import os
import re
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


def workdir(exclude):
    tmp = tempfile.mkdtemp(prefix='magmaker_lr_exclude_', dir=os.environ.get('TMPDIR'))
    os.symlink(os.path.join(REPO_ROOT, 'resources'), os.path.join(tmp, 'resources'))
    open(os.path.join(tmp, 'host.fna'), 'w').write('>c1\nACGTACGTACGGTTAACCGG\n')
    open(os.path.join(tmp, 'lambda.fa'), 'w').write('>NC_001416.1\nACGTTGCA\n')
    meta = os.path.join(tmp, 'metadata.txt')
    with open(meta, 'w') as fh:
        fh.write('Sample\tSequencing_Run\tPlatform\tR1_fp\tR2_fp\n')
        for s in ('John', 'Paul'):
            fh.write(f'{s}\tRun_1\tillumina\tresources/test/test_reads/{s}_R1.fastq.gz\t'
                     f'resources/test/test_reads/{s}_R2.fastq.gz\n')
        for s in ('John', 'Paul'):
            fh.write(f'{s}\tONT1\tont\t{tmp}/{s}.ont.fastq.gz\tNA\n')
            open(os.path.join(tmp, f'{s}.ont.fastq.gz'), 'w').close()
    over = {'metadata': meta, 'long_read_assemblers': ['flye'],
            'hybrid_assemblers': ['metaspades_hybrid'],
            'long_read_assembler_exclude': exclude,
            'params': {'long_reads': {'control_strand_fasta': os.path.join(tmp, 'lambda.fa')},
                       'prototypes': {'n': 2}},
            'host_filter': {'boardercollie': False,
                            'genome': os.path.join(tmp, 'host.fna'), 'db_dir': tmp}}
    path = os.path.join(tmp, 'over.yaml')
    with open(path, 'w') as fh:
        yaml.safe_dump(over, fh)
    return tmp, path


def snakemake(tmp, over, *args):
    p = subprocess.run(['snakemake', *args, '-s', os.path.join(REPO_ROOT, 'Snakefile'),
                        '--cores', '1', '-d', tmp,
                        '--configfile', os.path.join(REPO_ROOT, 'config.yaml'), over],
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def planned(out, rule):
    """Samples a printed dry run plans `rule` for."""
    found = set()
    for block in re.split(r"\n(?=\[[^\]]+\]\n(?:local)?rule )", out):
        m = re.match(r"\s*\[[^\]]+\]\n(?:local)?rule (\w+):", block)
        if m and m.group(1) == rule:
            w = re.search(r"^\s*wildcards: .*?sample=([^,\s]+)", block, re.M)
            if w:
                found.add(w.group(1))
    return found


def main():
    if not shutil.which('snakemake'):
        sys.exit("Needs snakemake on PATH")

    # Stage 1, Paul excluded.
    tmp, over = workdir(['Paul'])
    rc, out = snakemake(tmp, over, 'all', '-n')
    check("stage 1 plans with Paul excluded", rc == 0, out)
    if rc == 0:
        check("flye planned for John only", planned(out, 'flye') == {'John'},
              planned(out, 'flye'))
        check("long-read host removal still planned for both",
              planned(out, 'host_filter_ont') == {'John', 'Paul'},
              planned(out, 'host_filter_ont'))
        check("hybrid assembly still planned for both",
              planned(out, 'metaspades_hybrid') == {'John', 'Paul'},
              planned(out, 'metaspades_hybrid'))
    shutil.rmtree(tmp)

    # Stage 1, nothing excluded: the control.
    tmp, over = workdir([])
    rc, out = snakemake(tmp, over, 'all', '-n')
    check("stage 1 plans with nothing excluded", rc == 0, out)
    if rc == 0:
        check("control: flye planned for both", planned(out, 'flye') == {'John', 'Paul'},
              planned(out, 'flye'))
    shutil.rmtree(tmp)

    # generate_binning_config, Paul excluded, run for real.
    tmp, over = workdir(['Paul'])
    proto = os.path.join(tmp, 'output/prototype_selection/prototype_selection')
    os.makedirs(proto)
    with open(os.path.join(proto, 'selected_prototypes.yaml'), 'w') as fh:
        yaml.safe_dump({2: ['John', 'Paul']}, fh)
    rc, out = snakemake(tmp, over, 'output/config/auto_binning.txt',
                        '--allowed-rules', 'generate_binning_config',
                        '--rerun-triggers', 'mtime')
    table = os.path.join(tmp, 'output/config/auto_binning.txt')
    check("generate_binning_config runs with Paul excluded",
          rc == 0 and os.path.exists(table), out)
    if rc == 0 and os.path.exists(table):
        with open(table) as fh:
            rows = sorted(r['Sample'] for r in csv.DictReader(fh, delimiter='\t'))
        check("rows: John, Paul, John__flye and both hybrid rows, no Paul__flye",
              rows == sorted(['John', 'Paul', 'John__flye',
                              'John__metaspades_hybrid', 'Paul__metaspades_hybrid']), rows)
    shutil.rmtree(tmp)

    # Bad input stops the run at load, naming the problem.
    tmp, over = workdir(['Ringo'])
    rc, out = snakemake(tmp, over, 'all', '-n')
    check("a sample with no long reads is refused",
          rc != 0 and 'long_read_assembler_exclude' in out and 'Ringo' in out, out)
    shutil.rmtree(tmp)
    tmp, over = workdir('Paul')
    rc, out = snakemake(tmp, over, 'all', '-n')
    check("a bare string is refused",
          rc != 0 and 'must be a list' in out, out)
    shutil.rmtree(tmp)


if __name__ == '__main__':
    main()
    print("\n%d failure(s)" % len(FAILURES))
    sys.exit(1 if FAILURES else 0)
