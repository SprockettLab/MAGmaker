#!/usr/bin/env python3
"""Hybrid assemblies get binning rows, and are binned like a short-read row.

Part 1 runs generate_binning_config for real on a two-sample sheet (John
with long reads, Paul without) and checks the table it writes:
  - hybrid on: one John__metaspades_hybrid row, pointing at the hybrid
    contigs, no reads of its own, nothing for Paul;
  - hybrid off: no hybrid row at all, the table otherwise unchanged.
Part 2 dry-runs Snakefile-bin with a hybrid row and checks that it is
mapped with prototype reads, gets no long-read mapping and no circular
candidates, runs SemiBin2 in its default mode, and that geNomad reads the
hybrid contigs once and never the sample's megahit contigs under that row.

Every check here fails on the code before this change (no hybrid row is
written; geNomad is given megahit for the row).

Needs snakemake on PATH. Run from the repository root:

    python3 resources/test/test_hybrid_binning.py
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
ROW = 'John__metaspades_hybrid'


def check(name, condition, detail=''):
    print(("PASS  " if condition else "FAIL  ") + name)
    if not condition:
        FAILURES.append(name)
        if detail:
            print("      " + str(detail)[-1500:])


def workdir():
    tmp = tempfile.mkdtemp(prefix='magmaker_hybrid_bin_', dir=os.environ.get('TMPDIR'))
    os.symlink(os.path.join(REPO_ROOT, 'resources'), os.path.join(tmp, 'resources'))
    open(os.path.join(tmp, 'host.fna'), 'w').write('>c1\nACGTACGTACGGTTAACCGG\n')
    open(os.path.join(tmp, 'lambda.fa'), 'w').write('>NC_001416.1\nACGTTGCA\n')
    meta = os.path.join(tmp, 'metadata.txt')
    with open(meta, 'w') as fh:
        fh.write('Sample\tSequencing_Run\tPlatform\tR1_fp\tR2_fp\n')
        for s in ('John', 'Paul'):
            fh.write(f'{s}\tRun_1\tillumina\tresources/test/test_reads/{s}_R1.fastq.gz\t'
                     f'resources/test/test_reads/{s}_R2.fastq.gz\n')
        fh.write(f'John\tONT1\tont\t{tmp}/John.ont.fastq.gz\tNA\n')
    open(os.path.join(tmp, 'John.ont.fastq.gz'), 'w').close()
    return tmp, meta


def overrides(tmp, meta, hybrid, extra=None):
    over = {'metadata': meta, 'long_read_assemblers': ['flye'],
            'hybrid_assemblers': ['metaspades_hybrid'] if hybrid else [],
            'params': {'long_reads': {'control_strand_fasta': os.path.join(tmp, 'lambda.fa')},
                       'prototypes': {'n': 2}},
            'host_filter': {'boardercollie': False,
                            'genome': os.path.join(tmp, 'host.fna'), 'db_dir': tmp}}
    over.update(extra or {})
    path = os.path.join(tmp, 'over.yaml')
    with open(path, 'w') as fh:
        yaml.safe_dump(over, fh)
    return path


def binning_table(hybrid):
    """Run generate_binning_config and return its rows, or None on failure."""
    tmp, meta = workdir()
    proto = os.path.join(tmp, 'output/prototype_selection/prototype_selection')
    os.makedirs(proto)
    with open(os.path.join(proto, 'selected_prototypes.yaml'), 'w') as fh:
        yaml.safe_dump({2: ['John', 'Paul']}, fh)
    p = subprocess.run(['snakemake', 'output/config/auto_binning.txt',
                        '-s', os.path.join(REPO_ROOT, 'Snakefile'),
                        '--cores', '1', '-d', tmp,
                        '--allowed-rules', 'generate_binning_config',
                        '--rerun-triggers', 'mtime',
                        '--configfile', os.path.join(REPO_ROOT, 'config.yaml'),
                        overrides(tmp, meta, hybrid)],
                       capture_output=True, text=True)
    out = p.stdout + p.stderr
    table = os.path.join(tmp, 'output/config/auto_binning.txt')
    rows = None
    if p.returncode == 0 and os.path.exists(table):
        with open(table) as fh:
            rows = {r['Sample']: r for r in csv.DictReader(fh, delimiter='\t')}
    shutil.rmtree(tmp)
    return rows, out


def jobs(out):
    found = {}
    for block in re.split(r"\n(?=\[[^\]]+\]\n(?:local)?rule )", out):
        m = re.match(r"\s*\[[^\]]+\]\n(?:local)?rule (\w+):", block)
        if m:
            found.setdefault(m.group(1), []).append(block)
    return found


def row_of(block):
    m = re.search(r"^\s*wildcards: .*?contig_sample=([^,\s]+)", block, re.M)
    return m.group(1) if m else None


def main():
    if not shutil.which('snakemake'):
        sys.exit("Needs snakemake on PATH")

    # Part 1: the binning table generate_binning_config writes.
    on, out = binning_table(hybrid=True)
    check("generate_binning_config runs with hybrid on", on is not None, out)
    if on is not None:
        check("hybrid on: a John__metaspades_hybrid row", ROW in on, sorted(on))
        r = on.get(ROW, {})
        check("hybrid row points at the hybrid contigs",
              r.get('Contigs') == 'output/assemble/metaspades_hybrid/John.contigs.fasta', r)
        check("hybrid row contributes no reads of its own", r.get('Read_Groups', 'x') == '', r)
        check("hybrid row is binned with the prototypes", r.get('Contig_Groups') == 'Group1', r)
        check("hybrid row names its sample and assembler",
              r.get('Source_Sample') == 'John' and r.get('Assembler') == 'metaspades_hybrid', r)
        check("no hybrid row for Paul, who has no long reads",
              not any(k.startswith('Paul__') for k in on), sorted(on))
        check("rows: John, Paul, John__flye, John__metaspades_hybrid",
              sorted(on) == sorted(['John', 'Paul', 'John__flye', ROW]), sorted(on))
    off, out = binning_table(hybrid=False)
    check("generate_binning_config runs with hybrid off", off is not None, out)
    if off is not None:
        check("hybrid off: no hybrid row, everything else unchanged",
              sorted(off) == sorted(['John', 'Paul', 'John__flye']), sorted(off))

    # Part 2: how the binning stage plans a hybrid row.
    tmp, meta = workdir()
    binning = os.path.join(tmp, 'binning.txt')
    with open(binning, 'w') as fh:
        fh.write('Sample\tContigs\tRead_Groups\tContig_Groups\tSource_Sample\tAssembler\n')
        fh.write('John\toutput/assemble/megahit/John.contigs.fasta\tGroup1\tGroup1\tJohn\tmegahit\n')
        fh.write('Paul\toutput/assemble/megahit/Paul.contigs.fasta\tGroup1\tGroup1\tPaul\tmegahit\n')
        fh.write(f'{ROW}\toutput/assemble/metaspades_hybrid/John.contigs.fasta\t\tGroup1\tJohn\tmetaspades_hybrid\n')
    over = overrides(tmp, meta, hybrid=True, extra={'binning': binning})
    p = subprocess.run(['snakemake', '-s', os.path.join(REPO_ROOT, 'Snakefile-bin'),
                        '--configfile', os.path.join(REPO_ROOT, 'config.yaml'), over,
                        '--cores', '1', '-d', tmp, '-n', '-p', 'rename_mags', 'virus_all'],
                       capture_output=True, text=True)
    out = p.stdout + p.stderr
    check("binning stage plans with a hybrid row", p.returncode == 0, out)
    if p.returncode == 0:
        j = jobs(out)
        hy = lambda rule: [b for b in j.get(rule, []) if row_of(b) == ROW]
        mapped = [rule for rule, blocks in j.items()
                  if rule.startswith(('map_', 'minimap2', 'bowtie2'))
                  and rule != 'map_long_reads_to_contigs'
                  and any(row_of(b) == ROW for b in blocks)]
        check("prototype reads are mapped to the hybrid contigs", bool(mapped), sorted(j))
        check("no long-read mapping for the hybrid row", not hy('map_long_reads_to_contigs'))
        check("no circular candidates for the hybrid row", not hy('circular_candidates'))
        sb = hy('run_semibin2')
        check("SemiBin2 bins the hybrid contigs in its default mode",
              len(sb) == 1 and '--sequencing-type' not in sb[0]
              and 'output/assemble/metaspades_hybrid/John.contigs.fasta' in sb[0], sb)
        gn = hy('run_genomad')
        check("geNomad reads the hybrid contigs, once",
              len(gn) == 1 and 'output/assemble/metaspades_hybrid/John.contigs.fasta' in gn[0], gn)
        check("geNomad never reads megahit under the hybrid row",
              not any('assembler=megahit' in b for b in gn), gn)
        check("the short-read row John still gets megahit geNomad",
              any(row_of(b) == 'John' and 'assembler=megahit' in b
                  for b in j.get('run_genomad', [])))
    shutil.rmtree(tmp)


if __name__ == '__main__':
    main()
    print("\n%d failure(s)" % len(FAILURES))
    sys.exit(1 if FAILURES else 0)
