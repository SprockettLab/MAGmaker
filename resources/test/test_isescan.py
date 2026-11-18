#!/usr/bin/env python3
"""Checks the insertion sequence (ISEScan) module.

1. Snakemake dry runs of isescan_all on the test data:
   - off by default: no ISEScan job is planned;
   - isescan=True: one run_isescan per binning row (its own assembly) and
     one make_isescan_summary.

2. make_isescan_summary.py on a built output tree:
   - ISEScan's own column names (pred.py, ISEScan 1.7.3), with a MEGAHIT
     header carrying spaces in seqID, as ISEScan writes it;
   - a contig binned into a MAG is labelled with that MAG; one in no MAG
     gets NA;
   - the same contig name in two rows (MEGAHIT numbers every assembly from
     k141_0) is labelled only in the row whose bin holds it;
   - a contig in MAGs from two mappers gets both, ';'-joined;
   - an empty table (no IS found) adds no rows and is logged;
   - a MAG FASTA that is missing leaves its contigs unlabelled and is logged.

Run from the repository root (part 1 needs snakemake on PATH):

    python3 resources/test/test_isescan.py
"""

import csv
import os
import re
import subprocess
import sys
import tempfile
import types

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
SCRIPT = os.path.join(REPO_ROOT, 'resources', 'scripts', 'make_isescan_summary.py')
CHECKS = []

# ISEScan 1.7.3 .tsv header (pred.py titleLine4raw).
ISESCAN_COLS = ['seqID', 'family', 'cluster', 'isBegin', 'isEnd', 'isLen', 'ncopy4is',
                'start1', 'end1', 'start2', 'end2', 'score', 'irId', 'irLen', 'nGaps',
                'orfBegin', 'orfEnd', 'strand', 'orfLen', 'E-value', 'E-value4copy',
                'type', 'ov', 'tir']


def check(name, ok):
    CHECKS.append((name, ok))
    print(('  ok   ' if ok else '  FAIL ') + name)
    if not ok:
        raise AssertionError(name)


# ── part 1: the plan ────────────────────────────────────────────────────────

def plan(*extra_config):
    with tempfile.TemporaryDirectory(prefix='magmaker_isescan_') as tmp:
        os.symlink(os.path.join(REPO_ROOT, 'resources'), os.path.join(tmp, 'resources'))
        with open(os.path.join(tmp, 'host.fna'), 'w') as handle:
            handle.write('>c1\nACGTACGTACGGTTAACCGGTTAAGGCCAATTGG\n')
        over = os.path.join(tmp, 'over.yaml')
        with open(over, 'w') as handle:
            handle.write('host_filter:\n  boardercollie: false\n'
                         '  db_dir: {0}/host_db\n  genome: {0}/host.fna\n'.format(tmp))
        result = subprocess.run(
            ['snakemake', '-s', os.path.join(REPO_ROOT, 'Snakefile-bin'),
             '--configfile', os.path.join(REPO_ROOT, 'config.yaml'), over,
             '--cores', '1', '-d', tmp, '-n', 'isescan_all',
             '--config', 'binning=resources/config/binning.txt', *extra_config],
            capture_output=True, text=True)
        out = result.stdout + result.stderr
        if result.returncode != 0:
            return None, out
        if 'Job stats' not in out:
            return {}, out
        stats = out[out.index('Job stats'):]
        return ({m.group(1): int(m.group(2))
                 for m in re.finditer(r'^(\w+)\s+(\d+)$', stats, re.M)}, out)


def test_plan():
    print('dry runs:')
    rows = sum(1 for _ in open(os.path.join(REPO_ROOT, 'resources', 'config', 'binning.txt'))) - 1
    counts, out = plan()
    check('the default plan runs (exit 0)', counts is not None or print(out))
    check('off by default: no run_isescan job', counts.get('run_isescan', 0) == 0)
    counts, out = plan('isescan=True')
    check('isescan=True plans (exit 0)', counts is not None or print(out))
    check(f'isescan=True: one run_isescan per binning row ({rows})',
          counts.get('run_isescan') == rows)
    check('isescan=True: one make_isescan_summary', counts.get('make_isescan_summary') == 1)


# ── part 2: the summary ─────────────────────────────────────────────────────

def is_row(seqid, begin, end, typ='c', family='IS3'):
    row = dict.fromkeys(ISESCAN_COLS, '-')
    row.update(seqID=seqid, family=family, cluster=family + '_1', isBegin=str(begin),
               isEnd=str(end), isLen=str(end - begin + 1), ncopy4is='2', type=typ)
    return row


def write_tsv(path, rows, cols):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', newline='') as handle:
        w = csv.DictWriter(handle, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def write_fasta(path, names):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as handle:
        for n in names:
            handle.write('>{0}\nACGT\n'.format(n))


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def run_summary(tmp, pairs, tables, mag_summary):
    os.makedirs(os.path.join(tmp, 'output', 'logs'), exist_ok=True)
    stub = types.SimpleNamespace(
        input=types.SimpleNamespace(tables=tables, mag_summary=mag_summary),
        params=types.SimpleNamespace(pairs=pairs, row_sample={r: r.split('__')[0] for _, r in pairs}),
        output=types.SimpleNamespace(summary=os.path.join(tmp, 'output', 'isescan', 'isescan_summary.tsv')),
        log=[os.path.join(tmp, 'output', 'logs', 'summary.log')])
    saved_err, saved_cwd = sys.stderr, os.getcwd()
    os.chdir(tmp)
    try:
        with open(SCRIPT) as handle:
            exec(compile(handle.read(), SCRIPT, 'exec'), {'snakemake': stub, '__name__': '__main__'})
    finally:
        sys.stderr.close() if sys.stderr is not saved_err else None
        sys.stderr = saved_err
        os.chdir(saved_cwd)
    return read_tsv(stub.output.summary), open(stub.log[0]).read()


def test_summary():
    print('summary:')
    with tempfile.TemporaryDirectory(prefix='magmaker_isescan_sum_') as tmp:
        os.makedirs(os.path.join(tmp, 'output', 'isescan'))
        bins = 'output/selected_bins/{0}/Binette_Fastas/{1}/{2}.fa'
        # Row A: k141_1 in a MAG from both mappers; k141_2 in no MAG.
        write_fasta(os.path.join(tmp, bins.format('minimap2', 'A', 'A_bin1')), ['k141_1', 'k141_9'])
        write_fasta(os.path.join(tmp, bins.format('bowtie2', 'A', 'A_bin1')), ['k141_1'])
        # Row B: has its own k141_1, binned nowhere; its bin holds k141_5.
        write_fasta(os.path.join(tmp, bins.format('minimap2', 'B', 'B_bin3')), ['k141_5'])
        mag_cols = ['MAG_ID', 'New_Name', 'Original_Name', 'Original_Path', 'Sample_ID',
                    'GTDB_Classification', 'MIMAG_Quality']
        mags = [
            dict(MAG_ID='MAG_0001', New_Name='MAG_0001__Bacteroides_x', Original_Name='A_bin1',
                 Original_Path=bins.format('minimap2', 'A', 'A_bin1'), Sample_ID='A',
                 GTDB_Classification='g__Bacteroides', MIMAG_Quality='HQ'),
            dict(MAG_ID='MAG_0002', New_Name='MAG_0002__Bacteroides_x', Original_Name='A_bin1',
                 Original_Path=bins.format('bowtie2', 'A', 'A_bin1'), Sample_ID='A',
                 GTDB_Classification='g__Bacteroides', MIMAG_Quality='MQ'),
            dict(MAG_ID='MAG_0003', New_Name='MAG_0003__Prevotella_y', Original_Name='B_bin3',
                 Original_Path=bins.format('minimap2', 'B', 'B_bin3'), Sample_ID='B',
                 GTDB_Classification='g__Prevotella', MIMAG_Quality='MQ'),
            dict(MAG_ID='MAG_0004', New_Name='MAG_0004__Gone', Original_Name='B_bin4',
                 Original_Path=bins.format('minimap2', 'B', 'B_bin4'), Sample_ID='B',
                 GTDB_Classification='g__Gone', MIMAG_Quality='LQ'),
        ]
        ms = os.path.join(tmp, 'output', 'mag_qc', 'mag_summary.tsv')
        write_tsv(ms, mags, mag_cols)
        ta = os.path.join(tmp, 'output', 'isescan', 'megahit', 'A', 'is.tsv')
        tb = os.path.join(tmp, 'output', 'isescan', 'megahit', 'B', 'is.tsv')
        tc = os.path.join(tmp, 'output', 'isescan', 'megahit', 'C', 'is.tsv')
        write_tsv(ta, [is_row('k141_1 flag=1 multi=12.0 len=5000', 100, 1400),
                       is_row('k141_2 flag=0 multi=3.0 len=2000', 10, 900, typ='p')], ISESCAN_COLS)
        write_tsv(tb, [is_row('k141_1 flag=1 multi=8.0 len=3000', 50, 1300),
                       is_row('k141_5 flag=1 multi=8.0 len=9000', 500, 1800)], ISESCAN_COLS)
        os.makedirs(os.path.dirname(tc))
        open(tc, 'w').close()
        rows, log = run_summary(tmp, [('megahit', 'A'), ('megahit', 'B'), ('megahit', 'C')],
                                [ta, tb, tc], ms)
        check('4 IS copies from 2 non-empty tables', len(rows) == 4)
        got = {(r['Binning_Row'], r['Contig']): r for r in rows}
        check('Contig is the first word of a MEGAHIT seqID',
              ('A', 'k141_1') in got and got[('A', 'k141_1')]['seqID'].startswith('k141_1 flag='))
        check('ISEScan columns are kept as written',
              all(c in rows[0] for c in ISESCAN_COLS) and got[('A', 'k141_2')]['type'] == 'p')
        check('a contig in MAGs from two mappers gets both, joined',
              got[('A', 'k141_1')]['MAG_ID'] == 'MAG_0001;MAG_0002')
        check('a contig in no MAG gets NA', got[('A', 'k141_2')]['MAG_ID'] == 'NA')
        check("row B's own k141_1 is not labelled with row A's MAG",
              got[('B', 'k141_1')]['MAG_ID'] == 'NA')
        check("row B's binned contig gets its MAG and taxonomy",
              got[('B', 'k141_5')]['MAG_ID'] == 'MAG_0003'
              and got[('B', 'k141_5')]['MAG_GTDB_Classification'] == 'g__Prevotella')
        check('Sample_ID and Assembler are filled',
              got[('B', 'k141_5')]['Sample_ID'] == 'B' and got[('B', 'k141_5')]['Assembler'] == 'megahit')
        check('the empty table is logged as no IS elements', 'C (megahit): no IS elements' in log)
        check('the missing MAG FASTA is logged', 'MAG FASTA not found' in log and 'B_bin4' in log)
        check('the log gives the total with its denominator', 'from 2 of 3 tables' in log)


def test_conda_prefix_first():
    """The rule must put $CONDA_PREFIX/bin ahead of PATH before running
    ISEScan. Without it, another environment that reaches the job first
    supplies python3, and ISEScan dies on
    `OSError: libssw.so: cannot open shared object file` because dlopen
    resolves through the running interpreter's RUNPATH. This happened to all
    88 samples of a WF18 run on 2026-11-18. Asserted on the source so a
    rewrite of the rule cannot drop it silently."""
    print('PATH guard:')
    smk = open(os.path.join(REPO_ROOT, 'resources', 'snakefiles',
                            'isescan.smk')).read()
    export = 'export PATH="$CONDA_PREFIX/bin:$PATH"'
    check('run_isescan exports $CONDA_PREFIX/bin onto PATH', export in smk)
    if export in smk:
        check('it comes before isescan.py is invoked',
              smk.index(export) < smk.index('isescan.py --seqfile'))


if __name__ == '__main__':
    if '--summary-only' not in sys.argv:
        test_plan()
    test_conda_prefix_first()
    test_summary()
    print('test_isescan: {0} checks passed'.format(len(CHECKS)))
