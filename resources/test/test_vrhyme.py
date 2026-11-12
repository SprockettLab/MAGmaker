#!/usr/bin/env python3
"""Checks the two vRhyme scripts: input preparation and the summary tables.

vrhyme_prep.py:
  - headers are truncated to their first token. vRhyme keys scaffolds on the
    ENTIRE header line while the coverage table carries only the contig name,
    so with full headers nothing joins: vRhyme finds no pairs, reports
    "binning complete" with 0 bins and exits 0. WF30 lost a day to this.
  - proviruses are excluded, by geNomad's name OR CheckV's call. A dormant
    prophage's coverage tracks its host, and CheckV reports No for the ones
    geNomad already excised, so the CheckV column alone misses most of them.
  - contigs under the length floor are excluded, and the floor cannot go
    below vRhyme's own 2000.
  - every CheckV quality tier is kept, Low-quality and Not-determined
    included: the pieces of a fragmented phage each score badly alone, and
    joining them is the point.
  - every count reaches prep_params.txt with its denominator.

make_vrhyme_summary.py:
  - one row per vMAG, one per contig, with multi-contig bins counted;
  - circular scaffolds, which vRhyme sets aside before binning, are kept as
    their own vMAGs rather than going missing;
  - several membership files under vRhyme_alternate_bins/ do not make it
    report one iteration's numbers as the run's.

Run from the repository root:

    python3 resources/test/test_vrhyme.py
"""

import os
import sys
import tempfile
import types

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
PREP = os.path.join(REPO_ROOT, 'resources', 'scripts', 'vrhyme_prep.py')
SUMMARY = os.path.join(REPO_ROOT, 'resources', 'scripts', 'make_vrhyme_summary.py')
CHECKS = []


def check(name, ok):
    CHECKS.append((name, ok))
    print(('  ok   ' if ok else '  FAIL ') + name)
    if not ok:
        raise AssertionError(name)


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as handle:
        handle.write(text)


def read_tsv(path):
    with open(path) as handle:
        lines = handle.read().splitlines()
    columns = lines[0].split('\t')
    return [dict(zip(columns, line.split('\t'))) for line in lines[1:]]


def run_script(script, stub):
    saved_out, saved_err = sys.stdout, sys.stderr
    try:
        with open(script) as handle:
            code = compile(handle.read(), script, 'exec')
        exec(code, {'snakemake': stub, '__name__': '__main__'})
    finally:
        sys.stdout, sys.stderr = saved_out, saved_err


# ── vrhyme_prep ─────────────────────────────────────────────────────────────

# geNomad's real header shape, with a description after the name. The
# coverage table has only the first token, which is the whole trap.
VIRAL = [
    # name, length, what CheckV said, why it is or is not selected
    ('k141_1 flag=1 multi=10.0000 len=5000', 5000, ('High-quality', 'No'), True),
    ('k141_2 flag=1 multi=4.0000 len=2500', 2500, ('Low-quality', 'No'), True),
    # Not-determined is kept on purpose: a fragment of a phage scores badly
    # alone and joining it is what binning is for.
    ('k141_3 flag=0 multi=2.0000 len=2000', 2000, ('Not-determined', ''), True),
    # Excised by geNomad: CheckV says No because no host DNA is left, so only
    # the name reveals it.
    ('k141_4|provirus_100_5100 flag=1 multi=9.0 len=5000', 5000,
     ('Medium-quality', 'No'), False),
    # CheckV's own call, on a name without coordinates.
    ('k141_5 flag=1 multi=8.0000 len=4000', 4000, ('Medium-quality', 'Yes'), False),
    # Below vRhyme's floor.
    ('k141_6 flag=1 multi=3.0000 len=1500', 1500, ('Low-quality', 'No'), False),
]


def build(tmp, sample='S1', assembler='megahit'):
    stem = sample + '.contigs'
    summary_dir = os.path.join(tmp, 'output', 'virus', 'genomad', assembler,
                               sample, stem + '_summary')
    fasta = os.path.join(summary_dir, stem + '_virus.fna')
    write(fasta, ''.join('>{0}\n{1}\n'.format(header, 'ACGT' * 10)
                         for header, _, _, _ in VIRAL))
    write(os.path.join(summary_dir, stem + '_virus_summary.tsv'),
          'seq_name\tlength\ttopology\n' +
          ''.join('{0}\t{1}\tNo terminal repeats\n'.format(header.split()[0], length)
                  for header, length, _, _ in VIRAL))
    checkv = os.path.join(tmp, 'output', 'virus', 'checkv', assembler, sample)
    write(os.path.join(checkv, 'quality_summary.tsv'),
          'contig_id\tcheckv_quality\tprovirus\tproviral_start\tproviral_end\n' +
          ''.join('{0}\t{1}\t{2}\t\t\n'.format(header.split()[0], quality, provirus)
                  for header, _, (quality, provirus), _ in VIRAL))
    # The MetaBAT2 jgi table: contigName is the bare contig name.
    coverage = os.path.join(tmp, 'output', 'binning', 'metabat2', 'minimap2',
                            'coverage_tables', sample + '_coverage_table.txt')
    write(coverage,
          'contigName\tcontigLen\ttotalAvgDepth\tS1.bam\tS1.bam-var\n' +
          ''.join('{0}\t{1}\t10.0\t10.0\t1.0\n'.format(header.split()[0].split('|')[0], length)
                  for header, length, _, _ in VIRAL))
    return fasta, checkv, coverage, stem


def prep(tmp, min_length=2000):
    fasta, checkv, coverage, stem = build(tmp)
    out = os.path.join(tmp, 'output', 'virus', 'vrhyme', 'megahit', 'S1')
    os.makedirs(out, exist_ok=True)
    stub = types.SimpleNamespace(
        input=types.SimpleNamespace(coverage=coverage),
        output=types.SimpleNamespace(
            fasta=os.path.join(out, 'input.viral.fna'),
            coverage=os.path.join(out, 'input.coverage.jgi.tsv'),
            params_file=os.path.join(out, 'prep_params.txt')),
        params=types.SimpleNamespace(virus_fna=fasta, checkv_dir=checkv,
                                     min_length=min_length),
        log=[os.path.join(tmp, 'output', 'logs', 'prep.log')])
    run_script(PREP, stub)
    return stub


def params_of(path):
    values = {}
    for line in open(path):
        if line.startswith('#') or '\t' not in line:
            continue
        key, value = line.rstrip('\n').split('\t', 1)
        values[key] = value
    return values


def test_prep():
    print('vrhyme_prep:')
    with tempfile.TemporaryDirectory(prefix='magmaker_prep_') as tmp:
        stub = prep(tmp)
        headers = [line[1:].strip() for line in open(stub.output.fasta)
                   if line.startswith('>')]
        expected = [header.split()[0] for header, _, _, keep in VIRAL if keep]

        check('every header is truncated to its first token, so vRhyme can '
              'join the coverage table at all',
              headers == expected and all(' ' not in h for h in headers))
        check('the selected contigs are the non-provirus ones at or above '
              'the floor', sorted(headers) == sorted(expected))
        check("geNomad's excised provirus is excluded although CheckV says No",
              not any(h.startswith('k141_4') for h in headers))
        check("CheckV's own provirus call excludes a contig too",
              'k141_5' not in headers)
        check('a contig below the floor is excluded', 'k141_6' not in headers)
        check('a Low-quality and a Not-determined contig are both kept',
              'k141_2' in headers and 'k141_3' in headers)

        rows = [line.split('\t')[0] for line in
                open(stub.output.coverage).read().splitlines()[1:]]
        check('the coverage table is subset to exactly those contigs',
              sorted(rows) == sorted(expected))

        values = params_of(stub.output.params_file)
        check('every count is recorded with its denominator',
              values['genomad_virus_sequences'] == str(len(VIRAL))
              and values['selected'] == str(len(expected))
              and values['excluded_provirus'] == '2'
              and values['excluded_short'] == '1'
              and values['fasta_sequences_written'] == str(len(expected))
              and values['coverage_rows_kept'] == str(len(expected))
              and values['coverage_table_rows'] == str(len(VIRAL)))
        check('the params file records the settings and that nothing was '
              'subsampled',
              values['min_length_bp'] == '2000'
              and values['subsampling'] == 'none_all_sequences_read'
              and values['quality_filter'] == 'none_all_tiers'
              and values['fasta_headers'] == 'truncated_to_first_token')
        check('the quality tiers of the selected contigs are recorded',
              values['selected_tier_High-quality'] == '1'
              and values['selected_tier_Low-quality'] == '1'
              and values['selected_tier_Not-determined'] == '1')

    with tempfile.TemporaryDirectory(prefix='magmaker_prep2_') as tmp:
        try:
            prep(tmp, min_length=1000)
        except SystemExit as exit_error:
            check("a floor below vRhyme's own 2000 is refused",
                  'cannot be lowered' in str(exit_error))
        else:
            raise AssertionError('a min_length of 1000 was accepted')


# ── make_vrhyme_summary ─────────────────────────────────────────────────────

def summary(tmp, pairs, row_sample):
    out = os.path.join(tmp, 'output', 'virus')
    os.makedirs(out, exist_ok=True)
    stub = types.SimpleNamespace(
        input=types.SimpleNamespace(done=[]),
        output=types.SimpleNamespace(
            summary=os.path.join(out, 'vrhyme_summary.tsv'),
            membership=os.path.join(out, 'vrhyme_membership.tsv')),
        params=types.SimpleNamespace(pairs=pairs, row_sample=row_sample,
                                     base=os.path.join(out, 'vrhyme')),
        log=[os.path.join(tmp, 'output', 'logs', 'vrhyme_summary.log')])
    run_script(SUMMARY, stub)
    return (read_tsv(stub.output.summary), read_tsv(stub.output.membership),
            open(stub.log[0]).read())


def test_summary():
    print('make_vrhyme_summary:')
    with tempfile.TemporaryDirectory(prefix='magmaker_vsum_') as tmp:
        base = os.path.join(tmp, 'output', 'virus', 'vrhyme', 'megahit', 'S1',
                            'vrhyme_results')
        # Two bins: one holding three contigs, one holding a single contig.
        write(os.path.join(base, 'vRhyme_best_bins.1.membership.tsv'),
              'scaffold\tbin\nk141_1\t1\nk141_2\t1\nk141_3\t1\nk141_9\t2\n')
        write(os.path.join(base, 'vRhyme_best_bins.1.summary.tsv'),
              'bin\tmembers\tproteins\tredundancy\n1\t3\t40\t1\n2\t1\t12\t0\n')
        # A circular scaffold, pulled out before binning.
        write(os.path.join(base, 'vRhyme_circular.tsv'),
              'scaffold\tDTR\nk141_7\tyes\n')
        # vRhyme writes one of these per iteration; they must be ignored.
        write(os.path.join(base, 'vRhyme_alternate_bins',
                           'vRhyme_best_bins.7.membership.tsv'),
              'scaffold\tbin\nk141_1\t99\n')

        rows, members, log = summary(tmp, [('megahit', 'S1')], {'S1': 'S1'})
        check('one row per vMAG, the circular scaffold included',
              len(rows) == 3)
        multi = [r for r in rows if r['Multi_Contig'] == 'yes']
        check('the multi-contig bin is marked, with its contig count',
              len(multi) == 1 and multi[0]['N_Contigs'] == '3')
        check("vRhyme's own redundancy count is carried",
              multi[0]['Redundant_Proteins'] == '1')
        check('the single-contig bin is not marked multi',
              sum(1 for r in rows
                  if r['Multi_Contig'] == 'no' and r['Circular'] == 'no') == 1)
        check('the circular scaffold is kept as its own vMAG',
              sum(1 for r in rows if r['Circular'] == 'yes') == 1)
        check('an alternate iteration is ignored, not reported as the run',
              not any(r['vRhyme_Bin'] == '99' for r in rows))
        check('every contig has a membership row, the circular one included',
              len(members) == 5
              and {m['Seq_Name'] for m in members}
              == {'k141_1', 'k141_2', 'k141_3', 'k141_9', 'k141_7'})
        check('the log reports multi-contig bins with their denominators',
              '1 of 3 vMAGs hold more than one contig' in log
              and '3 of 5 binned contigs are in a multi-contig vMAG' in log)
        check('the vMAG ids carry the sample and assembler',
              all(r['vMAG_ID'].startswith('vMAG_S1_megahit') for r in rows))

    # An assembly whose geNomad calls were all excluded bins nothing; that is
    # an outcome, not an error.
    with tempfile.TemporaryDirectory(prefix='magmaker_vsum2_') as tmp:
        os.makedirs(os.path.join(tmp, 'output', 'virus', 'vrhyme', 'megahit',
                                 'S1', 'vrhyme_results'), exist_ok=True)
        rows, members, log = summary(tmp, [('megahit', 'S1')], {'S1': 'S1'})
        check('an assembly with no bins gives empty tables and says so',
              rows == [] and members == []
              and '0 vMAGs from 0 of 1 assemblies' in log)


if __name__ == '__main__':
    test_prep()
    test_summary()
    print('test_vrhyme: {0} checks passed'.format(len(CHECKS)))
