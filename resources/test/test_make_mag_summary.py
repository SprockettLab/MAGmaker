#!/usr/bin/env python3
"""Runs make_mag_summary.py end to end on a one-bin fixture.

The script reads its inputs from a `snakemake` object, so this builds a
stub with the same attributes and execs the script against a temporary
output tree. No snakemake install or cluster is needed.

Run from the repository root:

    python3 resources/test/test_make_mag_summary.py

Issue #1: the per-bin loop referenced full_classification before
unpacking it from the GTDB-Tk entry. A sample whose first bin had a
GTDB-Tk row raised NameError and no mag_summary.tsv was written. Where the
first bin had no row, the `or` short-circuited and the run went on, but
every later bin was checked against the previous bin's classification,
so the missing-taxonomy warning miscounted. One test covers each.
"""

import os
import sys
import tempfile
import types

import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
SCRIPT = os.path.join(REPO_ROOT, 'resources', 'scripts',
                      'make_mag_summary.py')

MAPPER = 'minimap2'
SAMPLE = 'George'
SELECTED_FASTAS = 'Binette_Fastas'


CLASSIFIED = ('d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;'
              'f__Lachnospiraceae;g__Blautia_A;s__Blautia_A obeum')


def build_fixture(root, bins):
    """
    in : empty directory; dict of bin name -> GTDB-Tk classification, or
         None for a bin GTDB-Tk has no row for
    out: stub snakemake object pointing at an output tree under root
    """
    bins_dir = os.path.join(root, 'selected_bins', MAPPER, SELECTED_FASTAS,
                            SAMPLE)
    os.makedirs(bins_dir)
    for bin_name in bins:
        with open(os.path.join(bins_dir, bin_name + '.fa'), 'w') as fh:
            fh.write('>contig_1\nACGTACGTACGTGGCC\n>contig_2\nGGCCAATT\n')

    gtdbtk_dir = os.path.join(root, 'gtdbtk', MAPPER, SAMPLE)
    os.makedirs(gtdbtk_dir)
    gtdbtk_df = pd.DataFrame(
        [{'user_genome': name, 'classification': tax, 'msa_percent': 90.0,
          'warnings': 'N/A'}
         for name, tax in bins.items() if tax is not None])
    gtdbtk_df.to_csv(os.path.join(gtdbtk_dir, 'gtdbtk.bac120.summary.tsv'),
                     sep='\t', index=False)

    params = types.SimpleNamespace(
        mappers=[MAPPER],
        contig_samples=[SAMPLE],
        bins_base=os.path.join(root, 'selected_bins'),
        gtdbtk_base=os.path.join(root, 'gtdbtk'),
        checkm2_base=os.path.join(root, 'checkm2'),
        gunc_base=os.path.join(root, 'gunc'),
        cmseq_base=os.path.join(root, 'cmseq'),
        assemblers=['megahit'],
        consolidation_tool='binette',
        selected_fastas=SELECTED_FASTAS,
    )
    return types.SimpleNamespace(
        params=params,
        log=[os.path.join(root, 'make_mag_summary.log')],
        output=types.SimpleNamespace(
            summary=os.path.join(root, 'mag_summary.tsv')),
    )


def gtdbtk_warning(log_path):
    """
    in : script log
    out: the GTDB-Tk missing-result warning line, or '' if none
    """
    with open(log_path) as fh:
        lines = [l for l in fh if 'no GTDB-Tk result' in l]
    return lines[0] if lines else ''


def run_script(stub):
    """
    in : stub snakemake object
    out: None; raises whatever the script raises

    The script redirects stdout and stderr to its log, so both are
    restored afterwards whether or not it succeeds.
    """
    saved = (sys.stdout, sys.stderr)
    try:
        with open(SCRIPT) as fh:
            code = compile(fh.read(), SCRIPT, 'exec')
        exec(code, {'snakemake': stub, '__name__': '__main__'})
    finally:
        for stream in (sys.stdout, sys.stderr):
            if stream not in saved:
                stream.close()
        sys.stdout, sys.stderr = saved


def test_classified_first_bin_does_not_crash():
    # The first classified bin is what raised NameError: with no GTDB-Tk
    # row the `or` short-circuits before reaching the unassigned name.
    with tempfile.TemporaryDirectory() as root:
        stub = build_fixture(root, {'bin.1': CLASSIFIED})
        run_script(stub)

        summary_df = pd.read_csv(stub.output.summary, sep='\t',
                                 keep_default_na=False)
        assert len(summary_df) == 1, summary_df
        row = summary_df.iloc[0]
        assert row['GTDB_Classification'] == CLASSIFIED, row
        assert row['Species'] == 'Blautia_A obeum', row
        assert int(row['Num_Contigs']) == 2, row
        assert gtdbtk_warning(stub.log[0]) == '', gtdbtk_warning(stub.log[0])


def test_missing_taxonomy_counts_only_unclassified_bins():
    # Before the fix the check read the PREVIOUS bin's classification, so
    # a classified bin after an unclassified one was also reported as
    # missing taxonomy, and the warning said ALL bins when only one was.
    with tempfile.TemporaryDirectory() as root:
        stub = build_fixture(root, {'bin.1': None, 'bin.2': CLASSIFIED})
        run_script(stub)

        warning = gtdbtk_warning(stub.log[0])
        assert 'for 1 of 2 bins (e.g. bin.1)' in warning, warning


if __name__ == '__main__':
    failed = 0
    for test in (test_classified_first_bin_does_not_crash,
                 test_missing_taxonomy_counts_only_unclassified_bins):
        try:
            test()
            print('ok:  ', test.__name__)
        except Exception as err:
            print('FAIL:', test.__name__, '->', type(err).__name__, err)
            failed += 1
    sys.exit(1 if failed else 0)
