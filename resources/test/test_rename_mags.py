#!/usr/bin/env python3
"""Runs rename_mags.py on two MAGs whose contigs share names.

Assemblers reuse contig names in every sample (MEGAHIT's k141_1), so two MAGs
from different samples can hold contigs with the same name. rename_mags.py
prefixes each contig with its MAG_ID so names are unique across all MAGs,
writes contig_map.tsv, and leaves the source FASTAs unchanged.

The script reads its inputs from a `snakemake` object, so this builds a stub
and execs the script. Run from the repository root:

    python3 resources/test/test_rename_mags.py
"""

import os
import sys
import tempfile
import types

import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
SCRIPT = os.path.join(REPO_ROOT, 'resources', 'scripts', 'rename_mags.py')


def run(root, rows):
    summary = os.path.join(root, 'mag_summary.tsv')
    pd.DataFrame(rows).to_csv(summary, sep='\t', index=False)
    stub = types.SimpleNamespace(
        input=types.SimpleNamespace(summary=summary),
        params=types.SimpleNamespace(renamed_dir=os.path.join(root, 'renamed_mags')),
        log=[os.path.join(root, 'rename.log')])
    with open(SCRIPT) as fh:
        code = compile(fh.read(), SCRIPT, 'exec')
    saved = sys.stdout, sys.stderr
    try:
        exec(code, {'snakemake': stub, '__name__': '__main__'})
    finally:
        sys.stdout.close() if sys.stdout is not saved[0] else None
        sys.stdout, sys.stderr = saved
    return os.path.join(root, 'renamed_mags')


def two_mags(root):
    a = os.path.join(root, 'George_bin.1.fa')
    b = os.path.join(root, 'Ringo_bin.3.fa')
    with open(a, 'w') as fh:
        fh.write('>k141_1 flag=1 multi=5.0\nACGT\n>k141_2\nGGCC\n')
    with open(b, 'w') as fh:
        fh.write('>k141_1 flag=0\nTTAA\n')
    rows = [
        {'MAG_ID': 'MAG_0001', 'New_Name': 'MAG_0001__Blautia_obeum', 'Original_Path': a},
        {'MAG_ID': 'MAG_0002', 'New_Name': 'MAG_0002__NA', 'Original_Path': b},
        {'MAG_ID': 'NONE', 'New_Name': 'NONE', 'Original_Path': ''},
    ]
    return a, b, rows


def test_contig_names_unique_and_sources_untouched():
    with tempfile.TemporaryDirectory() as root:
        a, b, rows = two_mags(root)
        before = open(a).read(), open(b).read()
        out = run(root, rows)
        first = open(os.path.join(out, 'MAG_0001__Blautia_obeum.fa')).read()
        second = open(os.path.join(out, 'MAG_0002__NA.fa')).read()
        assert first == '>MAG_0001__k141_1 flag=1 multi=5.0\nACGT\n>MAG_0001__k141_2\nGGCC\n', first
        assert second == '>MAG_0002__k141_1 flag=0\nTTAA\n', second
        assert (open(a).read(), open(b).read()) == before, 'source FASTAs changed'
        m = pd.read_csv(os.path.join(out, 'contig_map.tsv'), sep='\t')
        assert list(m['contig']) == ['MAG_0001__k141_1', 'MAG_0001__k141_2', 'MAG_0002__k141_1']
        assert list(m['original_contig']) == ['k141_1', 'k141_2', 'k141_1']
        assert sorted(f for f in os.listdir(out) if f.endswith('.fa')) == \
            ['MAG_0001__Blautia_obeum.fa', 'MAG_0002__NA.fa']


def test_duplicate_new_name_refused():
    with tempfile.TemporaryDirectory() as root:
        a, b, rows = two_mags(root)
        rows[1]['New_Name'] = rows[0]['New_Name']
        try:
            run(root, rows)
        except ValueError as e:
            assert 'New_Name must be unique' in str(e)
        else:
            raise AssertionError('two MAGs written to one file were not refused')


if __name__ == '__main__':
    test_contig_names_unique_and_sources_untouched()
    test_duplicate_new_name_refused()
    print('test_rename_mags: 2 passed')
