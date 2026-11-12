#!/usr/bin/env python3
"""Tests resources/scripts/long_read_contigs.py on hand-written Flye and
myloasm output (formats as written by Flye 2.9 and myloasm 0.7).

Run from the repository root:

    python3 resources/test/test_long_read_contigs.py
"""

import os
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO_ROOT, 'resources', 'scripts'))
import long_read_contigs as lc  # noqa: E402

FAILURES = []


def check(name, condition, detail=''):
    print(("PASS  " if condition else "FAIL  ") + name)
    if not condition:
        FAILURES.append(name)
        if detail:
            print("      " + str(detail))


def write(d, name, text):
    path = os.path.join(d, name)
    with open(path, 'w') as fh:
        fh.write(text)
    return path


FLYE_INFO = ("#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
             "contig_1\t12\t35\tY\tN\t1\t*\t1\n"
             "contig_2\t8\t12\tN\tN\t1\t*\t2\n")
FLYE_FASTA = ">contig_1\nACGTACGTAC\nGT\n>contig_2\nGGCCGGCC\n"


def main():
    with tempfile.TemporaryDirectory() as d:
        info, fa = write(d, 'info.txt', FLYE_INFO), write(d, 'a.fa', FLYE_FASTA)
        out_fa, out_tsv = os.path.join(d, 'o.fa'), os.path.join(d, 'o.tsv')
        lc.main(['--sample', 'v1928', '--assembler', 'flye', '--assembly', fa,
                 '--flye-info', info, '--out-fasta', out_fa, '--out-info', out_tsv])
        rows = [l.rstrip('\n').split('\t') for l in open(out_tsv)]
        check("flye: renamed <sample>__<assembler>__<original>",
              open(out_fa).read().startswith(">v1928__flye__contig_1\nACGTACGTACGT\n"))
        check("flye: table header", rows[0] == ['contig', 'original_name', 'length', 'circular', 'depth'])
        check("flye: Y is circular yes, N is no, depth and length kept",
              rows[1] == ['v1928__flye__contig_1', 'contig_1', '12', 'yes', '35']
              and rows[2] == ['v1928__flye__contig_2', 'contig_2', '8', 'no', '12'], rows)

        bad = write(d, 'bad.fa', ">contig_1\nACGT\n>contig_2\nGGCCGGCC\n")
        try:
            list(lc.contigs('s', 'flye', bad, info))
            check("flye: length mismatch stops", False)
        except lc.ContigError as e:
            check("flye: length mismatch stops", 'contig_1' in str(e))
        missing = write(d, 'missing.fa', ">contig_1\nACGTACGTACGT\n")
        try:
            list(lc.contigs('s', 'flye', missing, info))
            check("flye: a contig missing from the FASTA stops", False)
        except lc.ContigError as e:
            check("flye: a contig missing from the FASTA stops", 'contig_2' in str(e))

        mylo = write(d, 'm.fa',
                     ">u1ctg_len-6_circular-yes_depth-10-9-9_duplicated-no mult=1.00\nACGTAC\n"
                     ">u2ctg_len-4_circular-possibly_depth-3-3-3_duplicated-no mult=1.00\nGGCC\n"
                     ">u3ctg_len-5_circular-no_depth-7-7-7_duplicated-no mult=1.00\nTTTTT\n")
        rows = [(n, o, c, dp) for n, o, s, c, dp in lc.contigs('v2980', 'myloasm', mylo)]
        check("myloasm: circular yes / possibly / no read from the header",
              [r[2] for r in rows] == ['yes', 'possibly', 'no'], rows)
        check("myloasm: renamed from the name before the first space",
              rows[0][0] == 'v2980__myloasm__u1ctg_len-6_circular-yes_depth-10-9-9_duplicated-no')
        check("myloasm: depth from the header", rows[0][3] == '10-9-9')
        badm = write(d, 'bm.fa', ">u1ctg_len-9_circular-yes_depth-1_duplicated-no\nACGT\n")
        try:
            list(lc.contigs('s', 'myloasm', badm))
            check("myloasm: length mismatch stops", False)
        except lc.ContigError:
            check("myloasm: length mismatch stops", True)

    print("\n%d failure(s)" % len(FAILURES))
    return 1 if FAILURES else 0


if __name__ == '__main__':
    sys.exit(main())
