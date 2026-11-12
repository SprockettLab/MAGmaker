#!/usr/bin/env python3
"""Tests the parsers in resources/scripts/count_rrna_trna.py on hand-written
barrnap GFF and tRNAscan-SE output (formats of barrnap 0.9+ and tRNAscan-SE 2).

Run from the repository root:

    python3 resources/test/test_count_rrna_trna.py
"""

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO_ROOT, 'resources', 'scripts'))
import count_rrna_trna as c  # noqa: E402

FAILURES = []


def check(name, condition, detail=''):
    print(("PASS  " if condition else "FAIL  ") + name)
    if not condition:
        FAILURES.append(name)
        if detail:
            print("      " + str(detail))


BAC = """##gff-version 3
c1\tbarrnap:0.9\trRNA\t100\t1600\t0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA
c1\tbarrnap:0.9\trRNA\t1800\t4700\t0\t+\t.\tName=23S_rRNA;product=23S ribosomal RNA
c1\tbarrnap:0.9\trRNA\t4800\t4910\t0\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA
c2\tbarrnap:0.9\trRNA\t1\t700\t0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA (partial);note=aligned only 45 percent of the 16S ribosomal RNA
"""
ARC = """##gff-version 3
c1\tbarrnap:0.9\trRNA\t4800\t4910\t0\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA
c3\tbarrnap:0.9\trRNA\t10\t120\t0\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA
"""
TRNA = """Sequence\t\ttRNA\tBounds\ttRNA\tAnti\tIntron Bounds\tInf\t\t
Name    \ttRNA #\tBegin\tEnd\tType\tCodon\tBegin\tEnd\tScore\tNote
--------\t------\t-----\t------\t----\t-----\t-----\t----\t------\t------
c1\t1\t100\t172\tAla\tTGC\t0\t0\t70.1\t
c1\t2\t200\t272\tAla\tGGC\t0\t0\t69.0\t
c1\t3\t300\t385\tLeu\tCAG\t0\t0\t60.2\t
c2\t1\t50\t122\tSeC\tTCA\t0\t0\t55.0\t
c2\t2\t400\t472\tUndet\tNNN\t0\t0\t20.1\t
c2\t3\t600\t672\tPseudo\tAAA\t0\t0\t15.3\t
"""


PARTIAL_5S = ("c4\tbarrnap:0.9\trRNA\t1\t70\t0\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA (partial);"
              "note=aligned only 62 percent of the 5S ribosomal RNA\n")


def main():
    b = c.count_rrna(BAC, 0.5)
    check("at 0.5: complete genes counted, a 45% 16S not counted",
          b == {'5S_rRNA': 1, '16S_rRNA': 1, '23S_rRNA': 1}, b)
    a = c.count_rrna(ARC, 0.5)
    check("archaeal model counts", a == {'5S_rRNA': 2, '16S_rRNA': 0, '23S_rRNA': 0}, a)
    check("the larger count of the two models per gene",
          c.best_of(b, a) == {'5S_rRNA': 2, '16S_rRNA': 1, '23S_rRNA': 1})
    check("at 0.33 the 45% 16S counts (same rule for every gene)",
          c.count_rrna(BAC, 0.33)['16S_rRNA'] == 2)
    check("a 62% 5S counts at 0.5, not at 0.8",
          c.count_rrna(PARTIAL_5S, 0.5)['5S_rRNA'] == 1
          and c.count_rrna(PARTIAL_5S, 0.8)['5S_rRNA'] == 0)
    check("fraction from the barrnap note; no note means complete",
          c.aligned_fraction("Name=16S_rRNA;note=aligned only 45 percent of the 16S") == 0.45
          and c.aligned_fraction("Name=16S_rRNA;product=16S ribosomal RNA") == 1.0)
    n = c.trna_amino_acids(TRNA)
    check("distinct standard amino acids (Ala twice counts once; SeC, Undet, Pseudo not counted)", n == 2, n)
    print("\n%d failure(s)" % len(FAILURES))
    return 1 if FAILURES else 0


if __name__ == '__main__':
    sys.exit(main())
