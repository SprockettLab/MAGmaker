#!/usr/bin/env python3
"""Tests resources/scripts/polish_changes.py with known edits.

Run from the repository root:

    python3 resources/test/test_polish_changes.py
"""

import os
import random
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO_ROOT, 'resources', 'scripts'))
import polish_changes as pc  # noqa: E402

FAILURES = []


def check(name, condition, detail=''):
    print(("PASS  " if condition else "FAIL  ") + name)
    if not condition:
        FAILURES.append(name)
        if detail:
            print("      " + str(detail))


def fasta(d, name, seqs):
    path = os.path.join(d, name)
    with open(path, 'w') as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n{v}\n")
    return path


def main():
    rng = random.Random(3)
    a = "".join(rng.choice("ACGT") for _ in range(20000))   # synthetic contig, 20 kb
    b = "".join(rng.choice("ACGT") for _ in range(5000))
    sub = lambda s, i: s[:i] + ("A" if s[i] != "A" else "C") + s[i + 1:]
    # Polypolish: 3 substitutions in a; pypolca: 1 deletion and 1 insertion in a; b untouched.
    a1 = sub(sub(sub(a, 100), 5000), 15000)
    a2 = a1[:8000] + a1[8001:]
    a2 = a2[:12000] + "G" + a2[12000:]
    with tempfile.TemporaryDirectory() as d:
        draft = fasta(d, 'd.fa', {'c1': a, 'c2': b})
        pp = fasta(d, 'p.fa', {'c1': a1, 'c2': b})
        final = fasta(d, 'f.fa', {'c1': a2, 'c2': b})
        out = os.path.join(d, 'o.tsv')
        pc.main(['--draft', draft, '--polypolish', pp, '--final', final, '--out', out])
        rows = {l.split('\t')[0]: l.rstrip('\n').split('\t')
                for l in open(out) if not l.startswith('#') and not l.startswith('contig\t')}
        check("Polypolish substitutions counted", rows['c1'][2] == '3', rows['c1'])
        check("pypolca indels counted", rows['c1'][3] == '2', rows['c1'])
        check("total draft to final", rows['c1'][4] == '5', rows['c1'])
        check("an untouched contig has 0 changes", rows['c2'][2:5] == ['0', '0', '0'], rows['c2'])
        check("changes per kb", rows['c1'][5] == f"{5 / (len(a2) / 1000):.3f}", rows['c1'])
        check("header states every contig was measured", 'all 2 contigs measured' in open(out).readline())

        lost = fasta(d, 'lost.fa', {'c1': a2})
        try:
            pc.compare(draft, pp, lost)
            check("a lost contig stops", False)
        except pc.PolishError as e:
            check("a lost contig stops", 'c2' in str(e))
        renamed = fasta(d, 'ren.fa', {'c1': a1, 'c2_polished': b})
        try:
            pc.compare(draft, renamed, final)
            check("a renamed contig stops", False)
        except pc.PolishError as e:
            check("a renamed contig stops", 'c2_polished' in str(e))

        big = "".join(rng.choice("ACGT") for _ in range(20000))
        check("a contig changed past the bound is reported as >k",
              pc.changes(a, big, 0.05).startswith('>'))

    print("\n%d failure(s)" % len(FAILURES))
    return 1 if FAILURES else 0


if __name__ == '__main__':
    sys.exit(main())
