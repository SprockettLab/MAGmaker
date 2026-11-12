#!/usr/bin/env python3
"""Tests resources/scripts/control_strand.py: lambda contigs removed at 50%.

Hand-written PAF records (contig as query), so the expected answer is known;
then, if minimap2 is on PATH, a real alignment of contigs built from a random
"lambda" and random sequence.

Run from the repository root:

    python3 resources/test/test_control_strand.py
"""

import os
import random
import shutil
import subprocess
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO_ROOT, 'resources', 'scripts'))
import control_strand as cs  # noqa: E402

FAILURES = []


def check(name, condition, detail=''):
    print(("PASS  " if condition else "FAIL  ") + name)
    if not condition:
        FAILURES.append(name)
        if detail:
            print("      " + str(detail))


def paf(q, qlen, s, e):
    return "\t".join(map(str, [q, qlen, s, e, "+", "lambda", 48502, 0, e - s, e - s, e - s, 60])) + "\n"


def names(path):
    return [l[1:].split()[0] for l in open(path) if l.startswith(">")]


def main():
    with tempfile.TemporaryDirectory() as d:
        contigs = os.path.join(d, 'c.fa')
        with open(contigs, 'w') as fh:
            for n in ("lam_whole", "lam_50", "lam_49", "overlap", "clean"):
                fh.write(f">{n}\nACGT\n")
        p = os.path.join(d, 'x.paf')
        with open(p, 'w') as fh:
            fh.write(paf("lam_whole", 1000, 0, 1000))
            fh.write(paf("lam_50", 1000, 0, 500))
            fh.write(paf("lam_49", 1000, 0, 490))
            # two overlapping alignments of 400: union 600 (60%), not 800
            fh.write(paf("overlap", 1000, 0, 400) + paf("overlap", 1000, 200, 600))
        out = [os.path.join(d, x) for x in ("k.fa", "r.fa", "t.tsv")]
        n_in, n_kept, n_removed = cs.split(contigs, p, 0.5, *out)
        check("exactly 50% is removed, 49% kept", "lam_50" in names(out[1]) and "lam_49" in names(out[0]))
        check("removed set", sorted(names(out[1])) == ["lam_50", "lam_whole", "overlap"])
        check("kept set", sorted(names(out[0])) == ["clean", "lam_49"])
        check("counts add up", (n_in, n_kept, n_removed) == (5, 2, 3))
        rows = {l.split('\t')[0]: l.split('\t') for l in open(out[2]) if not l.startswith(('#', 'contig\t'))}
        check("overlapping alignments counted once", rows["overlap"][2] == "600", rows["overlap"])

    if not shutil.which("minimap2"):
        print("SKIP  real minimap2 run (minimap2 not on PATH)")
    else:
        rng = random.Random(5)
        seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
        lam = seq(48502)
        with tempfile.TemporaryDirectory() as d:
            ref = os.path.join(d, 'lambda.fa')
            open(ref, 'w').write(">NC_001416.1\n" + lam + "\n")
            contigs = os.path.join(d, 'c.fa')
            open(contigs, 'w').write(">lambda_contig\n" + lam[1000:30000] + "\n"
                                     ">bacterial\n" + seq(30000) + "\n"
                                     ">mostly_bacterial\n" + seq(20000) + lam[5000:9000] + "\n")
            p = os.path.join(d, 'x.paf')
            with open(p, 'w') as fh:
                subprocess.run(["minimap2", "-c", ref, contigs], stdout=fh, stderr=subprocess.DEVNULL, check=True)
            out = [os.path.join(d, x) for x in ("k.fa", "r.fa", "t.tsv")]
            cs.split(contigs, p, 0.5, *out)
            check("real minimap2: the lambda contig is removed", names(out[1]) == ["lambda_contig"])
            check("real minimap2: a contig with a short lambda piece is kept",
                  sorted(names(out[0])) == ["bacterial", "mostly_bacterial"])

    print("\n%d failure(s)" % len(FAILURES))
    return 1 if FAILURES else 0


if __name__ == '__main__':
    sys.exit(main())
