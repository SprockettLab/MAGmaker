#!/usr/bin/env python3
"""Per-contig changes made by short-read polishing of a long-read assembly.

Takes the draft, the Polypolish output and the pypolca output. Checks that all
three hold exactly the same contig names (polishing may change sequence, never
add, drop or rename contigs), then counts, per contig, the bases each stage
changed: the edit distance (substitutions, insertions, deletions) between the
stage's input and output, computed exactly with edlib. Contigs where polishing
changed many bases are where it may have mixed strains.

The distance is computed up to a bound (--max-fraction of the contig length,
default 0.05); a contig changed more than that is reported as ">N" rather than
spending quadratic time on it. Every contig is measured; none is skipped.
"""

from __future__ import annotations

import argparse
import sys

import edlib


class PolishError(ValueError):
    pass


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
                if name in seqs:
                    raise PolishError(f"{path}: contig {name} occurs twice")
            elif line:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def changes(a: str, b: str, max_fraction: float) -> str:
    """Edit distance between two versions of one contig, or '>k' past the bound."""
    if a == b:
        return "0"
    k = max(100, int(max_fraction * max(len(a), len(b))))
    d = edlib.align(a, b, mode="NW", task="distance", k=k)["editDistance"]
    return str(d) if d >= 0 else f">{k}"


def compare(draft, polypolish, final, max_fraction=0.05):
    d, p, f = read_fasta(draft), read_fasta(polypolish), read_fasta(final)
    for label, other in (("Polypolish", p), ("pypolca", f)):
        if set(other) != set(d):
            lost, added = sorted(set(d) - set(other)), sorted(set(other) - set(d))
            raise PolishError(f"{label} output does not hold the draft's contigs: "
                              f"lost {lost[:5]}, added {added[:5]}")
    rows = []
    for name in d:
        c1 = changes(d[name], p[name], max_fraction)
        c2 = changes(p[name], f[name], max_fraction)
        total = changes(d[name], f[name], max_fraction)
        per_kb = (f"{int(total) / (len(f[name]) / 1000):.3f}" if total.isdigit() and f[name] else "NA")
        rows.append((name, len(f[name]), c1, c2, total, per_kb))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--draft", required=True)
    ap.add_argument("--polypolish", required=True)
    ap.add_argument("--final", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-fraction", type=float, default=0.05)
    a = ap.parse_args(argv)
    rows = compare(a.draft, a.polypolish, a.final, a.max_fraction)
    with open(a.out, "w") as fh:
        fh.write(f"# edit distance per stage, exact up to {a.max_fraction} of the contig length; "
                 f"all {len(rows)} contigs measured\n")
        fh.write("contig\tlength\tpolypolish_changes\tpypolca_changes\ttotal_changes\tchanges_per_kb\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
