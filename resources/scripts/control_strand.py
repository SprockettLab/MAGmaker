#!/usr/bin/env python3
"""Remove DNA Control Strand (lambda phage) contigs from a long-read assembly.

The ligation kit adds lambda DNA (NC_001416.1) as a control strand; it
assembles into contigs. A contig is removed when the union of its bases that
align to lambda (minimap2 PAF, contig as query) is at least --min-fraction
(0.5) of its length; overlapping alignments count once.

Writes the kept contigs, the removed contigs, and a table of the removed ones
(contig, length, aligned bases, fraction), and checks kept + removed = in.
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ont_host_filter import union_length  # noqa: E402


def aligned_bases(paf_path):
    """contig -> (length, union of aligned contig bases) from a PAF with contigs as queries."""
    spans, lengths = defaultdict(list), {}
    with open(paf_path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 12:
                continue
            lengths[f[0]] = int(f[1])
            spans[f[0]].append((int(f[2]), int(f[3])))
    return {c: (lengths[c], union_length(s)) for c, s in spans.items()}


def split(contigs_path, paf_path, min_fraction, kept_path, removed_path, table_path):
    hits = aligned_bases(paf_path)
    drop = {c: v for c, v in hits.items() if v[1] >= min_fraction * v[0]}
    n_in = n_kept = n_removed = 0
    with open(contigs_path) as src, open(kept_path, "w") as kept, open(removed_path, "w") as removed:
        out = kept
        for line in src:
            if line.startswith(">"):
                name = line[1:].split()[0]
                n_in += 1
                if name in drop:
                    out, n_removed = removed, n_removed + 1
                else:
                    out, n_kept = kept, n_kept + 1
            out.write(line)
    if n_kept + n_removed != n_in or n_removed != len(drop):
        raise SystemExit(f"CHECK FAILED: {n_kept} kept + {n_removed} removed != {n_in} in "
                         f"(or {len(drop)} to remove)")
    with open(table_path, "w") as t:
        t.write(f"# contigs in {n_in}, removed {n_removed}, kept {n_kept}; removed when aligned "
                f"bases (union) >= {min_fraction} of length to the control strand\n")
        t.write("contig\tlength\taligned_bases\tfraction\n")
        for c, (L, a) in sorted(drop.items()):
            t.write(f"{c}\t{L}\t{a}\t{a / L:.4f}\n")
    return n_in, n_kept, n_removed


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    for k in ("contigs", "paf", "kept", "removed", "table"):
        p.add_argument(f"--{k}", required=True)
    p.add_argument("--min-fraction", type=float, default=0.5)
    a = p.parse_args(argv)
    split(a.contigs, a.paf, a.min_fraction, a.kept, a.removed, a.table)


if __name__ == "__main__":
    main()
