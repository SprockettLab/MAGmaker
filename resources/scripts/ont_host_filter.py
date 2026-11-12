#!/usr/bin/env python3
"""Long-read host calls from minimap2 SAM output.

minimap2 -a -x map-ont --secondary=no writes every read: unmapped, or one
primary alignment plus any supplementary ones, a read's records together.
A read is host when the union of the read positions its alignments cover is
at least --threshold of its length (overlapping alignments counted once).
Reads with some host alignment below the threshold are kept, and counted.

Reads SAM on stdin. Writes:
  stdout          SAM header plus every aligned record (for samtools: host BAM)
  --host-ids      gzipped list of host read IDs
  --counts        counts and every parameter, key<TAB>value
All reads are used; nothing is sampled.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
QUERY_OPS = set("MI=X")      # consume read bases inside the alignment
CLIP_OPS = set("SH")         # read bases outside it


def query_interval(flag: int, cigar: str) -> tuple[int, int, int]:
    """in : SAM flag and CIGAR of an aligned record
       out: (start, end, read length), start/end in the read's own
            orientation, 0-based half-open

    For a reverse-strand record the CIGAR runs along the reverse complement,
    so its interval is mirrored back onto the read. Hard and soft clips both
    count towards the read length (supplementary records are hard-clipped)."""
    ops = [(int(n), op) for n, op in CIGAR.findall(cigar)]
    lead = 0
    for n, op in ops:
        if op in CLIP_OPS:
            lead += n
        else:
            break
    aligned = sum(n for n, op in ops if op in QUERY_OPS)
    length = sum(n for n, op in ops if op in QUERY_OPS or op in CLIP_OPS)
    start, end = lead, lead + aligned
    if flag & 16:
        start, end = length - end, length - start
    return start, end, length


def union_length(intervals: list[tuple[int, int]]) -> int:
    """Bases covered by at least one interval."""
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                total += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s
    return total


def classify(aligned: int, length: int, threshold: float) -> str:
    """'host', 'partial' (some host alignment, below the threshold) or 'clean'."""
    if aligned == 0:
        return "clean"
    return "host" if aligned >= threshold * length else "partial"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--threshold", type=float, required=True)
    p.add_argument("--host-ids", required=True)
    p.add_argument("--counts", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--minimap2-version", required=True)
    a = p.parse_args(argv)
    if not 0 < a.threshold <= 1:
        sys.exit(f"--threshold must be in (0, 1], got {a.threshold}")

    n = {k: 0 for k in ("reads_in", "bases_in", "reads_host", "bases_host",
                        "reads_partial", "bases_partial")}
    out = sys.stdout
    ids = gzip.open(a.host_ids, "wt")

    def finish(name, length, intervals):
        if name is None:
            return
        n["reads_in"] += 1
        n["bases_in"] += length
        kind = classify(union_length(intervals), length, a.threshold)
        if kind == "host":
            n["reads_host"] += 1
            n["bases_host"] += length
            ids.write(name + "\n")
        elif kind == "partial":
            n["reads_partial"] += 1
            n["bases_partial"] += length

    name, length, intervals = None, 0, []
    for line in sys.stdin:
        if line.startswith("@"):
            out.write(line)
            continue
        f = line.split("\t", 10)
        qname, flag, cigar, seq = f[0], int(f[1]), f[5], f[9]
        if qname != name:
            finish(name, length, intervals)
            name, length, intervals = qname, 0, []
        if flag & 4:
            length = len(seq)
            continue
        s, e, L = query_interval(flag, cigar)
        if not flag & 0x900:          # the primary record defines the read length
            length = L
        intervals.append((s, e))
        out.write(line)
    finish(name, length, intervals)
    ids.close()

    kept = n["reads_in"] - n["reads_host"]
    with open(a.counts, "w") as fh:
        for k, v in [("sample", a.sample), ("reference", a.reference),
                     ("aligner", f"minimap2 {a.minimap2_version} -a -x map-ont --secondary=no"),
                     ("host_rule", f"union of aligned read bases >= {a.threshold} of read length"),
                     ("reads_used", "all reads, none sampled"),
                     ("reads_in", n["reads_in"]), ("bases_in", n["bases_in"]),
                     ("reads_host", n["reads_host"]), ("bases_host", n["bases_host"]),
                     ("reads_kept", kept), ("bases_kept", n["bases_in"] - n["bases_host"]),
                     ("reads_kept_with_partial_host_alignment", n["reads_partial"]),
                     ("bases_kept_with_partial_host_alignment", n["bases_partial"])]:
            fh.write(f"{k}\t{v}\n")


if __name__ == "__main__":
    main()
