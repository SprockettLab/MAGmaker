#!/usr/bin/env python3
"""Tests the long-read host rule in resources/scripts/ont_host_filter.py.

A read is host when the union of the read positions its primary and
supplementary alignments cover is at least the threshold (0.5) of its length.
These build SAM records by hand, so the expected answer is known exactly.

Run from the repository root:

    python3 resources/test/test_ont_host_filter.py
"""

import gzip
import io
import os
import sys
import tempfile
from contextlib import redirect_stdout

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO_ROOT, 'resources', 'scripts'))
import ont_host_filter as h  # noqa: E402

FAILURES = []


def check(name, condition):
    print(("PASS  " if condition else "FAIL  ") + name)
    if not condition:
        FAILURES.append(name)


def rec(name, flag, cigar, seq="*"):
    return "\t".join([name, str(flag), "chr1" if not flag & 4 else "*", "1", "60",
                      cigar if not flag & 4 else "*", "*", "0", "0", seq, "*"]) + "\n"


def run(sam_records, threshold=0.5):
    with tempfile.TemporaryDirectory() as d:
        ids, counts = os.path.join(d, "ids.gz"), os.path.join(d, "counts.tsv")
        sys.stdin = io.StringIO("@HD\tVN:1.6\n" + "".join(sam_records))
        out = io.StringIO()
        try:
            with redirect_stdout(out):
                h.main(["--threshold", str(threshold), "--host-ids", ids, "--counts", counts,
                        "--sample", "S", "--reference", "ref.fna", "--minimap2-version", "2.28"])
        finally:
            sys.stdin = sys.__stdin__
        host = gzip.open(ids, "rt").read().split()
        c = dict(l.rstrip("\n").split("\t", 1) for l in open(counts))
        return host, c, out.getvalue()


def main():
    # Interval arithmetic.
    check("forward interval after a soft clip",
          h.query_interval(0, "100S400M500S") == (100, 500, 1000))
    check("reverse-strand interval mirrored onto the read",
          h.query_interval(16, "100S400M500S") == (500, 900, 1000))
    check("hard clips count towards the read length",
          h.query_interval(2048, "600H400M") == (600, 1000, 1000))
    check("insertions are read bases, deletions are not",
          h.query_interval(0, "100M10I50D90M") == (0, 200, 200))
    check("overlapping intervals counted once",
          h.union_length([(0, 400), (300, 600)]) == 600)
    check("separate intervals summed",
          h.union_length([(0, 300), (700, 1000)]) == 600)
    check("exactly at the threshold is host", h.classify(500, 1000, 0.5) == "host")
    check("just below the threshold is partial", h.classify(499, 1000, 0.5) == "partial")
    check("no alignment is clean", h.classify(0, 1000, 0.5) == "clean")

    # Whole stream: a 1,000 bp read in two supplementary pieces (60%), one at
    # 40%, one with overlapping pieces (only 45% once overlaps are merged),
    # and an unmapped read.
    records = [
        rec("host60", 0, "300M700S"), rec("host60", 2048, "700H300M"),
        rec("part40", 16, "600S400M"),
        rec("overlap45", 0, "300M700S"), rec("overlap45", 2048, "150H300M550H"),
        rec("unmapped", 4, "*", "A" * 800),
    ]
    host, c, out = run(records)
    check("60% across two pieces is host", host == ["host60"])
    check("reads and bases counted",
          c["reads_in"] == "4" and c["bases_in"] == str(1000 + 1000 + 1000 + 800))
    check("host counts", c["reads_host"] == "1" and c["bases_host"] == "1000")
    check("kept counts", c["reads_kept"] == "3" and c["bases_kept"] == "2800")
    check("partial host reads counted (40% and the merged 45%)",
          c["reads_kept_with_partial_host_alignment"] == "2"
          and c["bases_kept_with_partial_host_alignment"] == "2000")
    check("parameters written into the counts",
          "0.5" in c["host_rule"] and "--secondary=no" in c["aligner"]
          and c["reads_used"].startswith("all reads"))
    check("aligned records passed on for the host BAM, unmapped dropped",
          out.count("\n") == 1 + 5 and "unmapped" not in out)

    real_minimap2()

    print("\n%d failure(s)" % len(FAILURES))
    return 1 if FAILURES else 0


def real_minimap2():
    """Real minimap2 output on reads cut from a random reference, so each
    read's host fraction is known by construction. Skipped without minimap2."""
    import random
    import shutil
    import subprocess
    if not shutil.which("minimap2"):
        print("SKIP  real minimap2 run (minimap2 not on PATH)")
        return
    rng = random.Random(7)
    seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
    ref = seq(100000)
    reads = {
        "whole_host": ref[1000:6000],                                   # 100% host
        "rc_host": rc(ref[20000:25000]),                                # 100%, reverse strand
        "host60": ref[30000:33000] + seq(2000),                         # 60%
        "part40": seq(3000) + ref[40000:42000],                         # 40%
        "split78": ref[50000:52000] + seq(1000) + ref[80000:81500],     # 3500/4500, two alignments
        "clean": seq(5000),                                             # 0%
    }
    with tempfile.TemporaryDirectory() as d:
        r, q = os.path.join(d, "ref.fa"), os.path.join(d, "reads.fq")
        open(r, "w").write(">ref\n" + ref + "\n")
        open(q, "w").write("".join(f"@{k}\n{v}\n+\n{'I' * len(v)}\n" for k, v in reads.items()))
        sam = subprocess.run(["minimap2", "-a", "-x", "map-ont", "--secondary=no", "-t", "2", r, q],
                             capture_output=True, text=True, check=True).stdout
        names = [l.split("\t")[0] for l in sam.splitlines() if not l.startswith("@")]
        grouped = all(names.index(n) + names.count(n) - 1 == len(names) - 1 - names[::-1].index(n)
                      for n in set(names))
        check("real minimap2: a read's records come together", grouped)
        check("real minimap2: split read gets a supplementary alignment",
              any(l.split("\t")[0] == "split78" and int(l.split("\t")[1]) & 2048
                  for l in sam.splitlines() if not l.startswith("@")))
        host, c, _ = run([l + "\n" for l in sam.splitlines() if not l.startswith("@")])
        check("real minimap2: host calls", sorted(host) == ["host60", "rc_host", "split78", "whole_host"])
        check("real minimap2: one partial read (40%)", c["reads_kept_with_partial_host_alignment"] == "1")
        check("real minimap2: read lengths recovered",
              c["bases_in"] == str(sum(len(v) for v in reads.values())))


if __name__ == '__main__':
    sys.exit(main())
