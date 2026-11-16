#!/usr/bin/env python3
"""Tests resources/scripts/host_genome_counts.py.

Part 1 builds SAM records by hand, so every expected count is known exactly.
Part 2 runs bowtie2 for real, as host_filter does, on a two-genome reference
built here from random sequence: reads drawn from genome B must be credited
to B, reads from A to A, and the removed total must match bowtie2's summary.
Part 2 is skipped (and says so) when bowtie2 or samtools is not on PATH.

Run from the repository root:

    python3 resources/test/test_host_genome_counts.py
"""

import io
import os
import random
import shutil
import subprocess
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
SCRIPT = os.path.join(REPO_ROOT, 'resources', 'scripts', 'host_genome_counts.py')
sys.path.insert(0, os.path.dirname(SCRIPT))
import host_genome_counts as h  # noqa: E402

FAILURES = []


def check(name, condition):
    print(("PASS  " if condition else "FAIL  ") + name)
    if not condition:
        FAILURES.append(name)


def raises(fn):
    try:
        fn()
    except ValueError:
        return True
    return False


HEADER = "@HD\tVN:1.0\n@SQ\tSN:a1\tLN:1000\n@SQ\tSN:b1\tLN:1000\n@SQ\tSN:bmt\tLN:100\n"
MAP = {"a1": "A", "b1": "B", "bmt": "B_MT"}
MEMBERS = ["A", "B", "B_MT"]


def rec(q, flag, rname="*", mapq=42):
    return "\t".join([q, str(flag), rname, "1", str(mapq), "*", "=", "1",
                      "0", "*", "*"]) + "\n"


def pair(q, rname1, rname2=None, proper=True, mapq1=42, mapq2=42,
         unmapped2=False):
    """Mate 1 then mate 2, as bowtie2 writes them."""
    f1 = 0x1 | 0x40 | (0x2 if proper else 0)
    f2 = 0x1 | 0x80 | (0x2 if proper else 0)
    if unmapped2:
        f2 |= 0x4
        f1 |= 0x8
    return (rec(q, f1, rname1, mapq1)
            + rec(q, f2, rname2 or rname1, 0 if unmapped2 else mapq2))


def run_count(sam, layout="PE", min_mapq=10):
    return h.count(io.StringIO(HEADER + sam), layout, MAP, MEMBERS, min_mapq)


def part1():
    sam = (pair("p1", "a1") + pair("p2", "a1", mapq1=1)
           + pair("p3", "b1") + pair("p4", "bmt", mapq2=5)
           + pair("p5", "a1", "b1", proper=False)       # discordant: not removed
           + pair("p6", "b1", proper=False, unmapped2=True))
    tally, sq = run_count(sam)
    check("concordant pairs credited to mate 1's member",
          [tally[m]["removed"] for m in MEMBERS] == [2, 1, 1])
    check("confident uses the lower MAPQ of the two mates",
          [tally[m]["confident"] for m in MEMBERS] == [1, 1, 0])
    check("aligned mates of non-concordant pairs counted, not removed",
          [tally[m]["mates"] for m in MEMBERS] == [1, 2, 0])
    check("@SQ names collected", sq == {"a1", "b1", "bmt"})

    check("mates out of order fail",
          raises(lambda: run_count(rec("x", 0x1 | 0x80 | 0x2, "a1")
                                   + rec("x", 0x1 | 0x40 | 0x2, "a1"))))
    check("mates of two different reads fail",
          raises(lambda: run_count(rec("x", 0x1 | 0x40 | 0x2, "a1")
                                   + rec("y", 0x1 | 0x80 | 0x2, "a1"))))
    check("a lone last record fails",
          raises(lambda: run_count(rec("x", 0x1 | 0x40 | 0x2, "a1"))))
    check("a secondary record fails",
          raises(lambda: run_count(rec("x", 0x1 | 0x40 | 0x100, "a1"))))
    check("a sequence missing from the map fails",
          raises(lambda: run_count(pair("x", "zz"))))

    tally, _ = run_count(rec("s1", 0, "a1") + rec("s2", 16, "b1", 0)
                         + rec("s3", 0, "bmt", 30), layout="SE")
    check("single-end: aligned reads removed per member",
          [tally[m]["removed"] for m in MEMBERS] == [1, 1, 1])
    check("single-end: confident by MAPQ",
          [tally[m]["confident"] for m in MEMBERS] == [1, 0, 1])
    check("single-end: a paired record fails",
          raises(lambda: run_count(pair("x", "a1"), layout="SE")))

    with tempfile.TemporaryDirectory() as d:
        log = os.path.join(d, "pe.log")
        with open(log, "w") as fh:
            fh.write("1000 reads; of these:\n"
                     "  1000 (100.00%) were paired; of these:\n"
                     "    100 (10.00%) aligned concordantly 0 times\n"
                     "    850 (85.00%) aligned concordantly exactly 1 time\n"
                     "    50 (5.00%) aligned concordantly >1 times\n"
                     "    ----\n"
                     "    100 pairs aligned concordantly 0 times; of these:\n"
                     "      7 (7.00%) aligned discordantly 1 time\n"
                     "90.00% overall alignment rate\n")
        check("bowtie2 PE summary parsed",
              h.bowtie2_summary(log, "PE") == (1000, 900))
        se = os.path.join(d, "se.log")
        with open(se, "w") as fh:
            fh.write("500 reads; of these:\n"
                     "  500 (100.00%) were unpaired; of these:\n"
                     "    300 (60.00%) aligned 0 times\n"
                     "    150 (30.00%) aligned exactly 1 time\n"
                     "    50 (10.00%) aligned >1 times\n"
                     "40.00% overall alignment rate\n")
        check("bowtie2 SE summary parsed",
              h.bowtie2_summary(se, "SE") == (500, 200))
        check("a log with no summary fails",
              raises(lambda: h.bowtie2_summary(se, "PE")))

        # contig-map from a members.tsv, mitochondrion credited to its genome
        fa, fb, fm = (os.path.join(d, n) for n in ("a.fna", "b.fna", "m.fna"))
        open(fa, "w").write(">a1 desc\nACGT\n>a2\nACGT\n")
        open(fb, "w").write(">b1\nACGT\n")
        open(fm, "w").write(">bmt\nACGT\n")
        members = os.path.join(d, "members.tsv")
        with open(members, "w") as fh:
            fh.write("accession\tclean_genome\tclean_md5\trun_record\tmember_of\n"
                     f"GCA_A\t{fa}\tx\tr\t-\nGCA_B\t{fb}\tx\tr\t-\n"
                     f"NC_MT\t{fm}\tx\tr\tGCA_B\n")
        rows = h.contig_map(members=members)
        check("contig map from members.tsv",
              rows == [("a1", "GCA_A", "GCA_A"), ("a2", "GCA_A", "GCA_A"),
                       ("b1", "GCA_B", "GCA_B"), ("bmt", "NC_MT", "GCA_B")])
        open(fm, "w").write(">b1\nACGT\n")
        check("a sequence name in two members fails",
              raises(lambda: h.contig_map(members=members)))
        os.remove(fm)
        check("a missing cleaned FASTA fails",
              raises(lambda: h.contig_map(members=members)))
        check("contig map for a single host FASTA",
              h.contig_map(reference=fa, single_name="Host")
              == [("a1", "Host", "Host"), ("a2", "Host", "Host")])


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def part2():
    if not (shutil.which("bowtie2") and shutil.which("samtools")):
        print("SKIP  part 2: bowtie2 and samtools are not both on PATH")
        return
    rng = random.Random(17)
    with tempfile.TemporaryDirectory() as d:
        ga, gb = rand_seq(rng, 20000), rand_seq(rng, 8000)
        mt = rand_seq(rng, 3000)
        fa, fb, fm = (os.path.join(d, n) for n in ("A.fna", "B.fna", "MT.fna"))
        open(fa, "w").write(f">chrA1\n{ga[:12000]}\n>chrA2\n{ga[12000:]}\n")
        open(fb, "w").write(f">chrB1\n{gb}\n")
        open(fm, "w").write(f">mtB\n{mt}\n")
        ref = os.path.join(d, "host.fna")
        with open(ref, "w") as out:
            for f in (fa, fb, fm):
                out.write(open(f).read())
        members = os.path.join(d, "host.fna.members.tsv")
        with open(members, "w") as fh:
            fh.write("accession\tclean_genome\tclean_md5\trun_record\tmember_of\n"
                     f"GCA_A\t{fa}\tx\tr\t-\nGCA_B\t{fb}\tx\tr\t-\n"
                     f"NC_MT\t{fm}\tx\tr\tGCA_B\n")
        subprocess.run(["bowtie2-build", "-q", ref, os.path.join(d, "host")],
                       check=True)

        # Known truth: 300 pairs from A, 120 from B, 40 from B's mito,
        # 200 from sequence in neither (all must stay nonhost).
        truth = {"GCA_A": 300, "GCA_B": 120, "NC_MT": 40}
        other = rand_seq(rng, 20000)
        r1, r2 = [], []
        n = 0

        def add(src, k):
            nonlocal n
            for _ in range(k):
                start = rng.randrange(0, len(src) - 400)
                frag = src[start:start + 300]
                n += 1
                r1.append(f"@r{n}\n{frag[:100]}\n+\n{'I' * 100}\n")
                r2.append(f"@r{n}\n{revcomp(frag)[:100]}\n+\n{'I' * 100}\n")
        add(ga[:11500], 150)
        add(ga[12000:19500], 150)
        add(gb, truth["GCA_B"])
        add(mt, truth["NC_MT"])
        add(other, 200)
        q1, q2 = os.path.join(d, "R1.fq"), os.path.join(d, "R2.fq")
        open(q1, "w").write("".join(r1))
        open(q2, "w").write("".join(r2))

        # host_filter's own bowtie2 call
        log = os.path.join(d, "S.log")
        bam = os.path.join(d, "S.bam")
        with open(log, "w") as lg, open(bam, "wb") as out:
            bt = subprocess.Popen(
                ["bowtie2", "-p", "2", "-x", os.path.join(d, "host"),
                 "-1", q1, "-2", q2, "--un-conc-gz",
                 os.path.join(d, "nonhost.R%.fastq.gz"), "--no-unal"],
                stdout=subprocess.PIPE, stderr=lg)
            subprocess.run(["samtools", "view", "-bS", "-"], stdin=bt.stdout,
                           stdout=out, check=True)
            bt.wait()

        cmap = os.path.join(d, "map.tsv")
        subprocess.run([sys.executable, SCRIPT, "contig-map", "--members",
                        members, "--out", cmap], check=True)
        res = os.path.join(d, "S.tsv")
        view = subprocess.Popen(["samtools", "view", "-h", bam],
                                stdout=subprocess.PIPE)
        p = subprocess.run([sys.executable, SCRIPT, "count", "--sample", "S",
                            "--layout", "PE", "--contig-map", cmap,
                            "--bowtie2-log", log, "--min-mapq", "10",
                            "--out", res], stdin=view.stdout,
                           capture_output=True, text=True)
        view.wait()
        check("real bowtie2 run: count exits 0 (its own check against the "
              "bowtie2 summary passed)", p.returncode == 0)
        if p.returncode != 0:
            print(p.stderr)
            return
        rows = [l.rstrip("\n").split("\t") for l in open(res)]
        col = {c: i for i, c in enumerate(rows[0])}
        got = {r[col["Member"]]: int(r[col["Removed"]]) for r in rows[1:]}
        print(f"      truth {truth}  counted {got}")
        check("real bowtie2 run: every pair credited to its true genome",
              got == truth)
        check("real bowtie2 run: Input is all 660 pairs",
              all(r[col["Input"]] == "660" for r in rows[1:]))
        check("real bowtie2 run: mitochondrion rolls up to its genome",
              [r[col["Genome"]] for r in rows[1:]
               if r[col["Member"]] == "NC_MT"] == ["GCA_B"])

        # The check can fail: a log claiming one more removed pair stops it.
        bad = os.path.join(d, "bad.log")
        text = open(log).read()
        import re
        m = re.search(r"(\d+) \(([\d.]+)%\) aligned concordantly exactly 1 time",
                      text)
        open(bad, "w").write(text.replace(
            m.group(0), m.group(0).replace(m.group(1), str(int(m.group(1)) + 1), 1)))
        view = subprocess.Popen(["samtools", "view", "-h", bam],
                                stdout=subprocess.PIPE)
        p = subprocess.run([sys.executable, SCRIPT, "count", "--sample", "S",
                            "--layout", "PE", "--contig-map", cmap,
                            "--bowtie2-log", bad, "--min-mapq", "10",
                            "--out", res + ".bad"], stdin=view.stdout,
                           capture_output=True, text=True)
        view.wait()
        check("real bowtie2 run: a total that disagrees with bowtie2 fails",
              p.returncode != 0 and "bowtie2's summary" in p.stderr)

        merged, params = os.path.join(d, "m.tsv"), os.path.join(d, "m_params.txt")
        subprocess.run([sys.executable, SCRIPT, "merge", "--out", merged,
                        "--params", params, "--host", "H", "--reference", ref,
                        "--members", members, "--min-mapq", "10", res],
                       check=True)
        check("merge: header plus one row per member",
              len(open(merged).read().splitlines()) == 4)
        check("merge: parameters file says nothing was sampled",
              "sampling\tnone" in open(params).read())


def main():
    part1()
    part2()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILED")
        sys.exit(1)
    print("all passed")


if __name__ == "__main__":
    main()
