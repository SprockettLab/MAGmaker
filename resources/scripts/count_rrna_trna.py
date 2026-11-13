#!/usr/bin/env python3
"""rRNA and tRNA genes per MAG, for the full MIMAG high-quality standard.

For every *.fa in a bin directory:
  rRNA  barrnap with the bacterial and with the archaeal model. A gene counts
        when at least --min-fraction of its expected length aligned, one rule
        for 5S, 16S and 23S. barrnap labels a gene partial below 80% of the
        expected length and writes "aligned only N percent" for it; genes
        without that note are complete (100%). barrnap itself drops genes below
        its --reject fraction (default 0.25), which is lowered to --min-fraction
        when that is smaller. Per gene, the larger of the two models' counts.
  tRNA  tRNAscan-SE, general model (-G): the number of the 20 standard amino
        acids with at least one tRNA.
Every MAG in the directory is measured. Writes one row per MAG; an empty
directory gives a header-only table.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import subprocess
import sys
import tempfile

RRNA = ("5S_rRNA", "16S_rRNA", "23S_rRNA")
AMINO_ACIDS = {"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
               "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"}


BARRNAP_REJECT = 0.25   # barrnap's own default --reject
ALIGNED = re.compile(r"aligned only (\d+) percent")


def aligned_fraction(attributes: str) -> float:
    """Fraction of the expected gene length barrnap aligned: from its
    'aligned only N percent' note, or 1.0 when there is none (not partial)."""
    m = ALIGNED.search(attributes)
    return int(m.group(1)) / 100 if m else 1.0


def count_rrna(gff_text: str, min_fraction: float) -> dict[str, int]:
    """5S, 16S and 23S genes in one barrnap GFF with at least min_fraction of
    their expected length aligned."""
    n = {g: 0 for g in RRNA}
    for line in gff_text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        name = f[8].split(";")[0].replace("Name=", "")
        if name in n and aligned_fraction(f[8]) >= min_fraction:
            n[name] += 1
    return n


def best_of(bacterial: dict, archaeal: dict) -> dict:
    return {g: max(bacterial.get(g, 0), archaeal.get(g, 0)) for g in RRNA}


def trna_amino_acids(trnascan_text: str) -> int:
    """Distinct standard amino acids with a tRNA in a tRNAscan-SE table
    (three header lines, then: sequence, number, begin, end, type, ...)."""
    seen = set()
    for line in trnascan_text.splitlines()[3:]:
        f = line.split()
        if len(f) > 4 and f[4] in AMINO_ACIDS:
            seen.add(f[4])
    return len(seen)


def run(cmd):
    """Run a tool and return its stdout; on failure, exit with the tool's own
    stderr so the reason reaches the rule log."""
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"{' '.join(cmd)}\nexited with status {p.returncode}; its stderr:\n{p.stderr}")
    return p.stdout


def check_barrnap_version():
    """The partial-gene rule reads barrnap 0.9's 'aligned only N percent' note
    and its --reject option; barrnap 1.x has neither."""
    p = subprocess.run(["barrnap", "--version"], capture_output=True, text=True)
    out = (p.stdout + p.stderr).strip()   # 0.9 prints it on stderr, 1.x on stdout
    if p.returncode != 0 or out.split()[:2] != ["barrnap", "0.9"]:
        sys.exit(f"barrnap 0.9 is required, found: {out.strip()!r}")


def barrnap(kingdom: str, fasta: str, threads: int, min_fraction: float) -> str:
    reject = min(BARRNAP_REJECT, min_fraction)
    return run(["barrnap", "--kingdom", kingdom, "--threads", str(threads),
                "--reject", str(reject), "--quiet", fasta])


def measure(fasta: str, threads: int, min_fraction: float) -> dict:
    rrna = best_of(count_rrna(barrnap("bac", fasta, threads, min_fraction), min_fraction),
                   count_rrna(barrnap("arc", fasta, threads, min_fraction), min_fraction))
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "trna.txt")
        run(["tRNAscan-SE", "-G", "-Q", "--thread", str(threads), "-o", out, fasta])
        aa = trna_amino_acids(open(out).read()) if os.path.exists(out) else 0
    return {**rrna, "tRNA": aa}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--bins-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--min-fraction", type=float, default=0.5,
                   help="fraction of a gene's expected length that must align (all rRNA genes)")
    a = p.parse_args(argv)
    if not 0 < a.min_fraction <= 1:
        sys.exit(f"--min-fraction must be in (0, 1], got {a.min_fraction}")
    check_barrnap_version()
    bins = sorted(glob.glob(os.path.join(a.bins_dir, "*.fa")))
    with open(a.out + ".partial", "w") as fh:
        fh.write(f"# barrnap (bac and arc, --reject {min(BARRNAP_REJECT, a.min_fraction)}; "
                 f"a gene counts at >= {a.min_fraction:g} of expected length aligned, "
                 f"same for 5S, 16S, 23S; larger count of the two models per gene) "
                 f"and tRNAscan-SE -G; all {len(bins)} MAGs measured\n")
        fh.write("bin\trRNA_5S\trRNA_16S\trRNA_23S\ttRNA_Amino_Acids\n")
        for f in bins:
            m = measure(f, a.threads, a.min_fraction)
            name = os.path.splitext(os.path.basename(f))[0]
            fh.write(f"{name}\t{m['5S_rRNA']}\t{m['16S_rRNA']}\t{m['23S_rRNA']}\t{m['tRNA']}\n")
    os.replace(a.out + ".partial", a.out)


if __name__ == "__main__":
    main()
