#!/usr/bin/env python3
"""rRNA and tRNA genes per MAG, for the full MIMAG high-quality standard.

For every *.fa in a bin directory:
  rRNA  barrnap with the bacterial and with the archaeal model; complete genes
        only (records marked partial are not counted), except 5S when
        --count-partial-5s is given (barrnap often marks real 5S genes partial,
        being short); per gene (5S, 16S, 23S) the larger of the two models' counts.
  tRNA  tRNAscan-SE, general model (-G): the number of the 20 standard amino
        acids with at least one tRNA.
Every MAG in the directory is measured. Writes one row per MAG; an empty
directory gives a header-only table.
"""

from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
import tempfile

RRNA = ("5S_rRNA", "16S_rRNA", "23S_rRNA")
AMINO_ACIDS = {"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
               "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"}


def complete_rrna(gff_text: str, partial_ok: tuple = ()) -> dict[str, int]:
    """Complete 5S, 16S and 23S genes in one barrnap GFF; genes named in
    partial_ok are counted even when barrnap marks them partial."""
    n = {g: 0 for g in RRNA}
    for line in gff_text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        name = f[8].split(";")[0].replace("Name=", "")
        if name in n and ("partial" not in f[8] or name in partial_ok):
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
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def measure(fasta: str, threads: int, partial_ok: tuple = ()) -> dict:
    rrna = best_of(
        complete_rrna(run(["barrnap", "--kingdom", "bac", "--threads", str(threads), "--quiet", fasta]), partial_ok),
        complete_rrna(run(["barrnap", "--kingdom", "arc", "--threads", str(threads), "--quiet", fasta]), partial_ok))
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "trna.txt")
        subprocess.run(["tRNAscan-SE", "-G", "-Q", "--thread", str(threads), "-o", out, fasta],
                       capture_output=True, text=True, check=True)
        aa = trna_amino_acids(open(out).read()) if os.path.exists(out) else 0
    return {**rrna, "tRNA": aa}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--bins-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--count-partial-5s", action="store_true")
    a = p.parse_args(argv)
    partial_ok = ("5S_rRNA",) if a.count_partial_5s else ()
    bins = sorted(glob.glob(os.path.join(a.bins_dir, "*.fa")))
    with open(a.out + ".partial", "w") as fh:
        rule = "complete genes, partial 5S also counted" if partial_ok else "complete genes only"
        fh.write(f"# barrnap (bac and arc, {rule}, larger count per gene) and tRNAscan-SE -G; "
                 f"all {len(bins)} MAGs measured\n")
        fh.write("bin\trRNA_5S\trRNA_16S\trRNA_23S\ttRNA_Amino_Acids\n")
        for f in bins:
            m = measure(f, a.threads, partial_ok)
            name = os.path.splitext(os.path.basename(f))[0]
            fh.write(f"{name}\t{m['5S_rRNA']}\t{m['16S_rRNA']}\t{m['23S_rRNA']}\t{m['tRNA']}\n")
    os.replace(a.out + ".partial", a.out)


if __name__ == "__main__":
    main()
