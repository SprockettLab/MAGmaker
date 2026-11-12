#!/usr/bin/env python3
"""A long-read assembly's contigs, renamed, with a table of what the assembler said.

Contigs are renamed <sample>__<assembler>__<original name>: Flye names contigs
contig_1, contig_2 ... in every sample, so names would collide when assemblies
are pooled. The table records, per contig: new name, original name, length
measured from the sequence, circular (yes, no, possibly) and depth.

  flye     circular and depth from assembly_info.txt (columns 4 and 3; Y = yes)
  myloasm  from the header, e.g.
           >u2658272ctg_len-2275750_circular-possibly_depth-10-9-9_duplicated-no mult=1.00
           only "yes" counts as circular downstream; "possibly" is kept as written.

Every contig must appear in both the FASTA and the assembler's own report, and
the length the assembler reports must equal the sequence length; otherwise
this stops.
"""

from __future__ import annotations

import argparse
import re
import sys


class ContigError(ValueError):
    pass


def read_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def flye_info(path) -> dict[str, dict]:
    """assembly_info.txt: #seq_name length cov. circ. repeat mult. alt_group graph_path"""
    info = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            info[f[0]] = {"length": int(f[1]), "depth": f[2], "circular": "yes" if f[3] == "Y" else "no"}
    return info


MYLOASM = re.compile(r"_len-(\d+)_circular-(yes|no|possibly)_depth-([^_\s]+)")


def myloasm_info(name: str) -> dict:
    m = MYLOASM.search(name)
    if not m:
        raise ContigError(f"myloasm contig name not understood: {name}")
    return {"length": int(m.group(1)), "circular": m.group(2), "depth": m.group(3)}


def contigs(sample, assembler, assembly, flye_info_path=None):
    """Yield (new name, original name, sequence, circular, depth), checked."""
    info = flye_info(flye_info_path) if assembler == "flye" else None
    seen = set()
    for name, seq in read_fasta(assembly):
        if assembler == "flye":
            if name not in info:
                raise ContigError(f"contig {name} is in the FASTA but not in assembly_info.txt")
            rec = info[name]
        else:
            rec = myloasm_info(name)
        if rec["length"] != len(seq):
            raise ContigError(f"length of {name}: assembler reports {rec['length']}, sequence has {len(seq)}")
        seen.add(name)
        yield f"{sample}__{assembler}__{name}", name, seq, rec["circular"], rec["depth"]
    if info is not None and set(info) - seen:
        raise ContigError(f"in assembly_info.txt but not the FASTA: {sorted(set(info) - seen)[:5]}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--sample", required=True)
    p.add_argument("--assembler", required=True, choices=["flye", "myloasm"])
    p.add_argument("--assembly", required=True)
    p.add_argument("--flye-info")
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-info", required=True)
    a = p.parse_args(argv)
    if a.assembler == "flye" and not a.flye_info:
        sys.exit("--flye-info is required for flye")
    n = 0
    with open(a.out_fasta, "w") as fa, open(a.out_info, "w") as tsv:
        tsv.write("contig\toriginal_name\tlength\tcircular\tdepth\n")
        for new, orig, seq, circ, depth in contigs(a.sample, a.assembler, a.assembly, a.flye_info):
            fa.write(f">{new}\n")
            for i in range(0, len(seq), 80):
                fa.write(seq[i:i + 80] + "\n")
            tsv.write(f"{new}\t{orig}\t{len(seq)}\t{circ}\t{depth}\n")
            n += 1
    if n == 0:
        sys.exit(f"{a.assembly} holds no contigs")


if __name__ == "__main__":
    main()
