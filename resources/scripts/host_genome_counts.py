#!/usr/bin/env python3
"""Host reads counted by the genome they came from.

A BoarderCollie host can join several genomes into one reference (a primate
and a parasite, say). host_filter removes reads that hit any of them; this
splits that removed count by member genome, from the host BAM host_filter
already keeps. Three subcommands:

  contig-map   which member genome each reference sequence belongs to,
               from BoarderCollie's members.tsv and each member's cleaned
               FASTA (or one name for every sequence, for a host FASTA not
               made by BoarderCollie)
  count        one sample: SAM text on stdin (samtools view -h of the host
               BAM), written as one row per member genome
  merge        every sample's rows in one table, plus a parameters file

What "removed" means here is exactly what host_filter removed:
  paired-end   pairs bowtie2 aligned concordantly (host_filter keeps
               --un-conc, the pairs that did NOT align concordantly)
  single-end   reads bowtie2 aligned (host_filter keeps --un)
Each removed unit is credited to the member holding mate 1's sequence.
The per-sample total is checked against bowtie2's own summary in the
host_filter log, and the run stops if they differ.

Mates that aligned without their pair aligning concordantly stayed in the
nonhost reads; they are counted separately (Mates_Aligned_Not_Removed), so
a genome's reads left in the assembly input are visible too.

Every record of the BAM is read; nothing is sampled.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import sys

FLAG_PAIRED = 0x1
FLAG_PROPER = 0x2
FLAG_UNMAPPED = 0x4
FLAG_MATE_UNMAPPED = 0x8
FLAG_READ1 = 0x40
FLAG_READ2 = 0x80
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800

COLUMNS = ["Sample", "Layout", "Unit", "Member", "Genome", "Input",
           "Removed", "Removed_Fraction", "Min_MAPQ", "Removed_Confident",
           "Mates_Aligned_Not_Removed"]


# ── contig-map ──────────────────────────────────────────────────────────────

def fasta_names(path):
    """in : FASTA path (plain text)
       out: list of sequence names (first word of each header), in order"""
    names = []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                names.append(line[1:].split(None, 1)[0].decode())
    return names


def read_members(path):
    """in : BoarderCollie members.tsv
       out: list of (member accession, genome, cleaned FASTA path)

    A mitochondrion added to a genome is its own member, with member_of
    naming that genome; its reads count towards the genome too, in the
    Genome column, while Member keeps them apart."""
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        need = {"accession", "clean_genome", "member_of"}
        if not need <= set(header):
            raise ValueError(f"{path}: needs columns {sorted(need)}, "
                             f"has {header}")
        col = {name: i for i, name in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            acc, fasta, of = (f[col["accession"]], f[col["clean_genome"]],
                              f[col["member_of"]])
            rows.append((acc, acc if of in ("", "-") else of, fasta))
    if not rows:
        raise ValueError(f"{path} lists no members")
    return rows


def contig_map(members=None, reference=None, single_name=None):
    """in : members.tsv path, or a reference FASTA plus one name for it
       out: list of (contig, member, genome) in reference order

    BoarderCollie refuses a sequence name that occurs in two members, but
    this checks again: a name in two members could not be assigned."""
    out = []
    if members:
        for acc, genome, fasta in read_members(members):
            if not os.path.exists(fasta):
                raise ValueError(f"member {acc}: cleaned FASTA {fasta} "
                                 "(from members.tsv) does not exist")
            out.extend((name, acc, genome) for name in fasta_names(fasta))
    else:
        out = [(name, single_name, single_name)
               for name in fasta_names(reference)]
    seen = set()
    for name, _, _ in out:
        if name in seen:
            raise ValueError(f"sequence name {name} occurs twice")
        seen.add(name)
    if not out:
        raise ValueError("the host reference has no sequences")
    return out


def cmd_contig_map(a):
    rows = contig_map(a.members, a.reference, a.single_name)
    with open(a.out, "w") as fh:
        fh.write("Contig\tMember\tGenome\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")


# ── count ───────────────────────────────────────────────────────────────────

def load_map(path):
    """in : contig map TSV from contig-map
       out: ({contig: member}, {member: genome}, [members in file order])"""
    contig_member, member_genome, order = {}, {}, []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            contig, member, genome = line.rstrip("\n").split("\t")
            contig_member[contig] = member
            if member not in member_genome:
                member_genome[member] = genome
                order.append(member)
    return contig_member, member_genome, order


BT2_TOTAL = re.compile(r"^(\d+) reads; of these:")
BT2_CONC = re.compile(r"^\s*(\d+) \([\d.]+%\) aligned concordantly "
                      r"(exactly 1 time|>1 times)")
BT2_SE = re.compile(r"^\s*(\d+) \([\d.]+%\) aligned (exactly 1 time|>1 times)")


def bowtie2_summary(path, layout):
    """in : host_filter log (bowtie2's stderr) and PE or SE
       out: (input units, units bowtie2 says it removed)

    PE: input pairs, and pairs aligned concordantly (exactly once + more
    than once). SE: input reads, and reads aligned. Only the first summary
    block is read; host_filter writes one."""
    total, removed, found = None, 0, 0
    pattern = BT2_CONC if layout == "PE" else BT2_SE
    with open(path) as fh:
        for line in fh:
            if total is None:
                m = BT2_TOTAL.match(line)
                if m:
                    total = int(m.group(1))
                continue
            m = pattern.match(line)
            if m:
                removed += int(m.group(1))
                found += 1
                if found == 2:
                    break
    if total is None or found != 2:
        raise ValueError(f"{path}: no bowtie2 {layout} alignment summary "
                         "(expected 'N reads; of these:' and the two "
                         "aligned-1-time / >1-times lines)")
    return total, removed


def sam_records(stream, sq_names):
    """in : SAM text lines; a set that receives every @SQ name
       out: yields (qname, flag, rname, mapq) for each alignment record"""
    for line in stream:
        if line.startswith("@"):
            if line.startswith("@SQ"):
                for field in line.rstrip("\n").split("\t")[1:]:
                    if field.startswith("SN:"):
                        sq_names.add(field[3:])
            continue
        f = line.split("\t", 6)
        flag = int(f[1])
        if flag & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY):
            raise ValueError(
                f"read {f[0]}: secondary or supplementary record (flag "
                f"{flag}); host_filter's bowtie2 call writes neither, so "
                "this BAM did not come from it")
        yield f[0], flag, f[2], int(f[4])


def count(stream, layout, contig_member, members, min_mapq):
    """in : SAM lines, PE or SE, {contig: member}, member names, MAPQ cutoff
       out: (per-member {removed, confident, mates_not_removed}, @SQ names)

    PE records come in mate pairs, mate 1 then mate 2 (bowtie2's order,
    unsorted, as host_filter writes it). A pair is removed when both mates
    carry the proper-pair flag, bowtie2's mark of a concordant alignment;
    it is confident when the lower of the two MAPQs reaches min_mapq.
    With --no-unal, bowtie2 leaves out the record of a mate that did not
    align, so a pair with one aligned mate is one record: aligned, flagged
    mate-unmapped, never proper. It counts as an aligned mate not removed."""
    tally = {m: {"removed": 0, "confident": 0, "mates": 0} for m in members}
    sq = set()

    def member_of(contig, qname):
        try:
            return contig_member[contig]
        except KeyError:
            raise ValueError(f"read {qname} aligned to {contig}, which is "
                             "not in the contig map: the map was made from "
                             "another reference") from None

    records = sam_records(stream, sq)
    if layout == "SE":
        for qname, flag, rname, mapq in records:
            if flag & FLAG_PAIRED:
                raise ValueError(f"read {qname} is paired in a single-end BAM")
            if flag & FLAG_UNMAPPED:
                continue
            t = tally[member_of(rname, qname)]
            t["removed"] += 1
            if mapq >= min_mapq:
                t["confident"] += 1
        return tally, sq

    def lone(rec):
        """A record whose mate's record bowtie2 left out (--no-unal)."""
        q, f, r, _ = rec
        if f & FLAG_UNMAPPED or not f & FLAG_MATE_UNMAPPED or f & FLAG_PROPER:
            raise ValueError(
                f"read {q} (flag {f}) has no mate record next to it, but is "
                "not an aligned mate of an unaligned one; the BAM must be "
                "bowtie2's unsorted output, as host_filter writes it")
        tally[member_of(r, q)]["mates"] += 1

    pending = None
    while True:
        first = pending if pending is not None else next(records, None)
        pending = None
        if first is None:
            break
        if not first[1] & FLAG_PAIRED:
            raise ValueError(f"read {first[0]} is unpaired in a paired-end BAM")
        second = next(records, None)
        if second is None or second[0] != first[0]:
            lone(first)
            pending = second
            continue
        q1, f1, r1, mq1 = first
        q2, f2, r2, mq2 = second
        if q1 != q2 or not (f1 & FLAG_READ1 and f2 & FLAG_READ2):
            raise ValueError(
                f"records {q1} (flag {f1}) and {q2} (flag {f2}) are not "
                "mate 1 and mate 2 of one pair; the BAM must be bowtie2's "
                "unsorted output, as host_filter writes it")
        if f1 & FLAG_PROPER and f2 & FLAG_PROPER:
            if f1 & FLAG_UNMAPPED or f2 & FLAG_UNMAPPED or r1 != r2:
                raise ValueError(f"pair {q1}: proper-pair flag set, but the "
                                 "mates are not both aligned to one sequence")
            t = tally[member_of(r1, q1)]
            t["removed"] += 1
            if min(mq1, mq2) >= min_mapq:
                t["confident"] += 1
        else:
            for qname, flag, rname in ((q1, f1, r1), (q2, f2, r2)):
                if not flag & FLAG_UNMAPPED:
                    tally[member_of(rname, qname)]["mates"] += 1
    return tally, sq


def fraction(part, whole):
    return "NA" if whole == 0 else f"{part / whole:.6g}"


def cmd_count(a):
    contig_member, member_genome, order = load_map(a.contig_map)
    total, bt2_removed = bowtie2_summary(a.bowtie2_log, a.layout)
    tally, sq = count(sys.stdin, a.layout, contig_member, order, a.min_mapq)

    # The BAM's own reference must be the one the map describes.
    if sq != set(contig_member):
        missing = len(set(contig_member) - sq)
        extra = len(sq - set(contig_member))
        raise ValueError(
            f"the BAM header's {len(sq)} sequences and the contig map's "
            f"{len(contig_member)} differ ({extra} only in the BAM, "
            f"{missing} only in the map): another reference")
    # The split must add up to what host_filter removed.
    removed = sum(t["removed"] for t in tally.values())
    if removed != bt2_removed:
        raise ValueError(
            f"{removed} {('pairs' if a.layout == 'PE' else 'reads')} counted "
            f"as removed in the BAM, but bowtie2's summary in "
            f"{a.bowtie2_log} says {bt2_removed}")

    unit = "pairs" if a.layout == "PE" else "reads"
    with open(a.out, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for m in order:
            t = tally[m]
            fh.write("\t".join(str(x) for x in [
                a.sample, a.layout, unit, m, member_genome[m], total,
                t["removed"], fraction(t["removed"], total), a.min_mapq,
                t["confident"],
                t["mates"] if a.layout == "PE" else "NA"]) + "\n")


# ── merge ───────────────────────────────────────────────────────────────────

def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def cmd_merge(a):
    n_rows = 0
    with open(a.out, "w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for path in a.counts:
            with open(path) as fh:
                header = fh.readline().rstrip("\n").split("\t")
                if header != COLUMNS:
                    raise ValueError(f"{path}: unexpected header {header}")
                for line in fh:
                    out.write(line)
                    n_rows += 1
    with open(a.params, "w") as fh:
        for key, value in [
            ("table", a.out),
            ("samples", len(a.counts)),
            ("rows", n_rows),
            ("host", a.host),
            ("reference", a.reference),
            ("members_tsv", a.members or "none (host FASTA not from "
                                         "BoarderCollie; one member)"),
            ("members_tsv_md5", md5(a.members) if a.members else "NA"),
            ("min_mapq", a.min_mapq),
            ("removed_pe", "pairs bowtie2 aligned concordantly (proper-pair "
                           "flag on both mates); what host_filter removed"),
            ("removed_se", "reads bowtie2 aligned; what host_filter_se "
                           "removed"),
            ("assignment", "the member holding mate 1's sequence"),
            ("removed_confident", "removed units whose MAPQ (the lower of "
                                  "the two mates for pairs) is >= min_mapq"),
            ("mates_aligned_not_removed", "PE only: aligned mates of pairs "
                                          "that did not align concordantly; "
                                          "these stayed in the nonhost reads"),
            ("input", "pairs (PE) or reads (SE) given to bowtie2, from its "
                      "summary in the host_filter log"),
            ("check", "per sample, the removed total equals bowtie2's own "
                      "summary, and the BAM's @SQ names equal the contig map"),
            ("sampling", "none; every record of every host BAM is read"),
        ]:
            fh.write(f"{key}\t{value}\n")


# ── command line ────────────────────────────────────────────────────────────

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("contig-map")
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--members", help="BoarderCollie members.tsv")
    g.add_argument("--reference", help="host FASTA not made by BoarderCollie")
    s.add_argument("--single-name", help="member name with --reference")
    s.add_argument("--out", required=True)

    s = sub.add_parser("count")
    s.add_argument("--sample", required=True)
    s.add_argument("--layout", choices=["PE", "SE"], required=True)
    s.add_argument("--contig-map", required=True)
    s.add_argument("--bowtie2-log", required=True)
    s.add_argument("--min-mapq", type=int, required=True)
    s.add_argument("--out", required=True)

    s = sub.add_parser("merge")
    s.add_argument("--out", required=True)
    s.add_argument("--params", required=True)
    s.add_argument("--host", required=True)
    s.add_argument("--reference", required=True)
    s.add_argument("--members")
    s.add_argument("--min-mapq", type=int, required=True)
    s.add_argument("counts", nargs="+")

    a = p.parse_args(argv)
    if a.cmd == "contig-map":
        if a.reference and not a.single_name:
            p.error("--reference needs --single-name")
        cmd_contig_map(a)
    elif a.cmd == "count":
        cmd_count(a)
    else:
        cmd_merge(a)


if __name__ == "__main__":
    try:
        main()
    except ValueError as err:
        print(f"host_genome_counts: {err}", file=sys.stderr)
        sys.exit(1)
