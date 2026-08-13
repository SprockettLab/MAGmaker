"""Roll CMSeq's per-reference output up to one row per MAG.

in : poly.py's raw table and a contig-to-MAG map
out: strain_heterogeneity.tsv, one row per MAG

CMSeq reports per reference sequence, and a MAG is many contigs, so the
rates have to be recombined weighted by how many positions each contig
contributed. An unweighted mean over contigs would let a 2 kb contig with
three evaluated positions count as much as a 200 kb one.

poly.py writes a HASH-separated table, which is not documented in its help
text, so the delimiter and columns are located rather than assumed and a
departure from the expected layout is logged rather than guessed through.
Anything unparseable is reported as NA rather than as zero: zero would read
as "no strain heterogeneity", which is a claim, while NA says the
measurement was not made.

Under-covered MAGs are NA for the same reason. Sanders 2023 and Pasolli
2019 require more than 100 positions covered at 10x with base quality
above 30 before reporting a value, and that threshold is what makes these
numbers comparable with the published databases.
"""

import csv
import os
import re
import sys


def log(handle, message):
    handle.write(message + "\n")


def find_column(header, *patterns):
    """First header index whose name matches any pattern, else None."""
    for i, name in enumerate(header):
        low = name.strip().lower()
        for p in patterns:
            if re.search(p, low):
                return i
    return None


def parse_poly(path, logfh):
    """in : poly.py's raw output
       out: {reference: (polymorphic_positions, evaluated_positions)}

    poly.py writes a HASH-separated table, not tab-separated, with a header:

        referenceID#total_covered_bases#total_polymorphic_bases#
        total_polymorphic_rate#dominant_allele_distr_mean#...

    References with no covered bases have empty trailing fields, which is
    why the counts are read individually rather than by splitting into a
    fixed number of columns.

    Returns an empty dict rather than raising when the layout cannot be
    established, so one unreadable sample cannot stop a whole run.
    """
    lines = [l.rstrip("\n") for l in open(path) if l.strip()]
    if not lines:
        log(logfh, "poly.py output is empty")
        return {}

    log(logfh, "first lines of poly.py output:")
    for l in lines[:3]:
        log(logfh, "  " + l[:160])

    # poly.py uses '#'. Fall back to tab or comma if a future version
    # changes it, rather than silently reading one giant column.
    delim = "#"
    for candidate in ("#", "\t", ","):
        if candidate in lines[0]:
            delim = candidate
            break
    if delim != "#":
        log(logfh, "WARNING: delimiter looks like %r, not '#'" % delim)

    header = [c.strip().lower() for c in lines[0].split(delim)]
    if "referenceid" not in header[0]:
        log(logfh, "WARNING: first column is %r, expected referenceID; "
                   "treating the file as headerless" % header[0])
        header = None
        body = lines
    else:
        body = lines[1:]

    if header:
        try:
            cov_i = header.index("total_covered_bases")
            pol_i = header.index("total_polymorphic_bases")
        except ValueError:
            cov_i = find_column(header, r"covered.*base", r"total.*covered")
            pol_i = find_column(header, r"polymorphic.*base")
            log(logfh, "exact column names absent; matched covered=%s "
                       "polymorphic=%s" % (cov_i, pol_i))
    else:
        cov_i, pol_i = 1, 2   # the documented order

    if cov_i is None or pol_i is None:
        log(logfh, "WARNING: could not locate the count columns; no rows read")
        return {}

    def as_float(cell):
        cell = (cell or "").strip()
        return float(cell) if cell else 0.0

    out = {}
    skipped = 0
    for line in body:
        parts = line.split(delim)
        if not parts or not parts[0].strip():
            continue
        try:
            out[parts[0].strip()] = (as_float(parts[pol_i]),
                                     as_float(parts[cov_i]))
        except (IndexError, ValueError):
            skipped += 1
    if skipped:
        log(logfh, "%d rows could not be read" % skipped)
    log(logfh, "%d references parsed" % len(out))
    return out


def main():
    raw = snakemake.input.raw               # noqa: F821
    cmap = snakemake.input.map              # noqa: F821
    out = snakemake.output.tsv              # noqa: F821
    min_positions = float(snakemake.params.min_positions)   # noqa: F821

    with open(snakemake.log[0], "w") as logfh:   # noqa: F821
        contig_to_mag = {}
        if os.path.exists(cmap):
            for line in open(cmap):
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0]:
                    contig_to_mag[parts[0]] = parts[1]
        log(logfh, "%d contigs mapped to %d MAGs"
            % (len(contig_to_mag), len(set(contig_to_mag.values()))))

        per_contig = parse_poly(raw, logfh) if os.path.exists(raw) else {}

        totals = {}
        for contig, (poly, tot) in per_contig.items():
            mag = contig_to_mag.get(contig)
            if mag is None:
                continue
            p, t = totals.get(mag, (0.0, 0.0))
            totals[mag] = (p + poly, t + tot)

        rows = []
        for mag in sorted(set(contig_to_mag.values())):
            poly, tot = totals.get(mag, (0.0, 0.0))
            if tot > min_positions:
                rows.append({
                    "MAG": mag,
                    "Strain_Heterogeneity": round(100.0 * poly / tot, 4),
                    "SH_Positions_Evaluated": int(tot),
                })
            else:
                # Measured but under-covered is not the same as clean.
                rows.append({
                    "MAG": mag,
                    "Strain_Heterogeneity": "NA",
                    "SH_Positions_Evaluated": int(tot),
                })

        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(
                fh, fieldnames=["MAG", "Strain_Heterogeneity",
                                "SH_Positions_Evaluated"],
                delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)

        scored = sum(1 for r in rows if r["Strain_Heterogeneity"] != "NA")
        log(logfh, "%d MAGs written, %d with enough coverage to score "
                   "(>%g positions)" % (len(rows), scored, min_positions))


main()
