"""Merge every binning row's ISEScan table into one, labelled with its MAG.

Inputs (from rule make_isescan_summary):
  input.tables       output/isescan/<assembler>/<row>/is.tsv, one per row;
                     ISEScan's own tab-separated table, or empty when its
                     HMM search found nothing (NO_IS_FOUND beside it)
  input.mag_summary  output/mag_qc/mag_summary.tsv
  params.pairs       [(assembler, row)] in the same order as input.tables
  params.row_sample  {row: sample}

Output: output/isescan/isescan_summary.tsv. One row per IS element copy, with
Sample_ID, Assembler, Binning_Row and Contig in front of ISEScan's own
columns (kept as written), then the MAG its contig was binned into:
MAG_ID, MAG_New_Name, MAG_GTDB_Classification, MAG_MIMAG_Quality (NA when
the contig is in no MAG; ';'-joined when it is in MAGs from more than one
mapper).

Contig names: ISEScan keys each sequence by its whole FASTA header, so seqID
can carry MEGAHIT's 'flag=... multi=... len=...'. Contig is the first word,
which is how the bins name the same contig.

MAG membership comes from the MAG FASTAs named in mag_summary.tsv's
Original_Path. Contig names repeat between samples (MEGAHIT numbers each
assembly from k141_0), so a MAG's contigs are attached only to the binning
row its FASTA sits under (the folder name in Original_Path).
"""
import os
import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], 'w')

MAG_COLS = ['MAG_ID', 'New_Name', 'GTDB_Classification', 'MIMAG_Quality']


def first_word(header):
    return str(header).split()[0] if str(header).strip() else ''


def mag_membership(mag_summary):
    """{(row, contig): {column: value}} for every contig of every MAG."""
    mags = pd.read_csv(mag_summary, sep='\t', dtype=str)
    out = {}
    n_mags = n_missing = 0
    for m in mags.itertuples(index=False):
        path = getattr(m, 'Original_Path', None)
        if not isinstance(path, str) or path in ('', 'NA'):
            continue
        if not os.path.exists(path):
            n_missing += 1
            print(f"MAG FASTA not found, its contigs stay unlabelled: {path}", file=sys.stderr)
            continue
        n_mags += 1
        row = os.path.basename(os.path.dirname(path))
        info = {c: getattr(m, c, 'NA') for c in MAG_COLS}
        with open(path) as h:
            for line in h:
                if line.startswith('>'):
                    key = (row, first_word(line[1:]))
                    if key in out:
                        for c in MAG_COLS:
                            out[key][c] = out[key][c] + ';' + str(info[c])
                    else:
                        out[key] = dict(info)
    print(f"MAGs read: {n_mags}; MAG FASTAs missing: {n_missing}; "
          f"contigs in MAGs: {len(out)}", file=sys.stderr)
    return out


def main():
    members = mag_membership(snakemake.input.mag_summary)
    frames = []
    for (assembler, row), path in zip(snakemake.params.pairs, snakemake.input.tables):
        if os.path.getsize(path) == 0:
            print(f"{row} ({assembler}): no IS elements (empty table)", file=sys.stderr)
            continue
        t = pd.read_csv(path, sep='\t', dtype=str)
        if 'seqID' not in t.columns:
            raise ValueError(f"{path}: no seqID column; columns are {list(t.columns)}")
        t.insert(0, 'Contig', t['seqID'].map(first_word))
        t.insert(0, 'Binning_Row', row)
        t.insert(0, 'Assembler', assembler)
        t.insert(0, 'Sample_ID', snakemake.params.row_sample[row])
        lab = [members.get((row, c)) for c in t['Contig']]
        for c in MAG_COLS:
            name = c if c == 'MAG_ID' else 'MAG_' + c
            t[name] = [x[c] if x else 'NA' for x in lab]
        n_in = sum(x is not None for x in lab)
        print(f"{row} ({assembler}): {len(t)} IS copies, {n_in} on a binned contig", file=sys.stderr)
        frames.append(t)
    if frames:
        out = pd.concat(frames, ignore_index=True)
    else:
        out = pd.DataFrame(columns=['Sample_ID', 'Assembler', 'Binning_Row', 'Contig', 'seqID']
                           + ['MAG_ID', 'MAG_New_Name', 'MAG_GTDB_Classification', 'MAG_MIMAG_Quality'])
    out.to_csv(snakemake.output.summary, sep='\t', index=False)
    print(f"total: {len(out)} IS copies from {len(frames)} of {len(snakemake.input.tables)} tables",
          file=sys.stderr)


main()
