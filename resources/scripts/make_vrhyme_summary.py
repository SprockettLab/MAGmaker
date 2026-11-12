#!/usr/bin/env python3
"""Merge every assembly's vRhyme bins into two tables.

vrhyme_summary.tsv   one row per vMAG: how many contigs it holds, their total
                     length, and vRhyme's own protein redundancy count.
vrhyme_membership.tsv one row per contig: which vMAG it went into.

The number that decides whether vRhyme earned its place is how many contigs
ended up in a MULTI-contig bin, since a bin of one contig is just that contig
renamed. Both tables carry it, and the log reports it with its denominator.

Circular scaffolds are set aside by vRhyme BEFORE binning (DTR/ITR detection).
Those are complete genomes, so they are listed in the membership table as
their own single-contig bins marked circular, rather than going missing
between binning and whatever reads these tables next.

Run by vrhyme.smk's vrhyme_summary rule through the `snakemake` object.
"""

import csv
import glob
import os
import sys

log_path = str(snakemake.log[0])
os.makedirs(os.path.dirname(log_path) or '.', exist_ok=True)
sys.stderr = open(log_path, 'w')
sys.stdout = sys.stderr

csv.field_size_limit(10 * 1024 * 1024)

pairs = list(snakemake.params.pairs or [])
row_sample = dict(snakemake.params.row_sample or {})
base = snakemake.params.base


def one_file(directory, pattern):
    """The single file matching pattern under directory, or None.

    vRhyme writes one membership and one summary per iteration under
    vRhyme_alternate_bins/, so the search excludes that folder; more than one
    match left over means the layout is not what this expects, and reporting
    the first as if it described the run would be wrong.
    """
    hits = [p for p in glob.glob(os.path.join(directory, '**', pattern),
                                 recursive=True)
            if 'alternate' not in p]
    if len(hits) > 1:
        sys.exit('ERROR: {0} files match {1} under {2}: {3}'.format(
            len(hits), pattern, directory, sorted(hits)))
    return hits[0] if hits else None


def rows_of(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def value(row, *names):
    for name in names:
        if name in row and str(row[name]).strip():
            return str(row[name]).strip()
    return ''


summary_rows, membership_rows = [], []
assemblies_with_bins = 0

for assembler, row in pairs:
    sample = row_sample.get(row, row)
    directory = os.path.join(base, assembler, row)
    membership_file = one_file(directory, 'vRhyme_best_bins.*.membership.tsv')
    summary_file = one_file(directory, 'vRhyme_best_bins.*.summary.tsv')
    circular_file = one_file(directory, '*circular.tsv')

    members = {}
    for entry in rows_of(membership_file):
        fields = list(entry.values())
        contig = value(entry, 'scaffold', 'sequence', 'contig') or \
            (fields[0] if fields else '')
        bin_id = value(entry, 'bin', 'cluster') or (fields[-1] if fields else '')
        if not contig or not bin_id:
            continue
        members.setdefault(bin_id, []).append(contig)

    if members:
        assemblies_with_bins += 1

    redundancy = {}
    for entry in rows_of(summary_file):
        fields = list(entry.values())
        bin_id = value(entry, 'bin', 'cluster') or (fields[0] if fields else '')
        redundancy[bin_id] = value(entry, 'redundancy', 'redundant_proteins')

    for bin_id, contigs in sorted(members.items()):
        name = 'vMAG_{0}_{1}_{2}'.format(sample, assembler, bin_id)
        summary_rows.append({
            'vMAG_ID': name,
            'Sample_ID': sample,
            'Assembler': assembler,
            'Binning_Row': row,
            'vRhyme_Bin': bin_id,
            'N_Contigs': len(contigs),
            'Multi_Contig': 'yes' if len(contigs) > 1 else 'no',
            'Circular': 'no',
            'Redundant_Proteins': redundancy.get(bin_id, ''),
        })
        for contig in contigs:
            membership_rows.append({
                'vMAG_ID': name, 'Sample_ID': sample, 'Assembler': assembler,
                'Binning_Row': row, 'Seq_Name': contig, 'Circular': 'no'})

    # Circular scaffolds, set aside before binning: complete genomes, each
    # its own vMAG so they are not lost.
    for entry in rows_of(circular_file):
        fields = list(entry.values())
        contig = value(entry, 'scaffold', 'sequence', 'contig') or \
            (fields[0] if fields else '')
        if not contig:
            continue
        name = 'vMAG_{0}_{1}_circular_{2}'.format(sample, assembler, contig)
        summary_rows.append({
            'vMAG_ID': name, 'Sample_ID': sample, 'Assembler': assembler,
            'Binning_Row': row, 'vRhyme_Bin': 'circular', 'N_Contigs': 1,
            'Multi_Contig': 'no', 'Circular': 'yes', 'Redundant_Proteins': ''})
        membership_rows.append({
            'vMAG_ID': name, 'Sample_ID': sample, 'Assembler': assembler,
            'Binning_Row': row, 'Seq_Name': contig, 'Circular': 'yes'})

summary_columns = ['vMAG_ID', 'Sample_ID', 'Assembler', 'Binning_Row',
                   'vRhyme_Bin', 'N_Contigs', 'Multi_Contig', 'Circular',
                   'Redundant_Proteins']
membership_columns = ['vMAG_ID', 'Sample_ID', 'Assembler', 'Binning_Row',
                      'Seq_Name', 'Circular']


def write(path, columns, rows):
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


write(snakemake.output.summary, summary_columns, summary_rows)
write(snakemake.output.membership, membership_columns, membership_rows)

multi = [r for r in summary_rows if r['Multi_Contig'] == 'yes']
in_multi = sum(r['N_Contigs'] for r in multi)
contigs = sum(r['N_Contigs'] for r in summary_rows)
print('{0} vMAGs from {1} of {2} assemblies'.format(
    len(summary_rows), assemblies_with_bins, len(pairs)))
if summary_rows:
    print('{0} of {1} vMAGs hold more than one contig ({2:.1f}%)'.format(
        len(multi), len(summary_rows), 100.0 * len(multi) / len(summary_rows)))
if contigs:
    print('{0} of {1} binned contigs are in a multi-contig vMAG ({2:.1f}%)'.format(
        in_multi, contigs, 100.0 * in_multi / contigs))
print('{0} circular vMAGs set aside before binning'.format(
    sum(1 for r in summary_rows if r['Circular'] == 'yes')))
