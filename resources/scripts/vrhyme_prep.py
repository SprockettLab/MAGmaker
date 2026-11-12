#!/usr/bin/env python3
"""Select geNomad's viral contigs for vRhyme and subset the coverage table.

Nothing is mapped and nothing is recomputed: the inputs are geNomad's viral
FASTA and the MetaBAT2 coverage table the MAG track already wrote.

Every count is written to prep_params.txt with its denominator, including how
many selected names were found in the FASTA and in the coverage table. A poor
join invalidates the binning, so it has to be visible rather than absorbed.

THE HEADERS ARE TRUNCATED TO THEIR FIRST TOKEN, and that is essential.
vRhyme keys scaffolds on the ENTIRE header line. geNomad writes
`>k141_602388 flag=1 multi=10.0000 len=2343` while the coverage table's
contigName is `k141_602388`, so with full headers not one coverage row matches
any scaffold: vRhyme finds no pairs, reports "binning complete" with 0 bins
and exits 0. WF30 lost a day to this on 2026-11-11, across 13 samples, six
--cov settings and four --penalty_n settings, all returning exactly 0 bins.
Truncating took one sample from 0 bins to 473. Upstream's example data hides
it, because its headers have no description.

Run by vrhyme.smk's vrhyme_prep rule through the `snakemake` object.
"""

import csv
import os
import sys

log_path = str(snakemake.log[0])
os.makedirs(os.path.dirname(log_path) or '.', exist_ok=True)
sys.stderr = open(log_path, 'w')
sys.stdout = sys.stderr

csv.field_size_limit(10 * 1024 * 1024)

QUALITY_ORDER = ['Not-determined', 'Low-quality', 'Medium-quality',
                 'High-quality', 'Complete']

virus_fna = snakemake.params.virus_fna
checkv_directory = snakemake.params.checkv_dir
coverage_in = snakemake.input.coverage
min_length = int(snakemake.params.min_length)
out_fasta = snakemake.output.fasta
out_coverage = snakemake.output.coverage
out_params = snakemake.output.params_file

if min_length < 2000:
    sys.exit('ERROR: params.vrhyme.min_length is {0}; vRhyme\'s floor is 2000 '
             'bp and upstream states it cannot be lowered.'.format(min_length))


def genomad_rows(summary_path):
    """seq_name -> (length, topology) from geNomad's virus summary."""
    rows = {}
    if not os.path.exists(summary_path) or os.path.getsize(summary_path) == 0:
        return rows
    with open(summary_path, newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            name = (row.get('seq_name') or '').strip()
            if not name:
                continue
            try:
                length = int(float(row.get('length') or 0))
            except ValueError:
                length = 0
            rows[name] = (length, (row.get('topology') or '').strip())
    return rows


def checkv_calls(directory):
    """contig_id -> (quality, provirus) from CheckV's quality summary."""
    calls = {}
    path = os.path.join(directory, 'quality_summary.tsv')
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return calls
    with open(path, newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            calls[(row.get('contig_id') or '').strip()] = (
                (row.get('checkv_quality') or '').strip() or 'Not-determined',
                (row.get('provirus') or '').strip())
    return calls


def records(path):
    """Yield (first-token name, record text with the header truncated)."""
    name, chunks = None, []
    if not os.path.exists(path):
        return
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                tokens = line[1:].strip().split()
                name = tokens[0] if tokens else 'unnamed'
                chunks = ['>' + name + '\n']
            else:
                chunks.append(line)
    if name is not None:
        yield name, ''.join(chunks)


summary_dir = os.path.dirname(virus_fna)
stem = os.path.basename(virus_fna)[:-len('_virus.fna')]
lengths = genomad_rows(os.path.join(summary_dir, stem + '_virus_summary.tsv'))
calls = checkv_calls(checkv_directory)

stats = {'genomad_virus_sequences': 0, 'excluded_provirus': 0,
         'excluded_short': 0, 'selected': 0}
tiers = {}
selected = set()

for name in lengths:
    stats['genomad_virus_sequences'] += 1
    quality, provirus = calls.get(name, ('Not-determined', ''))
    if '|provirus_' in name or provirus.strip().lower() in {'yes', 'true'}:
        stats['excluded_provirus'] += 1
        continue
    if lengths[name][0] < min_length:
        stats['excluded_short'] += 1
        continue
    selected.add(name)
    tiers[quality] = tiers.get(quality, 0) + 1
    stats['selected'] += 1

# geNomad's names carry |provirus_<start>_<end>; the coverage table and the
# assembly use the bare contig name, so both forms are matched.
wanted = {name: name.split('|')[0] for name in selected}

os.makedirs(os.path.dirname(out_fasta) or '.', exist_ok=True)
written, fasta_total = 0, 0
with open(out_fasta, 'w') as out:
    for name, record in records(virus_fna):
        fasta_total += 1
        if name in selected or name in wanted.values():
            out.write(record)
            written += 1

coverage_rows, coverage_total, columns = 0, 0, 0
contig_names = set(wanted.values())
with open(coverage_in) as src, open(out_coverage, 'w') as out:
    header = src.readline()
    out.write(header)
    columns = len(header.rstrip('\n').split('\t'))
    for line in src:
        coverage_total += 1
        if line.split('\t', 1)[0] in contig_names:
            out.write(line)
            coverage_rows += 1

with open(out_params, 'w') as handle:
    handle.write('# vRhyme input preparation (MAGmaker)\n')
    handle.write('genomad_virus_fasta\t{0}\n'.format(virus_fna))
    handle.write('checkv_dir\t{0}\n'.format(checkv_directory))
    handle.write('coverage_table\t{0}\n'.format(coverage_in))
    handle.write('min_length_bp\t{0}\n'.format(min_length))
    handle.write('quality_filter\tnone_all_tiers\n')
    handle.write('provirus_rule\tgenomad_name_or_checkv_call_excluded\n')
    handle.write('fasta_headers\ttruncated_to_first_token\n')
    handle.write('subsampling\tnone_all_sequences_read\n')
    for key, value in stats.items():
        handle.write('{0}\t{1}\n'.format(key, value))
    for tier in QUALITY_ORDER:
        handle.write('selected_tier_{0}\t{1}\n'.format(tier, tiers.get(tier, 0)))
    handle.write('fasta_sequences_available\t{0}\n'.format(fasta_total))
    handle.write('fasta_sequences_written\t{0}\n'.format(written))
    handle.write('coverage_table_rows\t{0}\n'.format(coverage_total))
    handle.write('coverage_rows_kept\t{0}\n'.format(coverage_rows))
    handle.write('coverage_sample_columns\t{0}\n'.format(max(0, (columns - 3) // 2)))

total = stats['genomad_virus_sequences']
print('selected {0} of {1} geNomad viral sequences'.format(stats['selected'], total))
print('  excluded: provirus {0}, shorter than {1} bp {2}'.format(
    stats['excluded_provirus'], min_length, stats['excluded_short']))
if stats['selected']:
    print('  FASTA: wrote {0} of {1} selected ({2:.1f}%)'.format(
        written, stats['selected'], 100.0 * written / stats['selected']))
if written:
    print('  coverage join: {0} of {1} written sequences ({2:.1f}%), '
          'from {3} table rows'.format(coverage_rows, written,
                                       100.0 * coverage_rows / written,
                                       coverage_total))
    if coverage_rows < written:
        print('  WARNING: {0} of {1} sequences have no coverage row; vRhyme '
              'bins only what is in both files'.format(written - coverage_rows,
                                                       written),
              file=sys.stderr)
