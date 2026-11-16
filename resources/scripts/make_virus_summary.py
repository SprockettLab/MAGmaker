#!/usr/bin/env python3
"""Merge every viral predictor's calls, with its own CheckV grades, into
output/virus/virus_summary.tsv, plus a name-mapping report.

ONE ROW PER CALL PER PREDICTOR, carrying Predictor and Source_Contig, so
nothing downstream ever has to parse a tool's own naming. With only geNomad
switched on, the rows are what they always were plus those columns.

NAME MAPPING IS THE PART THAT GOES WRONG QUIETLY. Every tool renames
sequences, so each needs its own rule back to the source contig, and a wrong
rule produces a small overlap that reads like genuine disagreement between
tools rather than like a bug. So virus_name_mapping.tsv reports, per
predictor and assembly: sequences seen, how many mapped to a source contig,
how many mapped to a name that is NOT in the assembly, and whether the
predictor's output file was there at all. A predictor that wrote no file is
reported as having none, never as zero calls.

Rules, each verified against real output (WF30, 2026-11-11):

  geNomad       first whitespace token, then cut at '|'. Excised proviruses
                carry |provirus_<start>_<end> on the name itself, so without
                the cut no provirus ever joins.
  VIBRANT       first whitespace token, then a trailing _fragment_N removed.
                VIBRANT appends _fragment_N to the END of the whole header.
                On MEGAHIT headers (k141_7 flag=1 multi=4.0 len=9000) that is
                after len=, so the first token is already the contig. On
                headers with no space (metaSPAdes NODE_..., flye contig_...,
                myloasm) it lands on the contig name itself: without the cut,
                WF18 left 1,528 of 66,546 metaSPAdes, 941 of 2,562 flye and
                1,061 of 2,508 myloasm VIBRANT calls on names not in the
                assembly. The suffix is cut only when what remains is a
                contig of the assembly (or the contig list is unavailable),
                so a real contig name ending in _fragment_N is kept.
  Cenote-Taker3 the first token of the input_name column of its own
                virus_summary.tsv. Its FASTA headers (ct3_99@C3) are renamed
                beyond recovery and are never parsed.
  VirSorter2    cut at '||' (it appends ||full, ||{i}_partial, ||lt2gene),
                then first whitespace token.

A __dup<N> suffix, added by unique_fasta_ids.py so CheckV would accept
VIBRANT's repeated IDs, is stripped first for every tool.

PROVIRUS COLUMNS. Five, not one boolean, because the provenance is what the
induction and prophage-host work needs: Provirus_Predictor (the predicted
name carried provirus coordinates), Provirus_Predictor_Coords (start and end
in the ORIGINAL contig, which locates the integration site), provirus
(CheckV's own call, kept under its old name), Provirus_CheckV_Coords
(CheckV's trimmed region, which is what was actually scored) and
Provirus_Any. CheckV reports No for a sequence geNomad already excised,
because no host DNA is left to detect, so Predictor=Yes with CheckV=No is the
EXPECTED combination for an excised provirus, not a disagreement. Collapsing
these to one flag loses that: on WF30's earlier catalog CheckV alone marked
339 of 3,910 medium-or-better genomes (8.7%) as proviruses, while geNomad OR
CheckV marked 1,157 of 3,910 (29.6%).

ASSEMBLER_CIRCULAR. Topology (geNomad) and CheckV's "Complete" detect a
circle only from a repeat at the two ends of the sequence (direct or inverted
terminal repeats, at least ~20 bp). MEGAHIT leaves about one k-mer of overlap
on a contig that closes a loop, so its circles show up that way. Flye and
myloasm decide circularity in their assembly graph and trim the overlap, so
their circular contigs carry no terminal repeat and read as linear to both
tools: in WF18_FullTest_2026, 6 of 468 flye virus calls on contigs flye marked
circular had Topology DTR, and 2 were CheckV Complete. Assembler_Circular is
the assembler's own call (yes, no, or myloasm's "possibly") from
<sample>.contig_info.tsv, for the call's Source_Contig; it is blank for
short-read assemblers, which report none.

Run by virus.smk's make_virus_summary rule through the `snakemake` object.
"""

import csv
import os
import re
import sys

log_path = str(snakemake.log[0])
os.makedirs(os.path.dirname(log_path) or '.', exist_ok=True)
sys.stderr = open(log_path, 'w')
sys.stdout = sys.stderr

csv.field_size_limit(10 * 1024 * 1024)

contig_samples = list(snakemake.params.contig_samples)
assemblers = list(snakemake.params.assemblers)
genomad_base = snakemake.params.genomad_base
checkv_base = snakemake.params.checkv_base
virus_base = getattr(snakemake.params, 'virus_base', 'output/virus')
predictors = list(getattr(snakemake.params, 'predictors', None) or ['genomad'])
pairs = list(getattr(snakemake.params, 'pairs', []) or
             [(a, r) for a in assemblers for r in contig_samples])
row_sample = dict(getattr(snakemake.params, 'row_sample', {}) or {})
row_stem = dict(getattr(snakemake.params, 'row_stem', {}) or {})
long_read_rows = set(getattr(snakemake.params, 'long_read_rows', []) or [])

DUP = re.compile(r'__dup\d+$')
PROVIRUS = re.compile(r'\|provirus_(\d+)_(\d+)')

CHECKV_FIELDS = ['checkv_quality', 'miuvig_quality', 'checkv_completeness',
                 'checkv_contamination', 'provirus', 'warnings']
CHECKV_BLANK = {field: '' for field in CHECKV_FIELDS}
CHECKV_BLANK['Provirus_CheckV_Coords'] = ''

COLUMNS = [
    'Sample_ID', 'Assembler', 'Binning_Row', 'Predictor', 'Element_Type',
    'Seq_Name', 'Source_Contig', 'Length_BP', 'Topology',
    'Assembler_Circular', 'N_Genes',
    'N_Hallmarks', 'Score', 'FDR', 'Taxonomy', 'Conjugation_Genes',
    'AMR_Genes', 'checkv_quality', 'miuvig_quality', 'checkv_completeness',
    'checkv_contamination', 'provirus', 'warnings', 'Provirus_Predictor',
    'Provirus_Predictor_Coords', 'Provirus_CheckV_Coords', 'Provirus_Any',
]

MAPPING_COLUMNS = [
    'Predictor', 'Sample_ID', 'Assembler', 'Binning_Row', 'Output_Present',
    'Sequences', 'Mapped_To_Contig', 'Unmapped', 'Not_In_Assembly',
    'Example_Not_In_Assembly', 'Source_File',
]


def strip_dup(name):
    """Remove a __dup<N> suffix, which unique_fasta_ids.py appends to the ID
    (the first token), not to the end of the whole header line."""
    tokens = name.split()
    if not tokens:
        return name
    return ' '.join([DUP.sub('', tokens[0])] + tokens[1:])


def first_token(name):
    tokens = strip_dup(name).split()
    return tokens[0] if tokens else ''


def map_genomad(name):
    # Cut at '|': an excised provirus carries |provirus_<start>_<end> on the
    # name itself, so without this no provirus joins its contig.
    return first_token(name).split('|')[0]


def map_vibrant(name):
    # First whitespace token only. _fragment_N is appended after len=, at the
    # end of the header, so splitting on '_' would mangle the contig name.
    return first_token(name)


VIBRANT_FRAGMENT = re.compile(r'_fragment_\d+$')


def vibrant_source(name, contigs):
    """VIBRANT's provirus pieces end in _fragment_N; on a header with no space
    that suffix is part of the first token. Cut it when the name is not a
    contig and the cut name is one (or the contig list is unavailable)."""
    if contigs and name in contigs:
        return name
    base = VIBRANT_FRAGMENT.sub('', name)
    if base != name and (not contigs or base in contigs):
        return base
    return name


def map_virsorter2(name):
    # It appends ||full, ||{i}_partial and ||lt2gene.
    return first_token(name).split('||')[0]


MAPPERS = {'genomad': map_genomad, 'vibrant': map_vibrant,
           'cenotetaker3': map_vibrant, 'virsorter2': map_virsorter2}


def load_tsv(path):
    """Rows of a TSV as dicts; an empty list when it is missing or empty."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    try:
        with open(path, newline='') as handle:
            return list(csv.DictReader(handle, delimiter='\t'))
    except Exception as exc:
        print('Warning: could not read {0}: {1}'.format(path, exc))
        return []


def fasta_lengths(path):
    """name -> sequence length, with the header kept as written."""
    lengths, name = {}, None
    if not path or not os.path.exists(path):
        return lengths
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                name = line[1:].rstrip('\n')
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def assembly_contigs(sample, assembler):
    path = os.path.join('output', 'assemble', assembler,
                        '{0}.contigs.fasta'.format(sample))
    names = set()
    if not os.path.exists(path):
        return names
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                tokens = line[1:].split()
                if tokens:
                    names.add(tokens[0])
    return names


def assembler_circular(sample, assembler):
    """contig -> the assembler's own circular call, from the
    <sample>.contig_info.tsv a long-read assembler writes; empty for an
    assembler that writes none."""
    path = os.path.join('output', 'assemble', assembler,
                        '{0}.contig_info.tsv'.format(sample))
    return {(entry.get('contig') or '').strip(): (entry.get('circular') or '').strip()
            for entry in load_tsv(path)}


def load_checkv(directory):
    """contig_id -> CheckV fields, with its provirus coordinates."""
    graded = {}
    for entry in load_tsv(os.path.join(directory, 'quality_summary.tsv')):
        key = strip_dup((entry.get('contig_id') or '').strip())
        start = (entry.get('proviral_start') or '').strip()
        end = (entry.get('proviral_end') or '').strip()
        graded[key] = {
            'checkv_quality': entry.get('checkv_quality', ''),
            'miuvig_quality': entry.get('miuvig_quality', ''),
            'checkv_completeness': entry.get('completeness', ''),
            'checkv_contamination': entry.get('contamination', ''),
            'provirus': entry.get('provirus', ''),
            'warnings': entry.get('warnings', ''),
            'Provirus_CheckV_Coords': ('{0}-{1}'.format(start, end)
                                       if start and end else ''),
        }
    return graded


def checkv_directory(predictor, assembler, row):
    if predictor == 'genomad':
        return os.path.join(checkv_base, assembler, row)
    return os.path.join(virus_base, 'checkv_' + predictor, assembler, row)


def predictor_fasta(predictor, assembler, row, stem):
    base = os.path.join(virus_base, predictor, assembler, row)
    if predictor == 'genomad':
        return os.path.join(genomad_base, assembler, row,
                            stem + '_summary', stem + '_virus.fna')
    if predictor == 'vibrant':
        return os.path.join(base, 'VIBRANT_' + stem,
                            'VIBRANT_phages_' + stem,
                            stem + '.phages_combined.fna')
    if predictor == 'cenotetaker3':
        return os.path.join(base, 'ct3', 'ct3_virus_sequences.fna')
    return os.path.join(base, 'final-viral-combined.fa')


def provirus_columns(seq_name, checkv):
    """The five provirus columns for one call."""
    match = PROVIRUS.search(seq_name)
    predictor_call = 'Yes' if match else 'No'
    coords = '{0}-{1}'.format(match.group(1), match.group(2)) if match else ''
    checkv_call = (checkv.get('provirus') or '').strip().lower()
    return {
        'Provirus_Predictor': predictor_call,
        'Provirus_Predictor_Coords': coords,
        'Provirus_Any': ('Yes' if match or checkv_call in {'yes', 'true'}
                         else 'No'),
    }


rows, mapping_rows = [], []


def output_present(predictor, assembler, row, path):
    """'no' when the file is missing, and also when the predictor rule left a
    NO_VIRAL_FASTA note: the empty FASTA it writes then is a placeholder, and
    counting it would report "wrote nothing" as zero calls."""
    if not path or not os.path.exists(path):
        return 'no'
    if predictor != 'genomad' and os.path.exists(
            os.path.join(virus_base, predictor, assembler, row, 'NO_VIRAL_FASTA')):
        return 'no'
    return 'yes'


def note_mapping(predictor, sample, assembler, row, path, sequences,
                 mapped, unmapped, not_in_assembly):
    mapping_rows.append({
        'Predictor': predictor, 'Sample_ID': sample, 'Assembler': assembler,
        'Binning_Row': row,
        'Output_Present': output_present(predictor, assembler, row, path),
        'Sequences': sequences, 'Mapped_To_Contig': mapped,
        'Unmapped': unmapped, 'Not_In_Assembly': len(not_in_assembly),
        'Example_Not_In_Assembly': ','.join(sorted(not_in_assembly)[:3]),
        'Source_File': path or '',
    })


def genomad_rows(sample, assembler, row, stem, contigs):
    """geNomad's own summary tables: the richest columns of any predictor."""
    summary_dir = os.path.join(genomad_base, assembler, row, stem + '_summary')
    virus_table = os.path.join(summary_dir, stem + '_virus_summary.tsv')
    plasmid_table = os.path.join(summary_dir, stem + '_plasmid_summary.tsv')
    graded = load_checkv(checkv_directory('genomad', assembler, row))

    mapped, unmapped, missing = 0, 0, set()
    for entry in load_tsv(virus_table):
        seq_name = str(entry.get('seq_name', ''))
        source = map_genomad(seq_name)
        if source:
            mapped += 1
            if contigs and source not in contigs:
                missing.add(source)
        else:
            unmapped += 1
        key = strip_dup(seq_name.split()[0]) if seq_name.split() else ''
        checkv = graded.get(key, graded.get(seq_name, CHECKV_BLANK))
        record = {
            'Sample_ID': sample, 'Assembler': assembler, 'Binning_Row': row,
            'Predictor': 'genomad', 'Element_Type': 'virus',
            'Seq_Name': seq_name, 'Source_Contig': source,
            'Length_BP': entry.get('length', ''),
            'Topology': entry.get('topology', ''),
            'N_Genes': entry.get('n_genes', ''),
            'N_Hallmarks': entry.get('n_hallmarks', ''),
            'Score': entry.get('virus_score', ''),
            'FDR': entry.get('fdr', ''),
            'Taxonomy': entry.get('taxonomy', ''),
            'Conjugation_Genes': '', 'AMR_Genes': '',
        }
        record.update(checkv)
        record.update(provirus_columns(seq_name, checkv))
        rows.append(record)

    # Plasmids. geNomad is the only one of the four that calls them, and they
    # come out of the same run, which is why geNomad is one switch for both.
    for entry in load_tsv(plasmid_table):
        seq_name = str(entry.get('seq_name', ''))
        source = map_genomad(seq_name)
        if source:
            mapped += 1
            if contigs and source not in contigs:
                missing.add(source)
        else:
            unmapped += 1
        record = {
            'Sample_ID': sample, 'Assembler': assembler, 'Binning_Row': row,
            'Predictor': 'genomad', 'Element_Type': 'plasmid',
            'Seq_Name': seq_name, 'Source_Contig': source,
            'Length_BP': entry.get('length', ''),
            'Topology': entry.get('topology', ''),
            'N_Genes': entry.get('n_genes', ''),
            'N_Hallmarks': entry.get('n_hallmarks', ''),
            'Score': entry.get('plasmid_score', ''),
            'FDR': entry.get('fdr', ''), 'Taxonomy': '',
            'Conjugation_Genes': entry.get('conjugation_genes', ''),
            'AMR_Genes': entry.get('amr_genes', ''),
        }
        record.update(CHECKV_BLANK)
        record.update(provirus_columns(seq_name, CHECKV_BLANK))
        rows.append(record)

    note_mapping('genomad', sample, assembler, row, virus_table,
                 mapped + unmapped, mapped, unmapped, missing)


def ct3_input_names(assembler, row):
    """Cenote-Taker3's own name -> the original header it came from."""
    path = os.path.join(virus_base, 'cenotetaker3', assembler, row, 'ct3',
                        'ct3_virus_summary.tsv')
    # Cenote-Taker3 3.4.4 sometimes lists a contig a second time as a blank
    # row: the name with a trailing space, no input_name, length 0 (WF18
    # metaSPAdes, e.g. "ct3_33 " right after ct3_33, all DTR genomes). The
    # strip makes the two names equal, so a blank row may never replace a
    # filled one: that silently unmapped 71 of 80,930 calls, 70 of them
    # medium quality or better.
    names = {}
    for entry in load_tsv(path):
        own = (entry.get('contig') or '').strip()
        original = (entry.get('input_name') or '').strip()
        if own and (original or own not in names):
            names[own] = original
    return names, path


def other_rows(predictor, sample, assembler, row, stem, contigs):
    """A predictor whose own tables we do not read: names come from its FASTA
    (through its mapping rule) and quality from its own CheckV."""
    fasta = predictor_fasta(predictor, assembler, row, stem)
    lengths = fasta_lengths(fasta)
    graded = load_checkv(checkv_directory(predictor, assembler, row))
    ct3_names, ct3_path = ({}, None)
    if predictor == 'cenotetaker3':
        ct3_names, ct3_path = ct3_input_names(assembler, row)

    mapped, unmapped, missing = 0, 0, set()
    for header, length in lengths.items():
        if predictor == 'cenotetaker3':
            # Its FASTA headers are unrecoverable, so the source comes from
            # input_name in its own summary table.
            own = strip_dup(header.split()[0]) if header.split() else ''
            source = map_vibrant(ct3_names.get(own, ''))
        else:
            source = MAPPERS[predictor](header)
            if predictor == 'vibrant':
                source = vibrant_source(source, contigs)
        if source:
            mapped += 1
            if contigs and source not in contigs:
                missing.add(source)
        else:
            unmapped += 1
        key = strip_dup(header.split()[0]) if header.split() else ''
        checkv = graded.get(key, graded.get(header, CHECKV_BLANK))
        record = {
            'Sample_ID': sample, 'Assembler': assembler, 'Binning_Row': row,
            'Predictor': predictor, 'Element_Type': 'virus',
            'Seq_Name': header, 'Source_Contig': source,
            'Length_BP': length, 'Topology': '', 'N_Genes': '',
            'N_Hallmarks': '', 'Score': '', 'FDR': '', 'Taxonomy': '',
            'Conjugation_Genes': '', 'AMR_Genes': '',
        }
        record.update(checkv)
        record.update(provirus_columns(header, checkv))
        rows.append(record)

    note_mapping(predictor, sample, assembler, row,
                 ct3_path if predictor == 'cenotetaker3' else fasta,
                 len(lengths), mapped, unmapped, missing)


for assembler, row in pairs:
    sample = row_sample.get(row, row)
    stem = row_stem.get(row, sample + '.contigs')
    contigs = assembly_contigs(sample, assembler)
    # A short-read row has one assembler that actually produced its contigs;
    # a long-read row has exactly one, its own.
    if row not in long_read_rows and not contigs:
        continue
    circular = assembler_circular(sample, assembler) if row in long_read_rows else {}
    first_new = len(rows)
    for predictor in predictors:
        if predictor == 'genomad':
            genomad_rows(sample, assembler, row, stem, contigs)
        else:
            other_rows(predictor, sample, assembler, row, stem, contigs)
    for record in rows[first_new:]:
        record['Assembler_Circular'] = circular.get(record['Source_Contig'], '')
    if row in long_read_rows:
        found = sum(1 for r in rows[first_new:] if r['Assembler_Circular'])
        print('{0} ({1}): assembler circular call for {2} of {3} rows'.format(
            row, assembler, found, len(rows) - first_new))


def write(path, columns, table):
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        for record in table:
            writer.writerow({column: record.get(column, '')
                             for column in columns})


write(snakemake.output.summary, COLUMNS, rows)
write(snakemake.output.mapping, MAPPING_COLUMNS, mapping_rows)

print('predictors: {0}'.format(', '.join(predictors)))
for predictor in predictors:
    own = [r for r in rows if r['Predictor'] == predictor]
    virus = sum(1 for r in own if r['Element_Type'] == 'virus')
    plasmid = sum(1 for r in own if r['Element_Type'] == 'plasmid')
    print('{0}: {1} viral + {2} plasmid calls'.format(predictor, virus, plasmid))

for entry in mapping_rows:
    if entry['Output_Present'] == 'no':
        print('NOTE: {0} has no output file for {1} ({2}): {3}'.format(
            entry['Predictor'], entry['Binning_Row'], entry['Assembler'],
            entry['Source_File']), file=sys.stderr)
    if entry['Unmapped']:
        print('WARNING: {0} {1}: {2} of {3} sequences did not map to a source '
              'contig'.format(entry['Predictor'], entry['Binning_Row'],
                              entry['Unmapped'], entry['Sequences']),
              file=sys.stderr)
    if entry['Not_In_Assembly']:
        print('WARNING: {0} {1}: {2} mapped names are not in the assembly '
              '({3}); the mapping rule is probably wrong'.format(
                  entry['Predictor'], entry['Binning_Row'],
                  entry['Not_In_Assembly'], entry['Example_Not_In_Assembly']),
              file=sys.stderr)

print('wrote {0} rows to {1}'.format(len(rows), snakemake.output.summary))
print('wrote {0} rows to {1}'.format(len(mapping_rows), snakemake.output.mapping))
