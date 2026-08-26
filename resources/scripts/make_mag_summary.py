import os
import sys
import glob
import gzip
import pandas as pd

log_path = str(snakemake.log[0])
sys.stderr = open(log_path, 'w')
sys.stdout = sys.stderr

mappers = snakemake.params.mappers
contig_samples = snakemake.params.contig_samples
bins_base = snakemake.params.bins_base
gtdbtk_base = snakemake.params.gtdbtk_base
checkm2_base = snakemake.params.checkm2_base
gunc_base = snakemake.params.gunc_base
cmseq_base = snakemake.params.cmseq_base
assemblers = snakemake.params.assemblers
# 'das_tool' or 'binette', and the directory that tool's selected bins are
# in. Defaulted so an older config that predates the switch still works.
consolidation_tool = getattr(snakemake.params, 'consolidation_tool', 'das_tool')
selected_fastas = getattr(snakemake.params, 'selected_fastas', 'DAS_Tool_Fastas')


def compute_mag_stats(fasta_path):
    lengths = []
    gc_count = 0
    total_bases = 0
    current_seq = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    seq = ''.join(current_seq)
                    lengths.append(len(seq))
                    gc_count += seq.upper().count('G') + seq.upper().count('C')
                    total_bases += len(seq)
                    current_seq = []
            else:
                current_seq.append(line)
        if current_seq:
            seq = ''.join(current_seq)
            lengths.append(len(seq))
            gc_count += seq.upper().count('G') + seq.upper().count('C')
            total_bases += len(seq)
    lengths.sort(reverse=True)
    total_length = sum(lengths)
    n50 = 0
    cumsum = 0
    for length in lengths:
        cumsum += length
        if cumsum >= total_length / 2:
            n50 = length
            break
    return {
        'total_length_bp': total_length,
        'num_contigs': len(lengths),
        'largest_contig': lengths[0] if lengths else 0,
        'gc_percent': round(gc_count / total_bases * 100, 2) if total_bases > 0 else 0,
        'N50': n50
    }


def parse_gtdbtk_classification(classification_str):
    rank_keys = {'d': 'domain', 'p': 'phylum', 'c': 'class',
                 'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species'}
    tax = {v: '' for v in rank_keys.values()}
    if not classification_str or str(classification_str) in ('N/A', 'nan', ''):
        return tax
    for part in str(classification_str).split(';'):
        part = part.strip()
        if '__' in part:
            prefix, value = part.split('__', 1)
            if prefix in rank_keys:
                tax[rank_keys[prefix]] = value.strip()
    return tax


def get_taxonomic_label(tax):
    for rank in ('species', 'genus', 'family', 'order', 'class', 'phylum', 'domain'):
        value = tax.get(rank, '').strip()
        if value:
            return value.replace(' ', '_')
    return 'unclassified'


def load_gtdbtk(gtdbtk_dir):
    """Taxonomy plus the two placement fields worth keeping.

    msa_percent is the share of the concatenated marker alignment the
    genome fills with residues rather than gaps. That is what decides how
    well a genome can be placed in a phylogeny, and it is not the same as
    completeness: a genome can recover most of its markers as fragments
    and still fill little of the alignment. GTDB-Tk already computes it,
    so surfacing it here saves every downstream analysis from re-deriving
    it from output that later runs overwrite.

    warnings carries GTDB-Tk's own caveats about a placement, which are
    otherwise lost once the summary is overwritten."""
    tax_map = {}
    for db_type in ('bac120', 'ar53'):
        summary = os.path.join(gtdbtk_dir, f'gtdbtk.{db_type}.summary.tsv')
        if os.path.exists(summary):
            df = pd.read_csv(summary, sep='\t')
            for _, row in df.iterrows():
                genome = str(row.get('user_genome', ''))
                classification = str(row.get('classification', ''))
                msa = row.get('msa_percent', '')
                warn = row.get('warnings', '')
                tax_map[genome] = (
                    parse_gtdbtk_classification(classification),
                    classification,
                    '' if pd.isna(msa) else msa,
                    '' if pd.isna(warn) else str(warn),
                )
    return tax_map


def load_msa_percent(gtdbtk_dir):
    """Recover msa_percent from the alignment `align` writes.

    in : a sample's gtdbtk output directory
    out: {genome: msa_percent as a float}, empty if no alignment is there

    GTDB-Tk defines the column as "the percentage of the MSA spanned by the
    genome (i.e. percentage of columns with an amino acid)", and `align`
    writes exactly that alignment as gtdbtk.<domain>.user_msa.fasta[.gz].
    Counting non-gap columns therefore reproduces GTDB-Tk's own number
    rather than approximating it.

    This exists so MSA_Percent survives --skip-gtdbtk-classify. classify is
    where pplacer and the ~140 GB memory requirement live, and it is also
    what writes summary.tsv, so without this the one column that decides
    whether a genome can be placed in a tree would be the one thing lost by
    skipping the step that places it.

    A genome present in the bin set but absent from the alignment is not
    recorded here. Whether that means "filtered by --min_perc_aa" or "no
    markers found" depends on the run, and guessing between them in a
    column that reads as a measurement would be worse than leaving it
    empty.
    """
    msa_map = {}
    align_dir = os.path.join(gtdbtk_dir, 'align')
    if not os.path.isdir(align_dir):
        return msa_map

    for db_type in ('bac120', 'ar53'):
        for suffix, opener in (('.gz', gzip.open), ('', open)):
            path = os.path.join(
                align_dir, f'gtdbtk.{db_type}.user_msa.fasta{suffix}')
            if not os.path.exists(path):
                continue
            try:
                with opener(path, 'rt') as fh:
                    name, seq = None, []
                    for line in fh:
                        line = line.strip()
                        if line.startswith('>'):
                            if name is not None:
                                msa_map[name] = _percent_aligned(''.join(seq))
                            # GTDB-Tk writes the genome id as the whole
                            # header, but split defensively in case a
                            # description is ever appended.
                            name = line[1:].split()[0] if len(line) > 1 else ''
                            seq = []
                        elif line:
                            seq.append(line)
                    if name:
                        msa_map[name] = _percent_aligned(''.join(seq))
            except Exception as e:
                print(f"Warning: could not read {path}: {e}")
            break   # one of .gz / plain per domain, whichever exists
    return msa_map


def _percent_aligned(sequence):
    """Share of alignment columns this genome fills with a residue."""
    if not sequence:
        return ''
    gaps = sequence.count('-') + sequence.count('.')
    return round(100.0 * (len(sequence) - gaps) / len(sequence), 2)


def load_cmseq(cmseq_dir):
    """in : a sample's cmseq directory
       out: {MAG: (strain heterogeneity, positions evaluated)}

    Missing directory means CMSeq did not run for this sample, which the
    caller renders as NA. That is deliberately different from a MAG CMSeq
    evaluated and found too poorly covered to score, which is also NA but
    carries a position count."""
    out = {}
    path = os.path.join(cmseq_dir, 'strain_heterogeneity.tsv')
    if os.path.exists(path):
        df = pd.read_csv(path, sep='\t', dtype=str).fillna('NA')
        for _, row in df.iterrows():
            out[str(row.get('MAG', ''))] = (
                str(row.get('Strain_Heterogeneity', 'NA')),
                str(row.get('SH_Positions_Evaluated', 'NA')),
            )
    return out


def load_checkm2(checkm2_dir):
    checkm2_map = {}
    report = os.path.join(checkm2_dir, 'quality_report.tsv')
    if not os.path.exists(report) or os.path.getsize(report) == 0:
        return checkm2_map
    df = pd.read_csv(report, sep='\t')
    if df.empty:
        return checkm2_map
    for _, row in df.iterrows():
        name = str(row.get('Name', ''))
        completeness = row.get('Completeness', None)
        contamination = row.get('Contamination', None)
        try:
            quality_score = round(float(completeness) - 5 * float(contamination), 2)
        except (TypeError, ValueError):
            quality_score = ''
        checkm2_map[name] = {
            'completeness': completeness,
            'contamination': contamination,
            'quality_score': quality_score,
            'coding_density': row.get('Coding_Density', ''),
            'total_coding_sequences': row.get('Total_Coding_Sequences', ''),
        }
    return checkm2_map


def load_gunc(gunc_dir):
    gunc_map = {}
    matches = glob.glob(os.path.join(gunc_dir, '*.maxCSS_level.tsv'))
    if matches:
        df = pd.read_csv(matches[0], sep='\t')
        for _, row in df.iterrows():
            name = str(row.get('genome', ''))
            gunc_map[name] = {
                'gunc_clade_separation_score': row.get('clade_separation_score', ''),
                'gunc_pass': row.get('pass.GUNC', '')
            }
    return gunc_map


def build_binner_map(bins_base, mapper, sample, binners=('metabat2', 'maxbin2', 'concoct')):
    """Map each winning bin name to the binner that produced it via scaffolds2bin files.
    DAS_Tool _sub bins are refined sub-bins of a parent bin — strip the suffix to find the parent."""
    bin_to_binner = {}
    for binner in binners:
        s2b = os.path.join(bins_base, binner, mapper, 'scaffolds2bin',
                           f'{sample}_scaffolds2bin.tsv')
        if not os.path.exists(s2b) or os.path.getsize(s2b) == 0:
            continue
        try:
            df = pd.read_csv(s2b, sep='\t', header=None, names=['contig', 'bin'])
            for bin_name in df['bin'].unique():
                bin_to_binner[str(bin_name)] = binner
        except Exception as e:
            print(f"Warning: could not read {s2b}: {e}")
    return bin_to_binner


def resolve_binner(bin_name, bin_to_binner,
                   binners=('metabat2', 'maxbin2', 'concoct')):
    """Look up the binner that produced a winning DAS_Tool bin.

    Primary: match the bin name against the scaffolds2bin map, falling back to
    the parent name for DAS_Tool `_sub` refined sub-bins.

    Fallback: run_DAS_Tool is called with `--labels metabat2,maxbin2,concoct`,
    so winning bins are prefixed with their source binner (e.g. `concoct_13`,
    `maxbin2_v2489_Mother_Day7_bin.088`). When the scaffolds2bin lookup does not
    resolve — e.g. a sample whose scaffolds2bin files are missing or empty, which
    otherwise sent every bin to 'unknown' — recover the binner from that prefix.
    """
    if bin_name in bin_to_binner:
        return bin_to_binner[bin_name]
    if bin_name.endswith('_sub'):
        parent = bin_name[:-4]
        if parent in bin_to_binner:
            return bin_to_binner[parent]
    # Recover from the DAS_Tool label prefix (labels are the binner names).
    for binner in binners:
        if bin_name.startswith(binner + '_'):
            return binner
    return 'unknown'


def build_binette_origin_map(bins_base, mapper, sample):
    """Map each Binette bin name to the bin set(s) it was built from.

    in : the selected_bins base, mapper and sample
    out: {bin_name: origin string}, empty if the report is missing

    Binette's final_bins_quality_reports.tsv carries an `origin` column
    holding the names of the input bin sets a bin derives from, semicolon
    separated. stage_binette_inputs names those sets after the binners, so
    the values come back as metabat2/maxbin2/concoct.
    """
    report = os.path.join(bins_base, mapper, 'run_binette', sample,
                          'final_bins_quality_reports.tsv')
    if not os.path.exists(report) or os.path.getsize(report) == 0:
        return {}
    try:
        df = pd.read_csv(report, sep='\t', dtype=str)
    except Exception as e:
        print(f"Warning: could not read {report}: {e}")
        return {}
    if 'name' not in df.columns or 'origin' not in df.columns:
        print(f"Warning: {report} lacks name/origin columns")
        return {}
    return {str(n): str(o) for n, o in zip(df['name'], df['origin'])}


def resolve_binette_binner(bin_name, origin_map,
                           binners=('metabat2', 'maxbin2', 'concoct')):
    """Turn a Binette origin into a Winning_Binner value.

    Unlike DAS_Tool, Binette can return a bin that no single binner
    produced: a hybrid built from the intersection, difference or union of
    bins from two binners. Reporting one binner for those would be wrong,
    so contributors are joined with '+' (e.g. `metabat2+concoct`), which
    stays greppable while being visibly not a single binner.

    Binette labels a bin it constructed with no surviving input origin as
    'binette'; that is passed through rather than called unknown, since it
    is a real answer about where the bin came from.
    """
    origin = origin_map.get(bin_name, '')
    if not origin:
        return 'unknown'
    found = [b for b in binners if b in origin]
    if found:
        return '+'.join(found)
    if 'binette' in origin:
        return 'binette'
    return 'unknown'


def get_assembler(sample, assemblers):
    """Return whichever configured assembler produced contigs for this sample."""
    for assembler in assemblers:
        contigs = os.path.join('output', 'assemble', assembler, f'{sample}.contigs.fasta')
        if os.path.exists(contigs):
            return assembler
    return 'unknown'


def mimag_quality(completeness, contamination):
    try:
        c = float(completeness)
        x = float(contamination)
    except (TypeError, ValueError):
        return ''
    if c >= 90 and x < 5:
        return 'HQ'
    elif c >= 50 and x < 10:
        return 'MQ'
    else:
        return 'LQ'


# Column order of mag_summary.tsv. Defined once because two kinds of row are
# written into it: one per MAG, and one per sample that produced none.
SUMMARY_COLUMNS = (
    'MAG_ID', 'New_Name', 'Original_Name', 'Original_Path', 'Sample_ID',
    'Assembler', 'Winning_Binner', 'Domain', 'Phylum', 'Class', 'Order',
    'Family', 'Genus', 'Species', 'GTDB_Classification', 'MSA_Percent',
    'GTDBTk_Warnings', 'Completeness', 'Contamination', 'Quality_Score',
    'MIMAG_Quality', 'GUNC_Clade_Separation_Score', 'GUNC_Pass',
    'Strain_Heterogeneity', 'SH_Positions_Evaluated', 'Total_Length_BP',
    'Num_Contigs', 'Largest_Contig', 'GC_Percent', 'N50', 'Coding_Density',
    'Total_Coding_Sequences', 'Notes',
)


def binners_with_bins(bins_base, mapper, sample):
    """
    in : selected_bins base, mapper, sample
    out: (contributed, declined) lists of binner names

    A binner's scaffolds2bin table exists for every sample it was asked
    about; an empty one means it declined that sample. Reading the tables
    rather than the logs keeps this independent of any binner's messages.
    """
    contributed, declined = [], []
    pattern = os.path.join(bins_base, '*', mapper, 'scaffolds2bin',
                           f'{sample}_scaffolds2bin.tsv')
    for path in sorted(glob.glob(pattern)):
        binner = path.split(os.sep)[-4]
        try:
            has_bins = os.path.getsize(path) > 0
        except OSError:
            continue
        (contributed if has_bins else declined).append(binner)
    return contributed, declined


# Collect all MAGs
all_mags = []

# Samples that finished but yielded no MAGs. Without these the summary is
# silent about them, which makes a sample nothing could bin look identical
# to one that was never run.
samples_without_mags = []

for mapper in mappers:
    for sample in contig_samples:
        bins_dir = os.path.join(bins_base, mapper, selected_fastas, sample)
        if not os.path.isdir(bins_dir):
            print(f"Warning: bins directory not found: {bins_dir}")
            samples_without_mags.append({
                'sample': sample,
                'mapper': mapper,
                'note': 'no bins directory; consolidation did not run',
            })
            continue

        if consolidation_tool == 'binette':
            bin_to_binner = {}
            origin_map = build_binette_origin_map(bins_base, mapper, sample)
        else:
            bin_to_binner = build_binner_map(bins_base, mapper, sample)
            origin_map = {}
        tax_map = load_gtdbtk(os.path.join(gtdbtk_base, mapper, sample))
        # Fallback for MSA_Percent when classify did not run and so wrote no
        # summary.tsv to read it from.
        msa_map = load_msa_percent(os.path.join(gtdbtk_base, mapper, sample))
        cmseq_map = load_cmseq(os.path.join(cmseq_base, mapper, sample))
        checkm2_map = load_checkm2(os.path.join(checkm2_base, mapper, sample))
        gunc_map = load_gunc(os.path.join(gunc_base, mapper, sample))
        assembler = get_assembler(sample, assemblers)

        n_before = len(all_mags)
        missing_tax = []
        missing_checkm2 = []
        missing_gunc = []

        for fa in sorted(glob.glob(os.path.join(bins_dir, '*.fa'))):
            bin_name = os.path.splitext(os.path.basename(fa))[0]
            if consolidation_tool == 'binette':
                binner = resolve_binette_binner(bin_name, origin_map)
            else:
                binner = resolve_binner(bin_name, bin_to_binner)

            tax_entry = tax_map.get(
                bin_name, (parse_gtdbtk_classification(''), '', '', ''))
            if bin_name not in tax_map:
                missing_tax.append(bin_name)
            tax_dict, full_classification, msa_percent, gtdbtk_warnings = tax_entry
            if msa_percent == '':
                msa_percent = msa_map.get(bin_name, '')
            tax_label = get_taxonomic_label(tax_dict)

            stats = compute_mag_stats(fa)
            qc = checkm2_map.get(bin_name, {
                'completeness': '', 'contamination': '', 'quality_score': '',
                'coding_density': '', 'total_coding_sequences': ''
            })
            if bin_name not in checkm2_map:
                missing_checkm2.append(bin_name)
            gunc = gunc_map.get(bin_name, {
                'gunc_clade_separation_score': '', 'gunc_pass': ''
            })
            if bin_name not in gunc_map:
                missing_gunc.append(bin_name)
            sh, sh_n = cmseq_map.get(bin_name, ('NA', 'NA'))

            all_mags.append({
                'original_name': bin_name,
                'original_path': fa,
                'sample_id': sample,
                'assembler': assembler,
                'winning_binner': binner,
                'tax_dict': tax_dict,
                'tax_label': tax_label,
                'gtdbtk_classification': full_classification,
                'strain_heterogeneity': sh,
                'sh_positions_evaluated': sh_n,
                'msa_percent': msa_percent,
                'gtdbtk_warnings': gtdbtk_warnings,
                **qc,
                **gunc,
                **stats
            })

        # A stage that produced no result for a bin is written to the summary
        # as an empty field, which is indistinguishable from a stage that ran
        # and had nothing to say. That ambiguity hid a GTDB-Tk failure across
        # four samples for four months: bac120 classify died, every bacterial
        # MAG got a blank taxonomy, and nothing in the output said so (#86).
        # Report the gaps here, while the sample that caused them is still
        # named, rather than leaving them to be noticed downstream.
        n_bins_here = len(all_mags) - n_before
        for stage, missing in (('GTDB-Tk', missing_tax),
                               ('CheckM2', missing_checkm2),
                               ('GUNC', missing_gunc)):
            if not missing:
                continue
            scope = 'ALL' if len(missing) == n_bins_here else str(len(missing))
            print(
                f"WARNING: {mapper}/{sample}: no {stage} result for {scope} of "
                f"{n_bins_here} bins (e.g. {', '.join(missing[:3])}). Those "
                f"fields are blank in mag_summary.tsv and are NOT a result.",
                file=sys.stderr,
            )

        # Record why, distinguishing "every binner declined this sample"
        # from "binners produced bins but none survived selection". Both
        # are real outcomes and they mean different things.
        if len(all_mags) == n_before:
            contributed, declined = binners_with_bins(bins_base, mapper, sample)
            if contributed:
                note = ('no MAGs passed selection; bins were offered by: '
                        + ', '.join(contributed))
            elif declined:
                note = ('no binner produced bins; declined by: '
                        + ', '.join(declined))
            else:
                note = 'no binner output found for this sample'
            samples_without_mags.append({
                'sample': sample, 'mapper': mapper, 'note': note,
            })

# Sort by taxonomy then MIMAG quality (HQ first) then original name
TAX_RANKS = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
MIMAG_ORDER = {'HQ': 0, 'MQ': 1, 'LQ': 2, '': 3}
all_mags.sort(
    key=lambda x: (
        tuple(x['tax_dict'].get(r, '\xff') or '\xff' for r in TAX_RANKS),
        MIMAG_ORDER.get(mimag_quality(x.get('completeness', ''), x.get('contamination', '')), 3),
        x['original_name'],
    )
)

rows = []
for i, mag in enumerate(all_mags, 1):
    mag_id = f"MAG_{i:04d}"
    tax = mag['tax_dict']
    rows.append({
        'MAG_ID': mag_id,
        'New_Name': f"{mag_id}__{mag['tax_label']}",
        'Original_Name': mag['original_name'],
        'Original_Path': mag['original_path'],
        'Sample_ID': mag['sample_id'],
        'Assembler': mag['assembler'],
        'Winning_Binner': mag['winning_binner'],
        'Domain': tax.get('domain', '') or 'NA',
        'Phylum': tax.get('phylum', '') or 'NA',
        'Class': tax.get('class', '') or 'NA',
        'Order': tax.get('order', '') or 'NA',
        'Family': tax.get('family', '') or 'NA',
        'Genus': tax.get('genus', '') or 'NA',
        'Species': tax.get('species', '') or 'NA',
        'GTDB_Classification': mag['gtdbtk_classification'] or 'NA',
        'MSA_Percent': mag.get('msa_percent', '') if str(
            mag.get('msa_percent', '')) != '' else 'NA',
        'GTDBTk_Warnings': mag.get('gtdbtk_warnings', '') or 'NA',
        'Completeness': mag.get('completeness', ''),
        'Contamination': mag.get('contamination', ''),
        'Quality_Score': mag.get('quality_score', ''),
        'MIMAG_Quality': mimag_quality(mag.get('completeness', ''), mag.get('contamination', '')),
        'GUNC_Clade_Separation_Score': mag.get('gunc_clade_separation_score', ''),
        'GUNC_Pass': mag.get('gunc_pass', ''),
        'Strain_Heterogeneity': mag.get('strain_heterogeneity', 'NA') or 'NA',
        'SH_Positions_Evaluated': mag.get('sh_positions_evaluated', 'NA') or 'NA',
        'Total_Length_BP': mag.get('total_length_bp', ''),
        'Num_Contigs': mag.get('num_contigs', ''),
        'Largest_Contig': mag.get('largest_contig', ''),
        'GC_Percent': mag.get('gc_percent', ''),
        'N50': mag.get('N50', ''),
        'Coding_Density': mag.get('coding_density', ''),
        'Total_Coding_Sequences': mag.get('total_coding_sequences', ''),
        'Notes': 'NA',
    })

# One row per sample that produced no MAGs, after the numbered MAGs. MAG_ID
# is NONE rather than a number so these never look like MAGs, and
# Original_Path is NA so rename_mags skips them.
for entry in sorted(samples_without_mags, key=lambda e: (e['sample'], e['mapper'])):
    row = dict.fromkeys(SUMMARY_COLUMNS, 'NA')
    row['MAG_ID'] = 'NONE'
    row['New_Name'] = 'NONE'
    row['Sample_ID'] = entry['sample']
    row['Assembler'] = get_assembler(entry['sample'], assemblers)
    row['Notes'] = entry['note']
    rows.append(row)

# Explicit columns so the placeholder rows and the MAG rows cannot drift
# apart, and so column order is stable across runs.
pd.DataFrame(rows, columns=list(SUMMARY_COLUMNS)).to_csv(
    snakemake.output.summary, sep='\t', index=False)

n_mags = len(all_mags)
print(f"Wrote summary for {n_mags} MAGs to {snakemake.output.summary}")
if samples_without_mags:
    print(f"{len(samples_without_mags)} sample(s) produced no MAGs:")
    for entry in samples_without_mags:
        print(f"  {entry['sample']} ({entry['mapper']}): {entry['note']}")
