import os
import sys
import glob
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
assemblers = snakemake.params.assemblers


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
    tax_map = {}
    for db_type in ('bac120', 'ar53'):
        summary = os.path.join(gtdbtk_dir, f'gtdbtk.{db_type}.summary.tsv')
        if os.path.exists(summary):
            df = pd.read_csv(summary, sep='\t')
            for _, row in df.iterrows():
                genome = str(row.get('user_genome', ''))
                classification = str(row.get('classification', ''))
                tax_map[genome] = (parse_gtdbtk_classification(classification), classification)
    return tax_map


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


def resolve_binner(bin_name, bin_to_binner):
    """Look up binner for a bin, falling back to the parent name for DAS_Tool _sub bins."""
    if bin_name in bin_to_binner:
        return bin_to_binner[bin_name]
    if bin_name.endswith('_sub'):
        parent = bin_name[:-4]
        if parent in bin_to_binner:
            return bin_to_binner[parent]
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


# Collect all MAGs
all_mags = []

for mapper in mappers:
    for sample in contig_samples:
        bins_dir = os.path.join(bins_base, mapper, 'DAS_Tool_Fastas', sample)
        if not os.path.isdir(bins_dir):
            print(f"Warning: bins directory not found: {bins_dir}")
            continue

        bin_to_binner = build_binner_map(bins_base, mapper, sample)
        tax_map = load_gtdbtk(os.path.join(gtdbtk_base, mapper, sample))
        checkm2_map = load_checkm2(os.path.join(checkm2_base, mapper, sample))
        gunc_map = load_gunc(os.path.join(gunc_base, mapper, sample))
        assembler = get_assembler(sample, assemblers)

        for fa in sorted(glob.glob(os.path.join(bins_dir, '*.fa'))):
            bin_name = os.path.splitext(os.path.basename(fa))[0]
            binner = resolve_binner(bin_name, bin_to_binner)

            tax_entry = tax_map.get(bin_name, (parse_gtdbtk_classification(''), ''))
            tax_dict, full_classification = tax_entry
            tax_label = get_taxonomic_label(tax_dict)

            stats = compute_mag_stats(fa)
            qc = checkm2_map.get(bin_name, {
                'completeness': '', 'contamination': '', 'quality_score': '',
                'coding_density': '', 'total_coding_sequences': ''
            })
            gunc = gunc_map.get(bin_name, {
                'gunc_clade_separation_score': '', 'gunc_pass': ''
            })

            all_mags.append({
                'original_name': bin_name,
                'original_path': fa,
                'sample_id': sample,
                'assembler': assembler,
                'winning_binner': binner,
                'tax_dict': tax_dict,
                'tax_label': tax_label,
                'gtdbtk_classification': full_classification,
                **qc,
                **gunc,
                **stats
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
        'Completeness': mag.get('completeness', ''),
        'Contamination': mag.get('contamination', ''),
        'Quality_Score': mag.get('quality_score', ''),
        'MIMAG_Quality': mimag_quality(mag.get('completeness', ''), mag.get('contamination', '')),
        'GUNC_Clade_Separation_Score': mag.get('gunc_clade_separation_score', ''),
        'GUNC_Pass': mag.get('gunc_pass', ''),
        'Total_Length_BP': mag.get('total_length_bp', ''),
        'Num_Contigs': mag.get('num_contigs', ''),
        'Largest_Contig': mag.get('largest_contig', ''),
        'GC_Percent': mag.get('gc_percent', ''),
        'N50': mag.get('N50', ''),
        'Coding_Density': mag.get('coding_density', ''),
        'Total_Coding_Sequences': mag.get('total_coding_sequences', ''),
    })

pd.DataFrame(rows).to_csv(snakemake.output.summary, sep='\t', index=False)
print(f"Wrote summary for {len(rows)} MAGs to {snakemake.output.summary}")
