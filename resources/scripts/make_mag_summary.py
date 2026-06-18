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
    for l in lengths:
        cumsum += l
        if cumsum >= total_length / 2:
            n50 = l
            break
    return {
        'total_length_bp': total_length,
        'num_contigs': len(lengths),
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
    if os.path.exists(report):
        df = pd.read_csv(report, sep='\t')
        for _, row in df.iterrows():
            name = str(row.get('Name', ''))
            checkm2_map[name] = {
                'completeness': row.get('Completeness', ''),
                'contamination': row.get('Contamination', ''),
                'quality_score': row.get('Completeness_General', row.get('quality_score', ''))
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


def load_dastool_summary(summary_path):
    binner_map = {}
    if os.path.exists(summary_path):
        df = pd.read_csv(summary_path, sep='\t')
        id_col = df.columns[0]
        for _, row in df.iterrows():
            bin_id = str(row[id_col])
            binner_map[bin_id] = str(row.get('tool_used', bin_id.split('.')[0]))
    return binner_map


# Collect all MAGs
all_mags = []

for mapper in mappers:
    for sample in contig_samples:
        bins_dir = os.path.join(bins_base, mapper, 'DAS_Tool_Fastas', sample)
        if not os.path.isdir(bins_dir):
            print(f"Warning: bins directory not found: {bins_dir}")
            continue

        das_tool_summary = os.path.join(
            bins_base, mapper, 'run_DAS_Tool', f'{sample}_DASTool_summary.tsv')
        binner_map = load_dastool_summary(das_tool_summary)
        tax_map = load_gtdbtk(os.path.join(gtdbtk_base, mapper, sample))
        checkm2_map = load_checkm2(os.path.join(checkm2_base, mapper, sample))
        gunc_map = load_gunc(os.path.join(gunc_base, mapper, sample))

        for fa in sorted(glob.glob(os.path.join(bins_dir, '*.fa'))):
            bin_name = os.path.splitext(os.path.basename(fa))[0]
            binner = binner_map.get(bin_name, bin_name.split('.')[0])

            tax_entry = tax_map.get(bin_name, (parse_gtdbtk_classification(''), ''))
            tax_dict, full_classification = tax_entry
            tax_label = get_taxonomic_label(tax_dict)

            stats = compute_mag_stats(fa)
            qc = checkm2_map.get(bin_name, {'completeness': '', 'contamination': '', 'quality_score': ''})
            gunc = gunc_map.get(bin_name, {'gunc_clade_separation_score': '', 'gunc_pass': ''})

            all_mags.append({
                'original_name': bin_name,
                'original_path': fa,
                'sample_id': sample,
                'mapper': mapper,
                'winning_binner': binner,
                'tax_dict': tax_dict,
                'tax_label': tax_label,
                'gtdbtk_classification': full_classification,
                **qc,
                **gunc,
                **stats
            })

# Sort for reproducible global numbering
all_mags.sort(key=lambda x: (x['sample_id'], x['original_name']))

rows = []
for i, mag in enumerate(all_mags, 1):
    mag_id = f"MAG_{i:04d}"
    tax = mag['tax_dict']
    rows.append({
        'mag_id': mag_id,
        'new_name': f"{mag_id}__{mag['tax_label']}",
        'original_name': mag['original_name'],
        'original_path': mag['original_path'],
        'sample_id': mag['sample_id'],
        'mapper': mag['mapper'],
        'winning_binner': mag['winning_binner'],
        'domain': tax.get('domain', ''),
        'phylum': tax.get('phylum', ''),
        'class': tax.get('class', ''),
        'order': tax.get('order', ''),
        'family': tax.get('family', ''),
        'genus': tax.get('genus', ''),
        'species': tax.get('species', ''),
        'gtdbtk_classification': mag['gtdbtk_classification'],
        'completeness': mag.get('completeness', ''),
        'contamination': mag.get('contamination', ''),
        'quality_score': mag.get('quality_score', ''),
        'gunc_clade_separation_score': mag.get('gunc_clade_separation_score', ''),
        'gunc_pass': mag.get('gunc_pass', ''),
        'total_length_bp': mag.get('total_length_bp', ''),
        'num_contigs': mag.get('num_contigs', ''),
        'gc_percent': mag.get('gc_percent', ''),
        'N50': mag.get('N50', '')
    })

pd.DataFrame(rows).to_csv(snakemake.output.summary, sep='\t', index=False)
print(f"Wrote summary for {len(rows)} MAGs to {snakemake.output.summary}")
