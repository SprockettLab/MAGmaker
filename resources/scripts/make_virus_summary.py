import os
import sys
import pandas as pd

log_path = str(snakemake.log[0])
sys.stderr = open(log_path, 'w')
sys.stdout = sys.stderr

contig_samples = snakemake.params.contig_samples
assemblers = snakemake.params.assemblers
genomad_base = snakemake.params.genomad_base
checkv_base = snakemake.params.checkv_base


def which_assembler(sample, assemblers):
    """Return whichever configured assembler produced contigs for this sample."""
    for assembler in assemblers:
        contigs = os.path.join('output', 'assemble', assembler,
                               f'{sample}.contigs.fasta')
        if os.path.exists(contigs):
            return assembler
    return 'unknown'


def load_tsv(path):
    """Read a geNomad/CheckV TSV, returning an empty frame if missing/empty."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep='\t')
    except Exception as e:
        print(f"Warning: could not read {path}: {e}")
        return pd.DataFrame()


def load_checkv(checkv_dir):
    """Map viral contig id -> CheckV quality fields."""
    checkv_map = {}
    df = load_tsv(os.path.join(checkv_dir, 'quality_summary.tsv'))
    if df.empty:
        return checkv_map
    for _, row in df.iterrows():
        checkv_map[str(row.get('contig_id', ''))] = {
            'checkv_quality': row.get('checkv_quality', ''),
            'miuvig_quality': row.get('miuvig_quality', ''),
            'checkv_completeness': row.get('completeness', ''),
            'checkv_contamination': row.get('contamination', ''),
            'provirus': row.get('provirus', ''),
            'warnings': row.get('warnings', ''),
        }
    return checkv_map


CHECKV_BLANK = {
    'checkv_quality': '', 'miuvig_quality': '', 'checkv_completeness': '',
    'checkv_contamination': '', 'provirus': '', 'warnings': '',
}

rows = []

for assembler in assemblers:
    for sample in contig_samples:
        # Only summarize the assembler that actually produced this sample.
        if which_assembler(sample, assemblers) != assembler:
            continue

        summary_dir = os.path.join(genomad_base, assembler, sample,
                                   f'{sample}.contigs_summary')
        virus_df = load_tsv(os.path.join(
            summary_dir, f'{sample}.contigs_virus_summary.tsv'))
        plasmid_df = load_tsv(os.path.join(
            summary_dir, f'{sample}.contigs_plasmid_summary.tsv'))
        checkv_map = load_checkv(os.path.join(checkv_base, assembler, sample))

        # Viral contigs (with CheckV quality joined on seq_name).
        for _, row in virus_df.iterrows():
            seq_name = str(row.get('seq_name', ''))
            checkv = checkv_map.get(seq_name, CHECKV_BLANK)
            rows.append({
                'Sample_ID': sample,
                'Assembler': assembler,
                'Element_Type': 'virus',
                'Seq_Name': seq_name,
                'Length_BP': row.get('length', ''),
                'Topology': row.get('topology', ''),
                'N_Genes': row.get('n_genes', ''),
                'N_Hallmarks': row.get('n_hallmarks', ''),
                'Score': row.get('virus_score', ''),
                'FDR': row.get('fdr', ''),
                'Taxonomy': row.get('taxonomy', ''),
                'Conjugation_Genes': '',
                'AMR_Genes': '',
                **checkv,
            })

        # Plasmid contigs (geNomad only; CheckV does not score plasmids).
        for _, row in plasmid_df.iterrows():
            rows.append({
                'Sample_ID': sample,
                'Assembler': assembler,
                'Element_Type': 'plasmid',
                'Seq_Name': str(row.get('seq_name', '')),
                'Length_BP': row.get('length', ''),
                'Topology': row.get('topology', ''),
                'N_Genes': row.get('n_genes', ''),
                'N_Hallmarks': row.get('n_hallmarks', ''),
                'Score': row.get('plasmid_score', ''),
                'FDR': row.get('fdr', ''),
                'Taxonomy': '',
                'Conjugation_Genes': row.get('conjugation_genes', ''),
                'AMR_Genes': row.get('amr_genes', ''),
                **CHECKV_BLANK,
            })

columns = [
    'Sample_ID', 'Assembler', 'Element_Type', 'Seq_Name', 'Length_BP',
    'Topology', 'N_Genes', 'N_Hallmarks', 'Score', 'FDR', 'Taxonomy',
    'Conjugation_Genes', 'AMR_Genes', 'checkv_quality', 'miuvig_quality',
    'checkv_completeness', 'checkv_contamination', 'provirus', 'warnings',
]
pd.DataFrame(rows, columns=columns).to_csv(
    snakemake.output.summary, sep='\t', index=False)
n_virus = sum(1 for r in rows if r['Element_Type'] == 'virus')
n_plasmid = sum(1 for r in rows if r['Element_Type'] == 'plasmid')
print(f"Wrote virus summary: {n_virus} viral + {n_plasmid} plasmid contigs "
      f"to {snakemake.output.summary}")
