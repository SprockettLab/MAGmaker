# Output

---

## Directory layout

```
output/
├── qc/
│   ├── fastqc/                           per-sample FastQC reports (pre-trim, post-trim, post-host)
│   ├── fastp/                            fastp JSON + HTML reports (if trimmer: fastp)
│   ├── cutadapt/                         Cutadapt logs (if trimmer: cutadapt)
│   ├── host_filter/nonhost/              host-filtered FASTQ files (input to assembly)
│   └── multiqc/multiqc.html             combined QC report
│
├── assemble/
│   ├── megahit/{sample}.contigs.fasta   MEGAHIT assemblies
│   ├── metaspades/{sample}/             metaSPAdes assemblies
│   └── multiqc_assemble/multiqc.html   assembly QC report (QUAST stats)
│
├── prototype_selection/
│   ├── sourmash_plot/                   pairwise MinHash similarity heatmap
│   └── prototype_selection/
│       └── selected_prototypes.yaml    representative sample IDs (input to binning config)
│
├── config/
│   └── auto_binning.txt               auto-generated binning config (from generate_binning_config)
│
├── profile/
│   └── metaphlan/
│       ├── {sample}.txt               per-sample MetaPhlAn 4 profiles
│       └── merged_abundance_table.txt merged taxonomy table (all samples)
│
├── selected_bins/                      (binning pipeline)
│   └── {mapper}/DAS_Tool_Fastas/{sample}/  DAS_Tool-selected MAG bins per sample
│
└── mag_qc/                             (binning pipeline)
    ├── checkm2/{mapper}/{sample}/      CheckM2 quality reports
    ├── gunc/{mapper}/{sample}/         GUNC chimera detection results
    ├── gtdbtk/{mapper}/{sample}/       GTDB-tk taxonomy outputs
    ├── mag_summary.tsv                 combined MAG table — editable
    └── renamed_mags/                   final MAG FASTAs with user-defined names
```

---

## MAG summary table

After `make_mag_summary` completes, `output/mag_qc/mag_summary.tsv` contains one row per MAG with the following columns:

| Column | Description |
|---|---|
| `mag_id` | Global sequential ID (`MAG_0001` … `MAG_N`), sorted by sample then bin name |
| `new_name` | Proposed FASTA filename — **edit this column to rename MAGs** |
| `original_name` | DAS_Tool bin name |
| `original_path` | Path to source FASTA |
| `sample_id` | Sample the MAG was assembled from |
| `winning_binner` | Which binner DAS_Tool selected (metabat2 / maxbin2 / concoct) |
| `domain` … `species` | GTDB-tk taxonomy in separate columns |
| `gtdbtk_classification` | Full GTDB-tk classification string |
| `completeness` | CheckM2 completeness (%) |
| `contamination` | CheckM2 contamination (%) |
| `quality_score` | CheckM2 quality score |
| `gunc_clade_separation_score` | GUNC chimera score |
| `gunc_pass` | Whether the bin passes GUNC QC |
| `total_length_bp` | Total assembly size |
| `num_contigs` | Number of contigs in the bin |
| `gc_percent` | GC content |
| `N50` | Assembly N50 |

---

## MAG renaming workflow

The default `new_name` values follow the pattern `MAG_0001__Genus_species` using GTDB-tk taxonomy. To use custom names:

1. Open `output/mag_qc/mag_summary.tsv` in a spreadsheet editor or text editor
2. Edit the `new_name` column as desired
3. Save the file
4. Re-run the rename step — Snakemake detects the table is newer than the output and re-runs automatically:

```bash
# Local
snakemake --cores 4 --use-conda rename_mags \
  --snakefile Snakefile-bin --config binning=output/config/auto_binning.txt

# demon
snakemake --snakefile Snakefile-bin --profile resources/profiles/demon rename_mags \
  --config binning=output/config/auto_binning.txt
```

The `rename_mags` rule clears `renamed_mags/` before copying so stale files from previous runs don't accumulate.

---

## Logs

- **Rule-level logs** (tool stdout/stderr): `output/logs/{rule}/{sample}.log`
- **SLURM job logs** (cluster submission details): `.snakemake/slurm_logs/{rule}/`

---

← [Back to README](../README.md)
