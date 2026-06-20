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

After `make_mag_summary` completes, `output/mag_qc/mag_summary.tsv` contains one row per MAG. MAGs are sorted by GTDB-tk taxonomy (domain → species) before numbering, so sequential IDs (`MAG_0001`, `MAG_0002`, …) group related organisms together regardless of which sample they came from. Empty taxonomy fields are written as `NA`.

| Column | Description |
|---|---|
| `MAG_ID` | Global sequential ID (`MAG_0001` … `MAG_N`), sorted by taxonomy |
| `New_Name` | Proposed FASTA filename — **edit this column to rename MAGs** |
| `Original_Name` | DAS_Tool bin name |
| `Original_Path` | Path to source FASTA |
| `Sample_ID` | Sample the MAG was assembled from |
| `Assembler` | Assembler that produced the contigs (megahit / metaspades) |
| `Winning_Binner` | Which binner DAS_Tool selected (metabat2 / maxbin2 / concoct) |
| `Domain` … `Species` | GTDB-tk taxonomy in separate columns; `NA` if unclassified |
| `GTDB_Classification` | Full GTDB-tk classification string |
| `Completeness` | CheckM2 completeness (%) |
| `Contamination` | CheckM2 contamination (%) |
| `Quality_Score` | Completeness − 5 × Contamination |
| `MIMAG_Quality` | MIMAG quality tier: `HQ` (≥90% complete, <5% contamination), `MQ` (≥50%, <10%), `LQ` (all else) |
| `GUNC_Clade_Separation_Score` | GUNC chimera score |
| `GUNC_Pass` | Whether the bin passes GUNC QC |
| `Total_Length_BP` | Total genome size (bp) |
| `Num_Contigs` | Number of contigs in the bin |
| `Largest_Contig` | Length of the longest contig (bp) |
| `GC_Percent` | GC content (%) |
| `N50` | Assembly N50 (bp) |
| `Coding_Density` | Fraction of genome that is coding sequence (from CheckM2) |
| `Total_Coding_Sequences` | Number of predicted coding sequences (from CheckM2) |

---

## MAG renaming workflow

The default `New_Name` values follow the pattern `MAG_0001__Genus_species` using the most resolved available GTDB-tk taxonomy rank. To use custom names:

1. Open `output/mag_qc/mag_summary.tsv` in a spreadsheet editor or text editor
2. Edit the `New_Name` column as desired
3. Save the file
4. Re-run the rename step — Snakemake detects the table is newer than the output and re-runs automatically:

```bash
# Local
snakemake --snakefile Snakefile-bin --cores 4 --use-conda rename_mags \
  --config binning=output/config/auto_binning.txt

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
