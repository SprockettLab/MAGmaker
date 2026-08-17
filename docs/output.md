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
│   ├── {mapper}/DAS_Tool_Fastas/{sample}/  DAS_Tool-selected MAG bins per sample
│   └── {mapper}/Binette_Fastas/{sample}/   Binette-selected MAG bins per sample
│
└── mag_qc/                             (binning pipeline)
    ├── checkm2/{mapper}/{sample}/      CheckM2 quality reports
    ├── gunc/{mapper}/{sample}/         GUNC chimera detection results
    ├── gtdbtk/{mapper}/{sample}/       GTDB-tk taxonomy outputs
    ├── mag_summary.tsv                 combined MAG table — editable
    └── renamed_mags/                   final MAG FASTAs with user-defined names
```

Only one of `DAS_Tool_Fastas/` and `Binette_Fastas/` is produced per run,
chosen by `params.consolidation.tool` or by `--binette` / `--das-tool`.
Everything under `mag_qc/` reads whichever one the run selected, so
`mag_summary.tsv` describes that tool's MAGs. Running an arm both ways
leaves both directories in place, which is how the two are compared.

For Binette, `Winning_Binner` in `mag_summary.tsv` can name more than one
binner (`metabat2+concoct`). Binette builds candidate bins from the
intersection, difference and union of overlapping input bins, so a selected
bin is not always one that a single binner produced.

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
| `MSA_Percent` | Share of the concatenated marker alignment this genome fills with residues rather than gaps, from GTDB-tk. What determines how well a genome can be placed in a phylogeny, and not the same as completeness: a genome can recover most markers as fragments and still fill little of the alignment. `NA` if GTDB-tk did not place it |
| `GTDBTk_Warnings` | GTDB-tk's own caveats about the placement, `NA` if none |
| `Strain_Heterogeneity` | Percentage of evaluated positions carrying more than one allele, from CMSeq, when the sample's own reads are mapped back to its own MAGs. Detects a MAG that is a clean consensus of several co-resident strains of one species, which CheckM2 and GUNC cannot see. `NA` when CMSeq did not run or the MAG had too few covered positions |
| `SH_Positions_Evaluated` | Positions covered at least 10x with base quality above 30. Fewer than 100 reports `NA` for the rate, following Pasolli 2019 and Sanders 2023 |
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
| `Notes` | `NA` for a MAG; on a `NONE` row, why that sample produced none |

### Samples that produced no MAGs

A sample can finish the pipeline without yielding a MAG, either because
every binner declined it or because nothing survived consolidation. These
appear after the numbered MAGs as rows with `MAG_ID` = `NONE`, carrying the
sample in `Sample_ID` and the reason in `Notes`:

```
NONE  NONE  ...  C59_R_TP4  ...  no binner produced bins; declined by: concoct, maxbin2, metabat2, semibin2
```

Without these rows such a sample is simply absent from the table, which is
indistinguishable from one that was never run. `rename_mags` skips them, so
they do not affect `renamed_mags/`. Filter them out with
`MAG_ID != "NONE"` before any per-MAG analysis.

The most common cause is an assembly too fragmented to bin. SemiBin2 needs a
contig of at least 4000 bp to form must-link pairs, and MetaBAT2, MaxBin2
and CONCOCT have their own minimum-length and marker-gene requirements.

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

## GTDB-Tk stages

`run_gtdbtk` runs `classify_wf` by default. Two flags shorten it:

| flag | stages run | taxonomy columns | `MSA_Percent` | memory reserved |
|---|---|---|---|---|
| *(none)* | identify, align, classify | populated | populated | 128-320 GB |
| `--skip-gtdbtk-classify` | identify, align | NA | populated | 32-64 GB |
| `--skip-gtdbtk` | none | NA | NA | rule does not run |

Nearly all of GTDB-Tk's cost is `classify`, which runs pplacer. GTDB-Tk's
documentation puts the bacterial requirement at about 140 GB and attributes
it there. `identify` and `align` do not run pplacer, so skipping classify
changes which nodes the rule can be scheduled on rather than merely
trimming its runtime.

`MSA_Percent` survives that, which is the point of the flag. GTDB-Tk defines
it as "the percentage of the MSA spanned by the genome (i.e. percentage of
columns with an amino acid)", and `align` writes that alignment as
`gtdbtk.<domain>.user_msa.fasta.gz`. `make_mag_summary` counts non-gap
columns there, so the value is GTDB-Tk's own number rather than an
approximation of it. An analysis that filters genomes on how much of the
marker alignment they fill, but takes taxonomy from elsewhere, never needs
`classify` at all.

The align-only path passes `--min_perc_aa 0`. GTDB-Tk's default of 10 drops
genomes below that threshold into `filtered.tsv` instead of reporting a low
value for them, and a genome absent from the output is not the same as one
measured at 4%.

### Adding classification later

The three stages are three rules, each with its own sentinel
(`.identify.done`, `.align.done`, `.done`). Re-running without
`--skip-gtdbtk-classify` therefore runs **classify only**: identify and
align are already done and Snakemake sees them as up to date. Nothing is
recomputed and nothing needs deleting.

Runs made before the split have `classify_wf`'s output but none of the
per-stage sentinels. Both upstream rules detect complete output and reuse
it rather than recomputing, and stamp the sentinel with that output's own
mtime, so a finished arm is not re-classified just because the sentinels
are new.

One consequence worth knowing: with `min_perc_aa_no_classify: 0`, an
align-only run followed later by classify will place genomes that a plain
`classify_wf` run would have filtered at the default floor of 10. GTDB-Tk
flags those in its warnings, which `mag_summary.tsv` carries as
`GTDBTk_Warnings`.
