# MAGmaker — Developer Guide

## What this is

A Snakemake pipeline for end-to-end processing of paired-end metagenomic shotgun data. It was developed in the Moeller Lab (CUMoellerLab) and is now maintained under SprockettLab. It will eventually be superseded by MAGforge but is actively used for ongoing projects.

The pipeline has two modes driven by two separate Snakefiles:

| Snakefile | Mode | Purpose |
|---|---|---|
| `Snakefile` → `resources/snakefiles/Snakefile.smk` | Main | QC → Assembly → Prototype Selection → Taxonomic Profiling |
| `Snakefile-bin` → `resources/snakefiles/Snakefile-bin.smk` | Binning | QC → Assembly → Mapping → Binning → Bin Selection (MAGs) |

Both Snakefiles are symlinked at the repository root for convenient invocation. `config.yaml` at root is also a symlink to `resources/config/config.yaml`.

---

## Repository layout

```
MAGmaker/
├── Snakefile              # symlink → resources/snakefiles/Snakefile.smk
├── Snakefile-bin          # symlink → resources/snakefiles/Snakefile-bin.smk
├── config.yaml            # symlink → resources/config/config.yaml
└── resources/
    ├── config/
    │   ├── config.yaml        # all pipeline parameters and tool settings
    │   ├── samples.txt        # sample metadata (tab-separated, Sample column required)
    │   ├── units.txt          # read file paths (Sample, Unit, R1, R2)
    │   ├── binning.txt        # binning groups (Sample, Contigs, Read_Groups, Contig_Groups)
    │   └── reference_list.txt
    ├── env/               # conda environment YAML files (one per tool/module)
    ├── snakefiles/        # modular rule files included by the Snakefiles
    ├── scripts/           # helper Python, R, and Perl scripts
    ├── notebooks/         # Jupyter notebooks for distributed job setup
    ├── test/
    │   ├── test_reads/    # tiny paired FASTQ files (amy, bob, carl, diane)
    │   └── test_dbs/      # pre-built bowtie2 index (GCA_000001635.9 mouse)
    └── db/
        ├── bt2/           # host genome FASTA (downloaded and indexed at runtime)
        └── metaquast/     # reference genome for MetaQUAST
```

---

## Pipeline steps

### Main pipeline (`Snakefile`)

Includes: `qc.smk`, `assemble.smk`, `prototype_selection.smk`, `profile.smk`

```
raw reads (units.txt)
  └─ fastqc_pre_trim
  └─ cutadapt_pe          (adapter trimming, quality trimming)
  └─ fastqc_post_trim
  └─ merge_units           (cat reads across sequencing units per sample)
  └─ host_bowtie2_build    (build index if needed)
  └─ host_filter           (bowtie2; nonhost reads retained)
  └─ fastqc_post_host
  └─ multiqc               → output/qc/multiqc/multiqc.html

nonhost reads
  └─ metaspades / megahit  (configurable; one or both)
  └─ quast / metaquast
  └─ multiqc_assemble      → output/assemble/multiqc_assemble/multiqc.html

nonhost reads (sourmash)
  └─ sourmash_sketch_reads
  └─ sourmash_dm
  └─ sourmash_plot         → output/prototype_selection/sourmash_plot/
  └─ prototype_selection   → output/prototype_selection/prototype_selection/selected_prototypes.yaml

nonhost reads (profiling)
  └─ metaphlan             (MetaPhlAn4 preferred)
  └─ merge_metaphlan_tables → output/profile/metaphlan/merged_abundance_table.txt
```

**Default `rule all` targets:** multiqc HTML, assembly multiqc HTML, sourmash plot, selected_prototypes.yaml, merged MetaPhlAn table.

### Binning pipeline (`Snakefile-bin`)

Includes: `qc.smk`, `assemble.smk`, `mapping.smk`, `binning.smk`, `selected_bins.smk`, `mag_qc.smk`

The `binning.txt` file defines which samples' reads get mapped to which samples' contigs. `Read_Groups` and `Contig_Groups` columns contain comma-separated group labels; all read-samples in a group are mapped to all contig-samples in the same group.

```
contigs (from assembly) + nonhost reads
  └─ index_contigs (bowtie2 or minimap2)
  └─ map_reads_[bt2|minimap2]
  └─ sort_index_bam

sorted BAMs + contigs
  ├─ metabat2 path:  make_metabat2_coverage_table → run_metabat2
  ├─ maxbin2 path:   make_maxbin2_coverage_table → make_maxbin2_abund_list → run_maxbin2
  └─ concoct path:   cut_up_fasta → make_concoct_coverage_table → run_concoct → merge_cutup_clustering → extract_fasta_bins

all bins → [Fasta_to_Scaffolds2Bin for each binner] → run_DAS_Tool → consolidate_DAS_Tool_bins
  → output/selected_bins/{mapper}/DAS_Tool_Fastas/{sample}/   (per-sample subdirectory)

MAGs (per sample)
  ├─ run_checkm2    → output/mag_qc/checkm2/{mapper}/{sample}/quality_report.tsv
  ├─ run_gunc       → output/mag_qc/gunc/{mapper}/{sample}/   (.done touch file)
  └─ run_gtdbtk     → output/mag_qc/gtdbtk/{mapper}/{sample}/ (.done touch file)

all samples complete
  └─ make_mag_summary  → output/mag_qc/mag_summary.tsv   (globally numbered, editable)
  └─ rename_mags       → output/mag_qc/renamed_mags/*.fa  (reads new_name from table)
```

**Top-level rules in Snakefile-bin:** `select_bins`, `bin_all`, `map_all`, `make_mag_summary`, `rename_mags`

**MAG renaming workflow:** `make_mag_summary` assigns global IDs (`MAG_0001`…`MAG_N`) sorted by sample then bin name, looks up GTDB-tk taxonomy to build the `new_name` label (`MAG_0001__Roseburia_intestinalis`), and writes the full table. Users can edit `new_name` in `mag_summary.tsv` and re-run `snakemake rename_mags` — Snakemake detects the edited file is newer and re-runs automatically. The `rename_mags` rule clears the `renamed_mags/` directory before copying so stale files don't accumulate.

---

## Configuration

### `config.yaml` — key sections

| Section | Purpose | Notes |
|---|---|---|
| `samples` / `units` / `binning` | paths to metadata TSVs | |
| `assemblers` | list: `metaspades`, `megahit` | comment out to skip one |
| `mappers` | list: `minimap2`, `bowtie2` | used in binning pipeline |
| `binners` | list: `concoct`, `metabat2`, `maxbin2` | |
| `params.cutadapt` | adapter sequences and trim settings | |
| `params.metaphlan` | `db_path` and `db_name` — **must be set per cluster** | MetaPhlAn4 uses different DB names |
| `params.gtdbtk.db_path` | **required for MAG taxonomy** — path to GTDB-tk reference data (~85 GB) | download manually from data.gtdb.ecogenomics.org; see README |
| `params.checkm2.db_path` | optional — CheckM2 diamond database; auto-downloaded if empty | |
| `params.gunc.db_path` | optional — GUNC database file; auto-downloaded if empty | |
| `params.kraken2` / `params.bracken` | DB paths — **must be set per cluster** | |
| `host_filter.genome` | path to host genome FASTA | varies by project |
| `host_filter.db_dir` | directory for bowtie2 index of host genome | |
| `threads.*` | per-rule thread counts | scale to cluster node sizes; checkm2=16, gunc=8, gtdbtk=16 |
| `mem_mb.*` | memory ceilings in MB | assembly rules auto-scale based on input size; checkm2=32000, gtdbtk=128000 |

### `samples.txt`

Tab-separated. First column must be named `Sample`. Additional columns are carried through but not used by the pipeline.

### `units.txt`

Tab-separated. Columns: `Sample`, `Unit`, `R1`, `R2`. One row per sequencing unit (e.g., lane or run). Multiple units per sample are concatenated by `merge_units`. Paths in R1/R2 can be absolute or relative to the working directory.

### `binning.txt`

Tab-separated. Columns: `Sample`, `Contigs`, `Read_Groups`, `Contig_Groups`. Controls which reads get mapped to which assemblies for binning. `Read_Groups` and `Contig_Groups` are comma-separated group labels. Samples with the same group label in both columns are paired for mapping. Samples with an empty `Read_Groups` or `Contig_Groups` are assembly-only or read-only (not mapped).

---

## Conda environments

Each module has its own environment YAML in `resources/env/`:

| File | Used by | Notes |
|---|---|---|
| `snakemake.yaml` | top-level | install this first |
| `qc.yaml` | qc.smk | bowtie2, samtools, pigz, tbb=2020.2 (pinned) |
| `assemble.yaml` | assemble.smk | spades, megahit, quast |
| `mapping.yaml` | mapping.smk | bowtie2, samtools, minimap2, tbb=2020.2 (pinned) |
| `binning.yaml` | binning.smk | metabat2, maxbin2, fraggenescan |
| `concoct_linux.yaml` | binning.smk (concoct rules) | Linux-specific |
| `concoct_osx.yaml` | — | macOS-specific (not used on clusters) |
| `selected_bins.yaml` | selected_bins.smk | das_tool |
| `mag_qc.yaml` | mag_qc.smk | checkm2, gunc, pandas, numpy |
| `gtdbtk.yaml` | mag_qc.smk (run_gtdbtk) | gtdbtk>=2.0; separate env due to dependency conflicts |
| `profile.yaml` | profile.smk | kraken2, bracken, krona, metaphlan, humann |
| `prototype_selection.yaml` | prototype_selection.smk | sourmash, scikit-bio |
| `sourmash.yaml` | sourmash.smk (standalone) | sourmash, scikit-bio |

**Important:** `tbb=2020.2` is pinned in `qc.yaml` and `mapping.yaml` because newer TBB versions break bowtie2 compatibility on some systems. Do not remove this pin without testing.

The `snakemake` environment uses unpinned `snakemake>=6.4.0` and `scikit-bio>=0.5`. The `prototype_selection.yaml` pins `scikit-bio=0.5` for stability.

---

## Installation

```bash
# Clone the repository
git clone https://github.com/SprockettLab/MAGmaker.git
cd MAGmaker

# Install mamba (if not already available)
conda install -n base -c conda-forge mamba

# Create and activate the snakemake environment
mamba env create -n snakemake -f resources/env/snakemake.yaml
conda activate snakemake
```

All other environments (qc, assemble, mapping, etc.) are created automatically by Snakemake via `--use-conda` on first run.

---

## Running the pipeline

### Local / interactive

```bash
conda activate snakemake

# Dry run (check DAG without executing)
snakemake --cores 8 --use-conda -n

# Main pipeline
snakemake --cores 8 --use-conda

# Binning pipeline
snakemake --snakefile Snakefile-bin --cores 8 --use-conda

# Specific target rules
snakemake --cores 8 --use-conda select_bins
snakemake --cores 8 --use-conda bin_all
snakemake --cores 8 --use-conda map_all
```

### SLURM (demon cluster)

A Snakemake profile for demon is in `resources/profiles/demon/config.yaml`. It uses the `defq` partition (the default, open to all users) with 8 GB RAM and 120-minute walltime as defaults.

```bash
conda activate snakemake
snakemake --profile resources/profiles/demon
```

For the binning pipeline:
```bash
snakemake --snakefile Snakefile-bin --profile resources/profiles/demon
```

Override resources for specific rules on the command line if needed:
```bash
snakemake --profile resources/profiles/demon \
  --set-resources megahit:mem_mb=256000 megahit:runtime=2880
```

**Demon `defq` partition limits:**
- Max 15 nodes per job
- Max walltime: 90 days
- Default memory: 4 GB/CPU (we override this explicitly via `--mem`)
- SLURM logs go to `output/slurm_logs/{rule}/`

**Note:** Assembly rules (megahit, metaspades) request 256 GB RAM via `mem_mb` in `config.yaml`. demon nodes in defq have sufficient RAM — confirm with `sinfo -N -o "%N %m"` if jobs fail with OOM errors.

---

## Test data

The repository includes a small test dataset:
- Reads: `resources/test/test_reads/{amy,bob,carl,diane}_{R1,R2}.fastq.gz`
- Pre-built bowtie2 index: `resources/test/test_dbs/GCA_000001635.9.*` (mouse genome, for testing host filter)
- The default `samples.txt`, `units.txt`, and `binning.txt` are configured to use this test data

To run a test, set `host_filter.genome` in `config.yaml` to the test genome FASTA and run with `--cores 4`.

---

## Known issues and quirks

- **MEGAHIT memory units:** The shell command passes `$(({resources.mem_mb}*1024*1024))` bytes. Since `mem_mb` is in megabytes, this is `MB × 1024 × 1024 = bytes`, which is correct for MEGAHIT's `--memory` flag (bytes). Do not change this. Assembly rules now auto-scale `mem_mb` based on input size (`max(16000, input.size_mb * 10)`) up to the config ceiling.
- **`sourmash.smk` is standalone:** This file is not included by either Snakefile. It is an older standalone version with different output paths. Use `prototype_selection.smk` instead.
- **`concoct` uses `concoct_linux.yaml`:** Hard-coded to the Linux environment file. Will not work on macOS without changing the `conda:` directive to `concoct_osx.yaml`.
- **Mixed Snakemake wrapper versions:** `qc.smk` uses `0.72.0/bio/fastqc` and `0.17.4/bio/cutadapt/pe` (old format) but `v3.1.0/bio/multiqc` (new format). Consider updating all wrappers to a consistent version.
- **GTDB-tk archaeal summary may be absent:** `run_gtdbtk` uses a `.done` touch file as output rather than the actual summary TSVs, because `gtdbtk.ar53.summary.tsv` is not produced when no archaeal MAGs are present. The `make_mag_summary` script handles this by checking for both `bac120` and `ar53` files with `os.path.exists()`.
- **GUNC output filename varies by DB version:** The script globs for `*.maxCSS_level.tsv` rather than hardcoding the filename, since the exact name includes the database version string.
- **TODO — Build Bracken database for Kraken2_db_Standard:** Bracken requires `.kmer_distrib` files generated separately from the Kraken2 build. Run on demon before using the Kraken2/Bracken profiling rules:
  ```bash
  bracken-build -d /isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/kraken2/Kraken2_db_Standard \
    -t <threads> -k 35 -l 150
  ```
  The `-l` flag should match your read length (150 bp is typical). Run once per read length. This is a prerequisite for `taxonomy_kraken` to succeed.
- **TODO — Download GTDB-tk reference data for demon:** `params.gtdbtk.db_path` is set in `config.yaml` but the data has not yet been downloaded. Note: `gtdbtk download_db` is not a valid subcommand — manual download is required. Browse `https://data.gtdb.ecogenomics.org/releases/` for the latest release, download the `gtdbtk_data_r*.tar.gz` with `wget -c`, extract into `/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/gtdbtk/`, then update `params.gtdbtk.db_path` to the extracted directory.

---

## GitHub migration (CUMoellerLab → SprockettLab)

To copy the repository to the SprockettLab GitHub organization while preserving full git history:

```bash
# 1. Create a new empty repository at github.com/SprockettLab/MAGmaker
#    (do not initialize with README/license)

# 2. Add the new remote
git remote add sprockettlab https://github.com/SprockettLab/MAGmaker.git

# 3. Push all branches and tags
git push sprockettlab --all
git push sprockettlab --tags

# 4. Make SprockettLab the default origin
git remote rename origin cumoellerlab
git remote rename sprockettlab origin
```

---

## Planned improvements

- Standardize Snakemake wrapper versions across all rule files (qc.smk mixes 0.72.0/bio/fastqc, 0.17.4/bio/cutadapt/pe, and v3.1.0/bio/multiqc)
- Evaluate whether `sourmash.smk` should be deleted (it is a standalone orphan not included by either Snakefile; `prototype_selection.smk` is the active version)
- Add co-assembly mode: pool reads across samples → single MEGAHIT run → map all samples back
- Add HUMAnN 3 functional profiling downstream of MetaPhlAn
- Add StrainPhlAn strain-level profiling (uses the .sam.bz2 files already produced by the metaphlan rule)
- Add MultiQC step for the binning pipeline mapping outputs
