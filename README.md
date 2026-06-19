# MAGmaker

A Snakemake pipeline for end-to-end processing of paired-end metagenomic shotgun sequencing data. MAGmaker takes raw FASTQ files through quality control, host read removal, assembly, taxonomic profiling, and — optionally — binning to produce metagenome-assembled genomes (MAGs).

Developed by the [Moeller Lab](https://moellerlab.com) at Cornell University and Princeton University, and maintained by the [Sprockett Lab](https://www.sprockettlab.com/) at Wake Forest University School of Medicine.

---

## Pipeline overview

MAGmaker has two modes, each driven by a separate Snakefile.

### Main pipeline (`Snakefile`)

Runs QC through taxonomic profiling and prototype (representative sample) selection.

```
raw reads (FASTQ)
  ├── FastQC (pre-trim)
  ├── fastp / Cutadapt      adapter and quality trimming (configurable; fastp default)
  ├── FastQC (post-trim)
  ├── merge across lanes    symlink if single unit; cat if multiple
  ├── bowtie2               host read removal
  ├── FastQC (post-host)
  └── MultiQC               → output/qc/multiqc/multiqc.html

non-host reads
  ├── MEGAHIT / metaSPAdes  de novo assembly (one or both, configurable)
  ├── QUAST                 assembly quality assessment
  └── MultiQC               → output/assemble/multiqc_assemble/multiqc.html

non-host reads (similarity)
  ├── sourmash sketch       MinHash sketches per sample
  ├── sourmash compare      pairwise distance matrix
  ├── sourmash plot         → output/prototype_selection/sourmash_plot/
  └── prototype selection   → output/prototype_selection/prototype_selection/selected_prototypes.yaml

non-host reads (taxonomy)
  ├── MetaPhlAn 4           marker-gene taxonomic profiling
  └── merge tables          → output/profile/metaphlan/merged_abundance_table.txt
```

### Binning pipeline (`Snakefile-bin`)

Extends the main pipeline with read mapping, binning, bin refinement, quality control, and taxonomy.

```
non-host reads + assemblies
  ├── bowtie2 / minimap2    map reads to contigs (configurable)
  └── samtools              sort and index BAMs

sorted BAMs + contigs
  ├── MetaBAT2              coverage-based binning
  ├── MaxBin2               coverage-based binning
  └── CONCOCT               composition + coverage binning

all bins
  └── DAS_Tool              bin refinement and dereplication → MAGs

MAGs
  ├── CheckM2               completeness and contamination assessment
  ├── GUNC                  chimeric bin detection
  ├── GTDB-tk               taxonomic classification
  └── make_mag_summary      → output/mag_qc/mag_summary.tsv

mag_summary.tsv  (editable)
  └── rename_mags           → output/mag_qc/renamed_mags/*.fa
```

---

## Requirements

- Linux or macOS
- [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/en/latest/)
- ~20 GB disk space for conda environments (built automatically on first run)
- A host genome FASTA for filtering (e.g., human GRCh38, mouse GRCm39)
- A [MetaPhlAn 4 database](https://huttenhower.sph.harvard.edu/metaphlan/) if running taxonomic profiling
- A [GTDB-tk reference package](https://ecogenomics.github.io/GTDBTk/) (~85 GB) if running MAG taxonomy — see [Database setup](#database-setup) below

All other tool dependencies (FastQC, Cutadapt, MEGAHIT, MetaPhlAn, sourmash, CheckM2, GUNC, etc.) are installed automatically by Snakemake into isolated conda environments.

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/SprockettLab/MAGmaker.git
cd MAGmaker

# 2. Install mamba if you don't already have it
conda install -n base -c conda-forge mamba

# 3. Create and activate the Snakemake environment
mamba env create -n snakemake -f resources/env/snakemake.yaml
conda activate snakemake
```

All other environments (qc, assemble, profile, etc.) are created automatically the first time a rule runs, via `--use-conda`.

---

## Database setup

### CheckM2 and GUNC

CheckM2 and GUNC will automatically download their databases on first use if no path is set in `config.yaml`. Downloads land in the tool's default cache directory (`~/.cache/checkm2/` and `~/.gunc/` respectively) and will repeat for every new user unless a shared path is configured.

To pre-download to a shared location:

```bash
mamba create -n db_setup -c conda-forge -c bioconda checkm2 gunc -y
conda activate db_setup

checkm2 database --download --path /your/shared/dbs/checkm2/
gunc download_db /your/shared/dbs/gunc/

conda deactivate

# Note the exact filenames that were created
ls /your/shared/dbs/checkm2/
ls /your/shared/dbs/gunc/
```

Then set `params.checkm2.db_path` to the downloaded `.dmnd` file and `params.gunc.db_path` to the downloaded `.db` file in `config.yaml`. Leaving either path empty restores auto-download behavior.

### GTDB-tk

GTDB-tk requires a reference data package (~85 GB) that must be downloaded before the taxonomy step will run. Set `params.gtdbtk.db_path` in `config.yaml` to the directory containing the extracted data.

Download the data package directly from the [GTDB data server](https://data.gtdb.ecogenomics.org/releases/). Browse to the latest release directory, find the `gtdbtk_data_r*.tar.gz` package, and download it with `wget`:

```bash
mkdir -p /your/dbs/gtdbtk
cd /your/dbs/gtdbtk

# Replace the URL with the current package from data.gtdb.ecogenomics.org/releases/
wget -c https://data.gtdb.ecogenomics.org/releases/release232/<path/to/gtdbtk_data_r232.tar.gz>
tar -xzf gtdbtk_data_r232.tar.gz

# Check what directory was created, then set that path in config.yaml
ls /your/dbs/gtdbtk/
```

The `-c` flag on `wget` allows resuming interrupted downloads. After extraction, set `params.gtdbtk.db_path` to the directory containing the unpacked reference data.

See the [GTDB-tk documentation](https://ecogenomics.github.io/GTDBTk/installing/index.html) and the [GTDB releases page](https://gtdb.ecogenomics.org/) for the current database version.

> **Sprockett Lab (demon cluster):** All three databases are pre-downloaded and configured in `config.yaml` by default:
> - CheckM2: `.../dbs/checkm2/CheckM2_database/uniref100.KO.1.dmnd`
> - GUNC: `.../dbs/gunc/gunc_db_progenomes2.1.dmnd`
> - GTDB-tk: `.../dbs/gtdbtk/release232/`

---

## Configuration

All input files and parameters live in `resources/config/`. The `config.yaml` at the repository root is a symlink to `resources/config/config.yaml`.

### `samples.txt`

Tab-separated. One row per sample. The first column must be named `Sample`. Additional columns are carried through but not used by the pipeline.

```
Sample    Subject    Timepoint
amy       A          T1
bob       B          T1
```

### `units.txt`

Tab-separated. One row per sequencing unit (e.g., lane or run). Columns: `Sample`, `Unit`, `R1`, `R2`. Multiple units per sample are automatically concatenated before assembly.

```
Sample    Unit      R1                        R2
amy       Run_1     /path/to/amy_R1.fastq.gz  /path/to/amy_R2.fastq.gz
bob       Run_1     /path/to/bob_R1.fastq.gz  /path/to/bob_R2.fastq.gz
```

### `config.yaml`

Key settings to review before running:

| Parameter | Description |
|---|---|
| `trimmer` | Trimmer to use: `fastp` (default) or `cutadapt` |
| `assemblers` | Which assembler(s) to use: `metaspades`, `megahit`, or both |
| `host_filter.genome` | Path to host genome FASTA (used to build bowtie2 index) |
| `host_filter.db_dir` | Directory for the bowtie2 host index |
| `params.metaphlan.db_path` | Path to MetaPhlAn 4 database directory |
| `params.metaphlan.db_name` | MetaPhlAn 4 database name (e.g., `mpa_vJan25_CHOCOPhlAnSGB_202503`) |
| `params.gtdbtk.db_path` | **Required for MAG taxonomy** — path to GTDB-tk reference data directory |
| `params.checkm2.db_path` | Optional — path to CheckM2 diamond database; auto-downloaded if empty |
| `params.gunc.db_path` | Optional — path to GUNC database file; auto-downloaded if empty |
| `threads.*` | Per-rule thread counts — scale to your hardware |
| `mem_mb.*` | Memory limits — assembly rules auto-scale based on input size (up to the configured ceiling) |

### `binning.txt` (binning pipeline only)

Tab-separated. Defines which reads get mapped to which assemblies. Columns: `Sample`, `Contigs`, `Read_Groups`, `Contig_Groups`. Samples sharing a group label in both `Read_Groups` and `Contig_Groups` are paired for mapping. See `resources/config/binning.txt` for an example.

---

## Running the pipeline

### Local / interactive

```bash
conda activate snakemake

# Dry run — check DAG without executing
snakemake --cores 8 --use-conda -n

# Main pipeline
snakemake --cores 8 --use-conda

# Binning pipeline (through DAS_Tool)
snakemake --snakefile Snakefile-bin --cores 8 --use-conda select_bins

# MAG QC and taxonomy (runs after select_bins)
snakemake --snakefile Snakefile-bin --cores 8 --use-conda make_mag_summary

# Rename MAGs using the summary table
# Edit output/mag_qc/mag_summary.tsv first if you want custom names, then:
snakemake --snakefile Snakefile-bin --cores 8 --use-conda rename_mags
```

### SLURM (demon cluster)

A SLURM profile for the demon cluster is included at `resources/profiles/demon/`. It submits jobs to the `defq` partition with 8 GB RAM and 120-minute walltime as defaults, and automatically manages conda environments in a shared location.

```bash
conda activate snakemake

# Main pipeline
snakemake --profile resources/profiles/demon

# Binning pipeline (through DAS_Tool)
snakemake --snakefile Snakefile-bin --profile resources/profiles/demon select_bins

# MAG QC, taxonomy, and renaming
snakemake --snakefile Snakefile-bin --profile resources/profiles/demon make_mag_summary
snakemake --snakefile Snakefile-bin --profile resources/profiles/demon rename_mags
```

To override resources for a specific rule:

```bash
snakemake --profile resources/profiles/demon \
  --set-resources megahit:mem_mb=256000 megahit:runtime=2880
```

SLURM logs are written to `.snakemake/slurm_logs/{rule}/`. Rule logs (stderr/stdout from the tools themselves) go to `output/logs/{rule}/`.

> **Note:** Assembly rules auto-scale memory based on input size (up to the `mem_mb` ceiling in `config.yaml`). GTDB-tk requests 128 GB by default; MetaPhlAn 4 requests 32 GB. Both can be overridden with `--set-resources`.

### Running on test data

The repository includes a small test dataset (four samples: amy, bob, carl, diane) with pre-built bowtie2 index files. The default `samples.txt`, `units.txt`, and `config.yaml` point to this test data. To run a quick test:

```bash
snakemake --cores 4 --use-conda -n   # dry run first
snakemake --cores 4 --use-conda
```

---

## Output

```
output/
├── qc/
│   ├── multiqc/multiqc.html              QC report (pre-trim through host removal)
│   └── host_filter/nonhost/             host-filtered FASTQ files
├── assemble/
│   ├── megahit/{sample}/                MEGAHIT assemblies
│   ├── metaspades/{sample}/             metaSPAdes assemblies
│   ├── assembly_stats.tsv               merged per-sample assembly stats (all assemblers)
│   └── multiqc_assemble/multiqc.html    assembly QC report
├── prototype_selection/
│   ├── sourmash_plot/                   pairwise similarity heatmap
│   └── prototype_selection/
│       └── selected_prototypes.yaml     representative sample selection
├── profile/
│   └── metaphlan/
│       ├── profiles/{sample}.txt        per-sample MetaPhlAn profiles
│       └── merged_abundance_table.txt   merged taxonomy table (all samples)
├── selected_bins/                        (binning pipeline only)
│   └── {mapper}/DAS_Tool_Fastas/{sample}/  per-sample DAS_Tool MAG bins
└── mag_qc/                               (binning pipeline only)
    ├── checkm2/{mapper}/{sample}/        CheckM2 quality reports
    ├── gunc/{mapper}/{sample}/           GUNC chimera detection results
    ├── gtdbtk/{mapper}/{sample}/         GTDB-tk taxonomy outputs
    ├── mag_summary.tsv                   combined MAG table (editable — see below)
    └── renamed_mags/                     MAG FASTAs with final names
```

### MAG renaming workflow

After `make_mag_summary` completes, `output/mag_qc/mag_summary.tsv` contains one row per MAG with columns including:

| Column | Description |
|---|---|
| `mag_id` | Global sequential ID (`MAG_0001` … `MAG_N`) |
| `new_name` | Proposed filename (`MAG_0001__Roseburia_intestinalis`) — **editable** |
| `original_name` | DAS_Tool bin name |
| `original_path` | Path to source FASTA |
| `sample_id` | Sample the MAG was assembled from |
| `winning_binner` | Which binner (metabat2 / maxbin2 / concoct) DAS_Tool selected |
| `domain` … `species` | GTDB-tk taxonomy in separate columns |
| `gtdbtk_classification` | Full GTDB-tk classification string |
| `completeness` / `contamination` / `quality_score` | CheckM2 metrics |
| `gunc_clade_separation_score` / `gunc_pass` | GUNC chimera metrics |
| `total_length_bp` / `num_contigs` / `gc_percent` / `N50` | Assembly stats |

To rename MAGs: edit the `new_name` column in the table, then re-run `rename_mags`. Snakemake detects the edited file is newer than the output and re-runs automatically.

---

## Citation

If you use MAGmaker in your research, please cite:

> Sanders JG, Sprockett DD, Li Y, Mjungu D, Lonsdorf EV, Ndjango JN, Georgiev AV, Hart JA, Sanz CM, Morgan DB, Peeters M, Hahn BH, Moeller AH. Widespread extinctions of co-diversified primate gut bacterial symbionts from humans. *Nat Microbiol.* 2023 Jun;8(6):1039-1050. doi: [10.1038/s41564-023-01388-w](https://doi.org/10.1038/s41564-023-01388-w). PMID: 37169918.
