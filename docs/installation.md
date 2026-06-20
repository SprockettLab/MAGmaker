# Installation

## Requirements

- Linux or macOS
- [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/en/latest/)
- ~20 GB disk space for conda environments (built automatically on first run)
- A host genome FASTA for host read removal (e.g., human GRCh38, mouse GRCm39)

All tool dependencies (FastQC, fastp, MEGAHIT, MetaPhlAn 4, sourmash, CheckM2, GUNC, GTDB-tk, etc.) are installed automatically by Snakemake into isolated conda environments on first run.

---

## Setup

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

All other environments (qc, assemble, profile, binning, mag_qc, etc.) are created automatically the first time a rule that needs them runs. Pass `--use-conda` when running locally; the demon profile handles this automatically.

---

## Test data

The repository includes a minimal test dataset to verify the pipeline is installed and configured correctly:

- **Reads:** `resources/test/test_reads/{John,Paul,George,Ringo}_{R1,R2}.fastq.gz` — tiny FASTQ files
- **Host index:** `resources/test/test_dbs/GCA_000001635.9.*` — pre-built bowtie2 index for mouse chromosome 1 (tiny subset, sufficient for testing host filter logic)

The default `samples.txt`, `units.txt`, and `binning.txt` point to this test data. To run a test:

```bash
conda activate snakemake

# Dry run first
snakemake --cores 4 --use-conda -n

# Execute
snakemake --cores 4 --use-conda
```

To use the test host index, set in `config.yaml`:

```yaml
host_filter:
  db_dir: resources/test/test_dbs
  genome: resources/db/bt2/GCA_000001635.9.fna
```

---

## Conda environments

Each pipeline module has its own environment YAML in `resources/env/`:

| File | Used by |
|---|---|
| `snakemake.yaml` | top-level (install this first) |
| `fastp.yaml` | `qc.smk` — default trimmer |
| `qc.yaml` | `qc.smk` — FastQC, bowtie2, samtools |
| `assemble.yaml` | `assemble.smk` — MEGAHIT, metaSPAdes, QUAST |
| `mapping.yaml` | `mapping.smk` — bowtie2, minimap2, samtools |
| `binning.yaml` | `binning.smk` — MetaBAT2, MaxBin2, FragGeneScan |
| `concoct_linux.yaml` | `binning.smk` — CONCOCT (Linux only) |
| `selected_bins.yaml` | `selected_bins.smk` — DAS_Tool |
| `mag_qc.yaml` | `mag_qc.smk` — CheckM2, GUNC |
| `gtdbtk.yaml` | `mag_qc.smk` — GTDB-tk (separate env due to dependency conflicts) |
| `profile.yaml` | `profile.smk` — MetaPhlAn 4, Kraken2, Bracken |
| `prototype_selection.yaml` | `prototype_selection.smk` — sourmash, scikit-bio |

> **Note:** `tbb=2020.2` is pinned in `qc.yaml` and `mapping.yaml` because newer TBB versions break bowtie2 on some systems. Do not remove this pin without testing.

---

← [Back to README](../README.md)
