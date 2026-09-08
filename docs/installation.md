# Installation

## Requirements

- Linux or macOS (the binning stage's CONCOCT env is Linux-first; see the note below)
- [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/en/latest/)
- ~30–40 GB disk space for conda environments (built automatically on first run)
- A host genome FASTA for host read removal (e.g., human GRCh38, mouse GRCm39)
- Reference databases are separate and much larger — GTDB-Tk ~100 GB, Kraken2
  Standard ~200 GB, geNomad ~1.5 GB, CheckV ~1.5 GB, CheckM2 / GUNC / MetaPhlAn
  a few GB each. See [Database setup](databases.md).

All tool dependencies (FastQC, fastp, MEGAHIT, MetaPhlAn 4, sourmash, MetaBAT2,
MaxBin2, CONCOCT, SemiBin2, DAS_Tool, Binette, CheckM2, GUNC, GTDB-Tk, CMSeq,
geNomad, CheckV, …) are installed automatically by Snakemake into isolated
conda environments on first run.

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
- **Single-end / mixed-layout examples:** `resources/test/metadata_single_end.txt` and `resources/test/metadata_mixed_layout.txt` run against the same test reads; `python3 resources/test/test_metadata_layout.py` exercises the metadata loader.

The default `metadata.txt` and `binning.txt` point to this test data. To run a test:

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
| `fastp.yaml` | `qc.smk` — default trimmer (paired + single-end) |
| `cutadapt.yaml` | `qc.smk` — non-default trimmer |
| `fastqc.yaml` | `qc.smk` — all FastQC rules |
| `qc.yaml` | `qc.smk` — host filter and index build (bowtie2, samtools, pigz) |
| `assemble.yaml` | `assemble.smk` — MEGAHIT, metaSPAdes, QUAST |
| `mapping.yaml` | `mapping.smk`, and `cmseq.smk`'s mapping step — minimap2, bowtie2, samtools |
| `prototype_selection.yaml` | `prototype_selection.smk` — sourmash, scikit-bio |
| `profile.yaml` | `profile.smk` (MetaPhlAn 4, Kraken2, Bracken, Krona) and `cmseq.smk`'s `poly.py` |
| `binning.yaml` | `binning.smk` — MetaBAT2, MaxBin2, FragGeneScan |
| `concoct_linux.yaml` | `binning.smk` — CONCOCT (Linux); on macOS point the `conda:` at `concoct_osx.yaml` |
| `semibin.yaml` | `binning.smk` — SemiBin2 |
| `selected_bins.yaml` | `selected_bins.smk` — DAS_Tool, `Fasta_to_Contig2Bin` |
| `binette.yaml` | `selected_bins.smk` — Binette (default consolidation tool) |
| `checkm2.yaml` | `mag_qc.smk` — CheckM2 (separate env; Python version gap) |
| `gunc.yaml` | `mag_qc.smk` — GUNC (separate env; conflicts with CheckM2) |
| `gtdbtk.yaml` | `mag_qc.smk` — GTDB-Tk (separate env; dependency conflicts) |
| `mag_qc.yaml` | `mag_qc.smk`, `cmseq.smk`, `virus.smk` summary/aggregation scripts (pandas) |
| `genomad.yaml` | `virus.smk` — geNomad |
| `checkv.yaml` | `virus.smk` — CheckV |

`resources/env/` also holds a few unwired alternates and single-tool
variants (`concoct_osx.yaml`, `bowtie2.yaml`, `metaphlan.yaml`,
`sourmash.yaml`, `binning_concoct_*.yaml`) that no rule references by default.

> **Note:** `tbb=2020.2` is pinned in `qc.yaml` and `mapping.yaml` because newer TBB versions break bowtie2 on some systems. Do not remove this pin without testing.

---

← [Documentation home](index.md)
