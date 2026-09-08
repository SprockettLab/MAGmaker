# Running MAGmaker on DEMON (WFUSM)

This page covers everything needed to run MAGmaker on the Wake Forest University School of Medicine DEMON cluster. If you have access to DEMON, start here — all databases and environments are already set up.

---

## Prerequisites

- An active DEMON account with access to the `defq` partition
- Access to the Sprockett Lab Isilon share (`/isilon/datalake/sprockett_lab/`)
- Git and conda available (both are available cluster-wide)

---

## Getting started

Clone MAGmaker to your project directory on the Isilon share:

```bash
cd /isilon/datalake/sprockett_lab/original/WF00SprockettLab/<your_project>/
git clone https://github.com/SprockettLab/MAGmaker.git
cd MAGmaker
```

Activate the shared Snakemake environment:

```bash
conda activate /isilon/datalake/sprockett_lab/original/WF00SprockettLab/envs/envs/snakemake
```

All pipeline tool environments (FastQC, MEGAHIT, MetaPhlAn, GTDB-tk, etc.) are pre-built and stored in the shared conda prefix — they do not need to be built on first run.

---

## Configuration

The default `config.yaml` is already configured for DEMON. All databases are pre-downloaded and their paths are set:

| Database | Path |
|---|---|
| CheckM2 | `/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/checkm2/CheckM2_database/uniref100.KO.1.dmnd` |
| GUNC | `/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/gunc/gunc_db_progenomes2.1.dmnd` |
| GTDB-Tk (release 232) | `/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/gtdbtk/release232/` |
| MetaPhlAn 4 | `/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/metaphlan/` |
| Kraken2 Standard | `/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/kraken2/Kraken2_db_Standard/` |
| geNomad | `/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/genomad/genomad_db/` |
| CheckV | `/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/checkv/checkv-db-v1.5/` |
| Human GRCh38 (bowtie2 index) | `/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/bt2/human_GCA_000001405.29_GRCh38.p14/` |
| Conda environments | `/isilon/datalake/sprockett_lab/original/WF00SprockettLab/envs/snakemake_envs/` |

The only things you need to change in `config.yaml` before running are:

1. **`host_filter.genome`** — path to your host genome FASTA (human GRCh38 bowtie2 index is already built; set `db_dir` to the path above)
2. **`metadata`** — path to your `metadata.txt` file

Everything else — CheckM2, GUNC, GTDB-Tk, MetaPhlAn, Kraken2, geNomad, CheckV
paths — is already set for DEMON.

See [Configuration](configuration.md) for the format of these files.

---

## Running the pipeline

The DEMON SLURM profile is at `resources/profiles/demon/`. It:

- Submits jobs to the `defq` partition
- Uses 8 GB RAM and 120-minute walltime as per-job defaults (overridden for memory-intensive rules)
- Passes `--use-conda` and the shared conda prefix automatically

### Full pipeline (recommended)

```bash
conda activate /isilon/datalake/sprockett_lab/original/WF00SprockettLab/envs/envs/snakemake

# Dry run first — checks the full DAG without submitting jobs
./run_magmaker.sh --profile resources/profiles/demon -n

# Execute
./run_magmaker.sh --profile resources/profiles/demon
```

### Stage by stage

```bash
# Stage 1 — Main pipeline
snakemake --profile resources/profiles/demon

# Stage 2 — Generate binning config
snakemake --profile resources/profiles/demon generate_binning_config

# Stage 3 — Binning pipeline through MAG renaming
snakemake --snakefile Snakefile-bin --profile resources/profiles/demon rename_mags \
  --config binning=output/config/auto_binning.txt

# Optional — viral / plasmid track
snakemake --snakefile Snakefile-bin --profile resources/profiles/demon virus_all \
  --config binning=output/config/auto_binning.txt
```

### Monitoring jobs

```bash
# Check running jobs
squeue -u $USER

# Watch the Snakemake log (run from the MAGmaker directory)
tail -f logs/$(ls -t logs/ | head -1)

# Cancel all your jobs if needed
scancel -u $USER
```

---

## SLURM profile defaults

The demon profile (`resources/profiles/demon/config.yaml`) submits to `defq`
with an 8 GB / 120-minute default, `--use-conda` with the shared conda prefix,
`latency-wait: 300`, and `resources: [kraken_slots=2]` (at most two Kraken2
jobs at once — each reads a ~200 GB index).

First-attempt memory overrides it sets (`set-resources`):

| Rule | First-attempt memory |
|---|---|
| `host_filter` | 16 GB |
| `metaspades` | 256 GB |
| `metaphlan` | 32 GB |
| `taxonomy_kraken` | 220 GB |
| `make_concoct_coverage_table`, `make_metabat2_coverage_table` | 64 GB (many-BAM references) |
| `run_DAS_Tool`, `run_gunc`, `run_genomad` | 32 GB (DIAMOND / marker search) |
| `run_checkv` | 16 GB |
| `run_gtdbtk` (classify) | 128 GB → 256 GB → 320 GB (escalates on retry) |

Runtimes are no longer set in the profile — as of 2026-08-19 they moved to
`runtime_escalate()` calls on the rules themselves, which grow the wall-clock
request on each retry (see [Running the pipeline](running.md#running-on-an-hpc-cluster-slurm)).
Several rules that use DIAMOND also carry `retries:` in `config.yaml` because
their parallel plumbing fails transiently under high cluster concurrency; a
fresh attempt lands in a lower-contention window and clears it.

To further override resources for a single rule:

```bash
snakemake --profile resources/profiles/demon \
  --set-resources megahit:mem_mb=512000 megahit:runtime=2880
```

### `defq` partition limits

- Max walltime: 90 days
- Default memory: 4 GB/CPU (overridden explicitly by the profile)
- Snakemake submits up to 100 jobs simultaneously

---

## Private analysis files

Project-specific `metadata.txt` and `binning.txt` files should be named with a project prefix (e.g., `MyProject_metadata.txt`) and stored in `resources/config/`. The `.gitignore` excludes `*_metadata.txt`, `*_binning.txt` and `*_config.yaml` to prevent accidental commits of data containing patient or project identifiers.

---

## Logs

- **Snakemake workflow logs:** `logs/<timestamp>.snakemake.log`
- **Rule-level tool logs:** `output/logs/{rule}/{sample}.log`
- **SLURM job logs:** `.snakemake/slurm_logs/{rule}/{sample}/`

---

← [Documentation home](index.md)
