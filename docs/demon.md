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
| GTDB-tk (release 232) | `/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/gtdbtk/release232/` |
| MetaPhlAn 4 | `/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/metaphlan/` |
| Kraken2 Standard | `/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/kraken2/Kraken2_db_Standard/` |
| Human GRCh38 (bowtie2 index) | `/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/bt2/human_GCA_000001405.29_GRCh38.p14/` |
| Conda environments | `/isilon/datalake/sprockett_lab/original/WF00SprockettLab/envs/snakemake_envs/` |

The only things you need to change in `config.yaml` before running are:

1. **`host_filter.genome`** — path to your host genome FASTA (human GRCh38 bowtie2 index is already built; set `db_dir` to the path above)
2. **`samples`** — path to your `samples.txt` file
3. **`units`** — path to your `units.txt` file

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

The demon profile (`resources/profiles/demon/config.yaml`) sets per-rule memory overrides for rules that need more than the 8 GB default:

| Rule | Memory | Walltime |
|---|---|---|
| `metaphlan` | 32 GB | 4 hours |
| `taxonomy_kraken` | 220 GB | 4 hours |
| all other rules | 8 GB | 2 hours |

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

Project-specific `samples.txt`, `units.txt`, and `binning.txt` files should be named with a project prefix (e.g., `MyProject_samples.txt`) and stored in `resources/config/`. The `.gitignore` excludes `*_samples.txt`, `*_units.txt`, and `*_binning.txt` to prevent accidental commits of data containing patient or project identifiers.

---

## Logs

- **Snakemake workflow logs:** `logs/<timestamp>.snakemake.log`
- **Rule-level tool logs:** `output/logs/{rule}/{sample}.log`
- **SLURM job logs:** `.snakemake/slurm_logs/{rule}/{sample}/`

---

← [Back to README](../README.md)
