# Running the pipeline

---

## Overview

MAGmaker runs in up to three stages. The `run_magmaker.sh` wrapper chains all three automatically; alternatively, each stage can be run individually.

| Stage | Command | Output |
|---|---|---|
| 1 — Main pipeline | `snakemake` | QC reports, assemblies, prototype selection, MetaPhlAn profiles |
| 2 — Binning config | `snakemake generate_binning_config` | `output/config/auto_binning.txt` |
| 3 — Binning pipeline | `snakemake --snakefile Snakefile-bin` | MAGs, QC reports, taxonomy, `mag_summary.tsv` |

---

## `run_magmaker.sh` — full pipeline in one command

The `run_magmaker.sh` wrapper at the repository root chains all three stages automatically. All arguments are passed through to every `snakemake` invocation, so `--profile`, `--cores`, `--use-conda`, `-n`, and any other Snakemake flags work as expected.

```bash
# On an HPC cluster with a Snakemake profile
./run_magmaker.sh --profile resources/profiles/your_cluster

# Local / interactive
./run_magmaker.sh --cores 8 --use-conda

# Dry run (checks all three DAGs without executing)
./run_magmaker.sh --cores 8 --use-conda -n
```

**How it works:**

1. Runs the main pipeline (`snakemake "$@"`)
2. Runs `generate_binning_config` to produce `output/config/auto_binning.txt`
3. Runs the binning pipeline through `rename_mags` (`snakemake --snakefile Snakefile-bin "$@" rename_mags --config binning=output/config/auto_binning.txt`)

**Dry-run behavior:** On a dry run (`-n`), stages 1 and 2 print their DAGs but `auto_binning.txt` is never written (nothing executes). Stage 3 requires `auto_binning.txt` to exist before Snakemake can parse `Snakefile-bin`, so the wrapper detects the missing file and prints a message instead of crashing. Run without `-n` to actually execute stages 1 and 2 first; subsequent runs (including further dry runs) will find the file and show the full stage 3 DAG.

**Stopping after the MAG summary table** to review/edit `mag_summary.tsv` before renaming: run the stages manually (see below).

---

## Running stages individually

If you prefer more control, each stage can be invoked directly.

### Stage 1 — Main pipeline

```bash
conda activate snakemake

# Dry run
snakemake --cores 8 --use-conda -n

# Execute
snakemake --cores 8 --use-conda
```

### Stage 2 — Binning config

After stage 1 completes, run `generate_binning_config` to write `output/config/auto_binning.txt`:

```bash
snakemake --cores 1 --use-conda generate_binning_config
```

This reads `output/prototype_selection/prototype_selection/selected_prototypes.yaml` and selects the `n` representative samples specified by `params.prototypes.n` in `config.yaml`. All samples contribute reads; only prototype samples contribute assemblies. See [Configuration](configuration.md) for details on `binning.txt` format.

### Stage 3 — Binning pipeline

```bash
# Through DAS_Tool bin selection only
snakemake --snakefile Snakefile-bin --cores 8 --use-conda select_bins \
  --config binning=output/config/auto_binning.txt

# Through MAG QC and taxonomy summary table
snakemake --snakefile Snakefile-bin --cores 8 --use-conda make_mag_summary \
  --config binning=output/config/auto_binning.txt

# Rename MAGs (re-run after editing mag_summary.tsv if desired)
snakemake --snakefile Snakefile-bin --cores 8 --use-conda rename_mags \
  --config binning=output/config/auto_binning.txt
```

---

## Running on an HPC cluster (SLURM)

Snakemake submits each rule as a separate SLURM job via a [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles). A profile is a directory containing a `config.yaml` that sets the executor, default resources, and any per-rule overrides.

An example profile for a generic SLURM cluster:

```yaml
# resources/profiles/my_cluster/config.yaml
executor: slurm

default-resources:
  slurm_partition: normal
  mem_mb: 8000
  runtime: 120        # minutes

jobs: 100
use-conda: true
conda-prefix: /shared/path/to/conda_envs/

set-resources:
  megahit:
    mem_mb: 256000
    runtime: 1440
  metaphlan:
    mem_mb: 32000
    runtime: 240
  taxonomy_kraken:
    mem_mb: 220000
    runtime: 240

latency-wait: 60
rerun-incomplete: true
keep-going: true
```

Then run with:

```bash
./run_magmaker.sh --profile resources/profiles/my_cluster
```

To override resources for a specific rule at runtime:

```bash
snakemake --profile resources/profiles/my_cluster \
  --set-resources megahit:mem_mb=512000 megahit:runtime=2880
```

SLURM job logs go to `.snakemake/slurm_logs/{rule}/`. Rule-level logs (tool stderr/stdout) go to `output/logs/{rule}/`.

> **WFUSM users:** See [Running on DEMON](demon.md) for a ready-to-use setup with all databases and environments pre-configured.

---

## Available top-level targets

### Main pipeline (`Snakefile`)

| Target | Description |
|---|---|
| *(default)* | MultiQC reports, assembly stats, sourmash plot, selected_prototypes.yaml, merged MetaPhlAn table |
| `generate_binning_config` | Auto-generate `output/config/auto_binning.txt` from prototype selection results |

### Binning pipeline (`Snakefile-bin`)

| Target | Description |
|---|---|
| `map_all` | All read mapping steps only |
| `bin_all` | All binning steps only (requires mapping) |
| `select_bins` | DAS_Tool bin dereplication (requires binning) |
| `make_mag_summary` | CheckM2 + GUNC + GTDB-tk + combined summary table |
| `rename_mags` | Copy MAGs to `renamed_mags/` using names from `mag_summary.tsv` |

---

← [Back to README](../README.md)
