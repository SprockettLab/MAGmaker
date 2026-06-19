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
# Demon cluster (SLURM)
./run_magmaker.sh --profile resources/profiles/demon

# Local / interactive
./run_magmaker.sh --cores 8 --use-conda

# Dry run (checks all three DAGs without executing)
./run_magmaker.sh --profile resources/profiles/demon -n
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

## SLURM (demon cluster)

A Snakemake profile for demon is at `resources/profiles/demon/`. It:
- Submits to the `defq` partition
- Uses 8 GB RAM and 120-minute walltime as per-job defaults
- Passes `--use-conda` automatically
- Sets the shared conda prefix (`/isilon/.../snakemake_envs/`)

```bash
conda activate /isilon/datalake/sprockett_lab/original/WF00SprockettLab/envs/snakemake_envs/<hash>

# Full pipeline
./run_magmaker.sh --profile resources/profiles/demon

# Or individual stages
snakemake --profile resources/profiles/demon
snakemake --profile resources/profiles/demon generate_binning_config
snakemake --snakefile Snakefile-bin --profile resources/profiles/demon rename_mags \
  --config binning=output/config/auto_binning.txt
```

To override resources for a specific rule:

```bash
snakemake --profile resources/profiles/demon \
  --set-resources megahit:mem_mb=256000 megahit:runtime=2880
```

SLURM job logs go to `.snakemake/slurm_logs/{rule}/`. Rule-level logs (tool stderr/stdout) go to `output/logs/{rule}/`.

**demon `defq` partition limits:**
- Max 15 nodes per job
- Max walltime: 90 days
- Default memory: 4 GB/CPU (overridden explicitly by the profile)

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
