# Running the pipeline

---

## Overview

MAGmaker runs in up to three stages. The `run_magmaker.sh` wrapper chains all three automatically; alternatively, each stage can be run individually.

| Stage | Command | Output |
|---|---|---|
| 1 — Main pipeline | `snakemake` | QC reports, assemblies, read-level profiles (MetaPhlAn / Kraken2), prototype selection |
| 2 — Binning config | `snakemake generate_binning_config` | `output/config/auto_binning.txt` |
| 3 — Binning pipeline | `snakemake --snakefile Snakefile-bin` | MAGs, QC reports, taxonomy, strain heterogeneity, `mag_summary.tsv` |

The optional viral / plasmid track (geNomad + CheckV) is a separate Stage 3
target, `virus_all` — see below and [Viral / plasmid track](viral.md).

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

## Per-run flags

These are consumed by `run_magmaker.sh` and turned into `--config` overrides
for one run, without editing a tracked file:

| flag | effect |
|---|---|
| `--skip-kraken` | drop the Kraken2 / Bracken targets for this run |
| `--skip-metaphlan` | drop the MetaPhlAn target for this run |
| `--skip-cmseq` | skip CMSeq strain heterogeneity; the two SH columns in `mag_summary.tsv` become `NA` |
| `--binette` / `--das-tool` | choose the bin consolidation tool, overriding `params.consolidation.tool` |
| `--skip-gtdbtk-classify` | run GTDB-Tk `identify` + `align` only — `MSA_Percent` is populated, taxonomy columns are `NA`, pplacer never runs (32–64 GB instead of 128–320 GB) |
| `--skip-gtdbtk` | skip GTDB-Tk entirely; taxonomy and `MSA_Percent` become `NA` |
| `--binners=concoct,metabat2,maxbin2` | override the binner set for this run |

`--skip-kraken` and `--skip-metaphlan` are separate because the tools are
useful independently: an analysis that takes its taxonomy from GTDB-Tk on the
MAGs needs neither, while one that wants community profiles may still want
MetaPhlAn without paying for Kraken2's ~200 GB per-sample index read.

`--config` cannot be passed **through** `run_magmaker.sh`. Snakemake keeps only
the last `--config` on a command line, so a caller-supplied one would silently
discard the overrides the wrapper sets. The wrapper refuses it and points you
at these flags or at calling `snakemake` directly.

---

## Viral / plasmid track

Not run by `run_magmaker.sh`. Request it explicitly as a Stage 3 target:

```bash
snakemake --snakefile Snakefile-bin --profile resources/profiles/demon virus_all \
  --config binning=output/config/auto_binning.txt
```

This runs geNomad (viral + plasmid contig classification, straight from the
assemblies — no mapping needed) and CheckV (viral completeness / contamination),
and writes `output/virus/virus_summary.tsv`. It reads only the assemblies, so
it can run alongside the MAG track. See [Viral / plasmid track](viral.md).

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

This reads `output/prototype_selection/prototype_selection/selected_prototypes.yaml`
and selects the `n` representative samples set by `params.prototypes.n` in
`config.yaml`. In the generated `auto_binning.txt`, **every sample contributes
its assembly** (every assembly is binned) and **only the `n` prototype samples
contribute reads**, which are mapped to every assembly to give each one a
differential-coverage profile. See [Configuration](configuration.md) for the
`binning.txt` format.

### Stage 3 — Binning pipeline

```bash
# Through bin consolidation only (Binette or DAS_Tool)
snakemake --snakefile Snakefile-bin --cores 8 --use-conda select_bins \
  --config binning=output/config/auto_binning.txt

# Through MAG QC and the taxonomy summary table
snakemake --snakefile Snakefile-bin --cores 8 --use-conda make_mag_summary \
  --config binning=output/config/auto_binning.txt

# Rename MAGs (re-run after editing mag_summary.tsv if desired)
snakemake --snakefile Snakefile-bin --cores 8 --use-conda rename_mags \
  --config binning=output/config/auto_binning.txt

# Viral / plasmid track (optional, independent of the above)
snakemake --snakefile Snakefile-bin --cores 8 --use-conda virus_all \
  --config binning=output/config/auto_binning.txt
```

`make_mag_summary` also runs **CMSeq** strain heterogeneity (unless
`--skip-cmseq`) and, for GTDB-Tk, only the stage the run needs — see the
GTDB-Tk stages section of [Output](output.md).

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

resources:
  - kraken_slots=2      # max concurrent kraken2 jobs per workflow (each reads a ~200 GB index)

set-resources:
  metaspades:
    mem_mb: 256000
  metaphlan:
    mem_mb: 32000
  taxonomy_kraken:
    mem_mb: 220000

latency-wait: 300        # a busy shared filesystem can lag well past 60s when a run's outputs land at once
rerun-incomplete: true
keep-going: true
```

Then run with:

```bash
./run_magmaker.sh --profile resources/profiles/my_cluster
```

**Most rules size their own memory and runtime and grow the request on each
retry** (`mem_escalate` / `runtime_escalate` in `resources/snakefiles/common.smk`):
the value in the profile or `config.yaml` is the *first attempt*, a
`<rule>_max` key is the ceiling, and `retries:` in `config.yaml` controls how
many escalations happen. You normally only touch a profile entry to raise a
first-attempt value for your cluster. `run_gtdbtk` (GTDB-Tk `classify`) is the
one to watch: it requests 128 GB → 256 GB → 320 GB across attempts.

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
| *(default)* | MultiQC reports, assembly stats, sourmash plot, `selected_prototypes.yaml`, and the merged MetaPhlAn / Kraken2 targets that the `biobakery` and `profilers` lists enable |
| `generate_binning_config` | Auto-generate `output/config/auto_binning.txt` from prototype selection results |

### Binning pipeline (`Snakefile-bin`)

| Target | Description |
|---|---|
| `map_all` | All read mapping steps only |
| `bin_all` | All binning steps only (requires mapping) |
| `select_bins` | Bin consolidation — Binette or DAS_Tool (requires binning) |
| `make_mag_summary` | CheckM2 + GUNC + GTDB-Tk + CMSeq + combined summary table |
| `rename_mags` | Copy MAGs to `renamed_mags/` using names from `mag_summary.tsv` |
| `virus_all` | geNomad + CheckV → `virus_summary.tsv` (independent of the MAG track) |

---

← [Documentation home](index.md)
