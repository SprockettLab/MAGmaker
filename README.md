# MAGmaker

A Snakemake pipeline for end-to-end processing of metagenomic shotgun sequencing data (paired-end or single-end). MAGmaker takes raw FASTQ files through quality control, host read removal, assembly, and taxonomic profiling, then — in a second stage — through mapping, binning, bin consolidation, and MAG QC to produce metagenome-assembled genomes (MAGs). An optional viral/plasmid discovery track (geNomad + CheckV) runs alongside the MAG track.

Developed by the [Moeller Lab](https://moellerlab.com) at Cornell University and Princeton University, and maintained by the [Sprockett Lab](https://www.sprockettlab.com/) at Wake Forest University School of Medicine.

---

## Pipeline overview

MAGmaker runs in two stages, each driven by a separate Snakefile. Both stages can be chained automatically using the `run_magmaker.sh` wrapper.

**Stage 1 — Main pipeline (`Snakefile`)**

```
raw reads → FastQC → fastp / Cutadapt → host removal (Bowtie2) → FastQC → MultiQC
         → MEGAHIT / metaSPAdes → QUAST → MultiQC
         → MetaPhlAn 4 → merged abundance table          (optional)
         → Kraken 2 → Bracken → per-sample profiles       (optional)
         → sourmash sketch / compare → prototype selection
```

**Stage 2 — Binning pipeline (`Snakefile-bin`)**

```
contigs + prototype reads → minimap2 / Bowtie2 → sorted BAMs (differential coverage)
BAMs + contigs → MetaBAT2 + MaxBin2 + CONCOCT + SemiBin2
              → Binette (default) or DAS_Tool → one MAG set per sample
MAGs → CheckM2 + GUNC + GTDB-Tk + CMSeq → mag_summary.tsv → renamed_mags/

Viral / plasmid track (opt-in, `virus_all` target):
contigs → geNomad → CheckV → virus_summary.tsv
```

A `generate_binning_config` rule bridges the two stages: it reads
`selected_prototypes.yaml`, picks `params.prototypes.n` representative samples,
and writes `auto_binning.txt`. **Every sample's assembly is binned; the
prototype samples supply the reads that form the differential-coverage basis.**

---

## Quick start

```bash
# Clone and install
git clone https://github.com/SprockettLab/MAGmaker.git
cd MAGmaker
mamba env create -n snakemake -f resources/env/snakemake.yaml
conda activate snakemake

# Edit resources/config/config.yaml and metadata.txt

# Run Stage 1 + Stage 2 through MAG renaming — passes all args through to snakemake
./run_magmaker.sh --cores 8 --use-conda

# On an HPC cluster with a Snakemake SLURM profile
./run_magmaker.sh --profile resources/profiles/your_cluster
```

`run_magmaker.sh` runs Stage 1 and Stage 2 through MAG renaming. Read-level
profiling (MetaPhlAn, Kraken2) is toggled by the `biobakery` / `profilers`
config lists or per-run flags, and the viral track is opt-in — see
[Running the pipeline](docs/running.md).

> **WFUSM users on DEMON:** See [Running on DEMON](docs/demon.md) for a self-contained setup guide with all databases and environments pre-configured.

See the documentation below for details on each step.

---

## Documentation

| Page | Contents |
|---|---|
| [Pipeline overview](docs/pipeline.md) | What each stage does and why, tool choices, the two-Snakefile structure |
| [Installation](docs/installation.md) | Requirements, conda setup, test data |
| [Configuration](docs/configuration.md) | `config.yaml`, `metadata.txt`, `binning.txt` |
| [Database setup](docs/databases.md) | CheckM2, GUNC, GTDB-Tk, MetaPhlAn, Kraken2, geNomad, CheckV, host genome |
| [Running the pipeline](docs/running.md) | Local, HPC/SLURM, `run_magmaker.sh` wrapper, per-run flags, stage-by-stage |
| [Output](docs/output.md) | Directory layout, MAG summary table, renaming workflow, GTDB-Tk stages |
| [Viral / plasmid track](docs/viral.md) | geNomad + CheckV, `virus_summary.tsv` |
| [Running on DEMON (WFUSM)](docs/demon.md) | DEMON cluster setup, pre-configured databases and environments |

---

## Citation

If you use MAGmaker in your research, please cite:

> Sanders JG, Sprockett DD, Li Y, Mjungu D, Lonsdorf EV, Ndjango JN, Georgiev AV, Hart JA, Sanz CM, Morgan DB, Peeters M, Hahn BH, Moeller AH. Widespread extinctions of co-diversified primate gut bacterial symbionts from humans. *Nat Microbiol.* 2023 Jun;8(6):1039-1050. doi: [10.1038/s41564-023-01388-w](https://doi.org/10.1038/s41564-023-01388-w). PMID: 37169918.
