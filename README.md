# MAGmaker

A Snakemake pipeline for end-to-end processing of paired-end metagenomic shotgun sequencing data. MAGmaker takes raw FASTQ files through quality control, host read removal, assembly, taxonomic profiling, and — optionally — binning to produce metagenome-assembled genomes (MAGs).

Developed by the [Moeller Lab](https://moellerlab.com) at Cornell University and Princeton University, and maintained by the [Sprockett Lab](https://www.sprockettlab.com/) at Wake Forest University School of Medicine.

---

## Pipeline overview

MAGmaker runs in two stages, each driven by a separate Snakefile. Both stages can be chained automatically using the `run_magmaker.sh` wrapper.

**Stage 1 — Main pipeline (`Snakefile`)**

```
raw reads → FastQC → fastp/Cutadapt → host removal (bowtie2) → MultiQC
         → MEGAHIT/metaSPAdes → QUAST → MultiQC
         → sourmash sketch/compare → prototype selection
         → MetaPhlAn 4 → merged abundance table
```

**Stage 2 — Binning pipeline (`Snakefile-bin`)**

```
non-host reads + assemblies → bowtie2/minimap2 → sorted BAMs
BAMs + contigs → MetaBAT2 + MaxBin2 + CONCOCT → DAS_Tool → MAGs
MAGs → CheckM2 + GUNC + GTDB-tk → mag_summary.tsv → renamed_mags/
```

A `generate_binning_config` rule bridges the two stages by reading `selected_prototypes.yaml` and automatically writing the binning configuration file.

---

## Quick start

```bash
# Clone and install
git clone https://github.com/SprockettLab/MAGmaker.git
cd MAGmaker
mamba env create -n snakemake -f resources/env/snakemake.yaml
conda activate snakemake

# Edit resources/config/config.yaml, samples.txt, and units.txt

# Run everything (all three stages) — passes all args through to snakemake
./run_magmaker.sh --profile resources/profiles/demon

# Or run interactively
./run_magmaker.sh --cores 8 --use-conda
```

See the documentation below for details on each step.

---

## Documentation

| Page | Contents |
|---|---|
| [Installation](docs/installation.md) | Requirements, conda setup, test data |
| [Configuration](docs/configuration.md) | `config.yaml`, `samples.txt`, `units.txt`, `binning.txt` |
| [Database setup](docs/databases.md) | CheckM2, GUNC, GTDB-tk, MetaPhlAn, host genome |
| [Running the pipeline](docs/running.md) | Local, SLURM (demon), `run_magmaker.sh` wrapper, stage-by-stage |
| [Output](docs/output.md) | Directory layout, MAG summary table, renaming workflow |

---

## Citation

If you use MAGmaker in your research, please cite:

> Sanders JG, Sprockett DD, Li Y, Mjungu D, Lonsdorf EV, Ndjango JN, Georgiev AV, Hart JA, Sanz CM, Morgan DB, Peeters M, Hahn BH, Moeller AH. Widespread extinctions of co-diversified primate gut bacterial symbionts from humans. *Nat Microbiol.* 2023 Jun;8(6):1039-1050. doi: [10.1038/s41564-023-01388-w](https://doi.org/10.1038/s41564-023-01388-w). PMID: 37169918.
