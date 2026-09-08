# MAGmaker

A Snakemake pipeline for end-to-end processing of metagenomic shotgun
sequencing data (paired-end or single-end). MAGmaker takes raw FASTQ files
through quality control, host read removal, assembly, and taxonomic profiling,
then — in a second stage — through mapping, binning, bin consolidation, and MAG
QC to produce metagenome-assembled genomes (MAGs). An optional viral / plasmid
discovery track (geNomad + CheckV) runs alongside the MAG track.

Developed by the [Moeller Lab](https://moellerlab.com) at Cornell University
and Princeton University, and maintained by the
[Sprockett Lab](https://www.sprockettlab.com/) at Wake Forest University School
of Medicine.

---

## Quick start

```bash
git clone https://github.com/SprockettLab/MAGmaker.git
cd MAGmaker
mamba env create -n snakemake -f resources/env/snakemake.yaml
conda activate snakemake

# Edit resources/config/config.yaml and metadata.txt

# Run Stage 1 + Stage 2 through MAG renaming
./run_magmaker.sh --cores 8 --use-conda

# On an HPC cluster with a Snakemake SLURM profile
./run_magmaker.sh --profile resources/profiles/your_cluster
```

WFUSM users on DEMON: start at [Running on DEMON](demon.md) — all databases and
environments are pre-configured.

---

## Where to go next

| Page | Contents |
|---|---|
| [Pipeline overview](pipeline.md) | What each stage does and why, tool choices, the two-Snakefile structure |
| [Installation](installation.md) | Requirements, conda setup, test data |
| [Configuration](configuration.md) | `config.yaml`, `metadata.txt`, `binning.txt` |
| [Database setup](databases.md) | CheckM2, GUNC, GTDB-Tk, MetaPhlAn, Kraken2, geNomad, CheckV, host genome |
| [Running the pipeline](running.md) | Local, HPC/SLURM, `run_magmaker.sh`, per-run flags, stage-by-stage |
| [Output](output.md) | Directory layout, MAG summary table, renaming workflow, GTDB-Tk stages |
| [Viral / plasmid track](viral.md) | geNomad + CheckV, `virus_summary.tsv` |
| [Running on DEMON (WFUSM)](demon.md) | DEMON cluster setup, pre-configured databases and environments |

---

## Citation

> Sanders JG, Sprockett DD, Li Y, Mjungu D, Lonsdorf EV, Ndjango JN,
> Georgiev AV, Hart JA, Sanz CM, Morgan DB, Peeters M, Hahn BH, Moeller AH.
> Widespread extinctions of co-diversified primate gut bacterial symbionts
> from humans. *Nat Microbiol.* 2023 Jun;8(6):1039-1050.
> doi: [10.1038/s41564-023-01388-w](https://doi.org/10.1038/s41564-023-01388-w).
> PMID: 37169918.
