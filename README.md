# MAGmaker

A Snakemake pipeline for end-to-end processing of paired-end metagenomic shotgun sequencing data. MAGmaker takes raw FASTQ files through quality control, host read removal, assembly, taxonomic profiling, and — optionally — binning to produce metagenome-assembled genomes (MAGs).

Developed in the [Moeller Lab](https://moellerlab.com) at Princeton University and maintained by the [Sprockett Lab](https://www.sprockettlab.com/) at Wake Forest University School of Medicine.

---

## Pipeline overview

MAGmaker has two modes, each driven by a separate Snakefile.

### Main pipeline (`Snakefile`)

Runs QC through taxonomic profiling and prototype (representative sample) selection.

```
raw reads (FASTQ)
  ├── FastQC (pre-trim)
  ├── Cutadapt              adapter and quality trimming
  ├── FastQC (post-trim)
  ├── merge across lanes    cat reads per sample across sequencing units
  ├── bowtie2               host read removal
  ├── FastQC (post-host)
  └── MultiQC               → output/qc/multiqc/multiqc.html

non-host reads
  ├── MEGAHIT / metaSPAdes  de novo assembly (one or both, configurable)
  ├── QUAST                 assembly quality assessment
  └── MultiQC               → output/assemble/multiqc_assemble/multiqc.html

non-host reads (similarity)
  ├── sourmash sketch       MinHash sketches per sample
  ├── sourmash compare      pairwise distance matrix
  ├── sourmash plot         → output/prototype_selection/sourmash_plot/
  └── prototype selection   → output/prototype_selection/prototype_selection/selected_prototypes.yaml

non-host reads (taxonomy)
  ├── MetaPhlAn 4           marker-gene taxonomic profiling
  └── merge tables          → output/profile/metaphlan/merged_abundance_table.txt
```

### Binning pipeline (`Snakefile-bin`)

Extends the main pipeline with read mapping, binning, and bin refinement.

```
non-host reads + assemblies
  ├── bowtie2 / minimap2    map reads to contigs (configurable)
  └── samtools              sort and index BAMs

sorted BAMs + contigs
  ├── MetaBAT2              coverage-based binning
  ├── MaxBin2               coverage-based binning
  └── CONCOCT               composition + coverage binning

all bins
  └── DAS_Tool              bin refinement and dereplication → MAGs
```

---

## Requirements

- Linux or macOS
- [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/en/latest/)
- ~20 GB disk space for conda environments (built automatically on first run)
- A host genome FASTA for filtering (e.g., human GRCh38, mouse GRCm39)
- A [MetaPhlAn 4 database](https://huttenhower.sph.harvard.edu/metaphlan/) if running taxonomic profiling

All other tool dependencies (FastQC, Cutadapt, MEGAHIT, MetaPhlAn, sourmash, etc.) are installed automatically by Snakemake into isolated conda environments.

---

## Installation

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

All other environments (qc, assemble, profile, etc.) are created automatically the first time a rule runs, via `--use-conda`.

---

## Configuration

All input files and parameters live in `resources/config/`. The `config.yaml` at the repository root is a symlink to `resources/config/config.yaml`.

### `samples.txt`

Tab-separated. One row per sample. The first column must be named `Sample`. Additional columns are carried through but not used by the pipeline.

```
Sample    Subject    Timepoint
amy       A          T1
bob       B          T1
```

### `units.txt`

Tab-separated. One row per sequencing unit (e.g., lane or run). Columns: `Sample`, `Unit`, `R1`, `R2`. Multiple units per sample are automatically concatenated before assembly.

```
Sample    Unit      R1                        R2
amy       Run_1     /path/to/amy_R1.fastq.gz  /path/to/amy_R2.fastq.gz
bob       Run_1     /path/to/bob_R1.fastq.gz  /path/to/bob_R2.fastq.gz
```

### `config.yaml`

Key settings to review before running:

| Parameter | Description |
|---|---|
| `assemblers` | Which assembler(s) to use: `metaspades`, `megahit`, or both |
| `host_filter.genome` | Path to host genome FASTA (used to build bowtie2 index) |
| `host_filter.db_dir` | Directory for the bowtie2 host index |
| `params.metaphlan.db_path` | Path to MetaPhlAn 4 database directory |
| `params.metaphlan.db_name` | MetaPhlAn 4 database name (e.g., `mpa_vJan25_CHOCOPhlAnSGB_202503`) |
| `threads.*` | Per-rule thread counts — scale to your hardware |

### `binning.txt` (binning pipeline only)

Tab-separated. Defines which reads get mapped to which assemblies. Columns: `Sample`, `Contigs`, `Read_Groups`, `Contig_Groups`. Samples sharing a group label in both `Read_Groups` and `Contig_Groups` are paired for mapping. See `resources/config/binning.txt` for an example.

---

## Running the pipeline

### Local / interactive

```bash
conda activate snakemake

# Dry run — check DAG without executing
snakemake --cores 8 --use-conda -n

# Main pipeline
snakemake --cores 8 --use-conda

# Binning pipeline
snakemake --snakefile Snakefile-bin --cores 8 --use-conda
```

### SLURM (demon cluster)

A SLURM profile for the demon cluster is included at `resources/profiles/demon/`. It submits jobs to the `defq` partition with 8 GB RAM and 120-minute walltime as defaults, and automatically manages conda environments in a shared location.

```bash
conda activate snakemake

# Main pipeline
snakemake --profile resources/profiles/demon

# Binning pipeline
snakemake --snakefile Snakefile-bin --profile resources/profiles/demon
```

To override resources for a specific rule:

```bash
snakemake --profile resources/profiles/demon \
  --set-resources megahit:mem_mb=256000 megahit:runtime=2880
```

SLURM logs are written to `.snakemake/slurm_logs/{rule}/`. Rule logs (stderr/stdout from the tools themselves) go to `output/logs/{rule}/`.

> **Note:** Assembly rules (MEGAHIT, metaSPAdes) request 256 GB RAM by default. Adjust `mem_mb` in `config.yaml` if your nodes have different memory limits.

### Running on test data

The repository includes a small test dataset (four samples: amy, bob, carl, diane) with pre-built bowtie2 index files. The default `samples.txt`, `units.txt`, and `config.yaml` point to this test data. To run a quick test:

```bash
snakemake --cores 4 --use-conda -n   # dry run first
snakemake --cores 4 --use-conda
```

---

## Output

```
output/
├── qc/
│   ├── multiqc/multiqc.html              QC report (pre-trim through host removal)
│   └── host_filter/nonhost/              host-filtered FASTQ files
├── assemble/
│   ├── megahit/{sample}/                 MEGAHIT assemblies
│   ├── metaspades/{sample}/              metaSPAdes assemblies
│   └── multiqc_assemble/multiqc.html     assembly QC report
├── prototype_selection/
│   ├── sourmash_plot/                    pairwise similarity heatmap
│   └── prototype_selection/
│       └── selected_prototypes.yaml      representative sample selection
├── profile/
│   └── metaphlan/
│       ├── profiles/{sample}.txt         per-sample MetaPhlAn profiles
│       └── merged_abundance_table.txt    merged taxonomy table (all samples)
└── bins/                                 (binning pipeline only)
    └── das_tool/{sample}/                refined MAG bins
```

---

## Citation

If you use MAGmaker in your research, please cite:

> Sanders JG, Sprockett DD, Li Y, Mjungu D, Lonsdorf EV, Ndjango JN, Georgiev AV, Hart JA, Sanz CM, Morgan DB, Peeters M, Hahn BH, Moeller AH. Widespread extinctions of co-diversified primate gut bacterial symbionts from humans. *Nat Microbiol.* 2023 Jun;8(6):1039-1050. doi: [10.1038/s41564-023-01388-w](https://doi.org/10.1038/s41564-023-01388-w). PMID: 37169918.
