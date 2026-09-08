# Database setup

MAGmaker needs several reference databases. Only the ones for the stages you
run are required: read-level profiling needs MetaPhlAn and/or Kraken2, MAG QC
needs CheckM2 + GUNC + GTDB-Tk, and the viral track needs geNomad + CheckV.

| Database | For | Size (approx.) | `config.yaml` key |
|---|---|---|---|
| Host genome FASTA | host read removal | varies | `host_filter.genome` |
| MetaPhlAn 4 | read-level profiling | ~20 GB | `params.metaphlan.db_path` |
| Kraken2 + Bracken | read-level profiling | ~200 GB (Standard) | `params.kraken2.db` / `bracken-db` |
| CheckM2 | MAG completeness / contamination | ~3 GB | `params.checkm2.db_path` |
| GUNC | MAG chimerism | ~13 GB | `params.gunc.db_path` |
| GTDB-Tk | MAG taxonomy | ~100 GB | `params.gtdbtk.db_path` |
| geNomad | viral / plasmid discovery | ~1.5 GB | `params.genomad.db_path` |
| CheckV | viral completeness / contamination | ~1.5 GB | `params.checkv.db_path` |

---

## CheckM2

CheckM2 assesses MAG completeness and contamination using a diamond protein database.

If `params.checkm2.db_path` is empty in `config.yaml`, CheckM2 will auto-download the database on first use to `~/.cache/checkm2/`. This download repeats for every new user unless a shared path is configured.

To pre-download to a shared location:

```bash
mamba create -n db_setup -c conda-forge -c bioconda checkm2 -y
conda activate db_setup

checkm2 database --download --path /your/shared/dbs/checkm2/

conda deactivate
ls /your/shared/dbs/checkm2/   # note the exact .dmnd filename
```

Then set `params.checkm2.db_path` to the full path of the downloaded `.dmnd` file in `config.yaml`.

---

## GUNC

GUNC detects chimeric and contaminated MAG bins.

If `params.gunc.db_path` is empty in `config.yaml`, GUNC will auto-download the database on first use to `~/.gunc/`. This download repeats for every new user unless a shared path is configured.

To pre-download to a shared location:

```bash
mamba create -n db_setup -c conda-forge -c bioconda gunc -y
conda activate db_setup

gunc download_db /your/shared/dbs/gunc/

conda deactivate
ls /your/shared/dbs/gunc/   # note the exact .dmnd filename
```

Then set `params.gunc.db_path` to the full path of the downloaded `.dmnd` file in `config.yaml`.

---

## GTDB-Tk

GTDB-Tk classifies MAGs against the Genome Taxonomy Database. It requires a reference data package (~100 GB extracted for r220+) that must be downloaded before the taxonomy step will run.

Download the current release from the [GTDB data server](https://data.gtdb.ecogenomics.org/releases/). Browse to the latest release directory and download the `gtdbtk_data_r*.tar.gz` package:

```bash
mkdir -p /your/dbs/gtdbtk
cd /your/dbs/gtdbtk

# Check data.gtdb.ecogenomics.org/releases/ for the current package URL
wget -c https://data.gtdb.ecogenomics.org/releases/release232/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz

tar -xzf gtdbtk_r232_data.tar.gz

# Confirm what directory was extracted
ls /your/dbs/gtdbtk/
```

The `-c` flag on `wget` allows resuming interrupted downloads. After extraction, set `params.gtdbtk.db_path` to the directory containing the unpacked reference data in `config.yaml`.

See the [GTDB-Tk documentation](https://ecogenomics.github.io/GTDBTk/installing/index.html) and [GTDB releases page](https://gtdb.ecogenomics.org/) for the current database version.

> **Note:** `gtdbtk download_db` is not a valid subcommand — the database must be downloaded manually as shown above.

GTDB-Tk's archaeal summary (`gtdbtk.ar53.summary.tsv`) is only written when a
sample has archaeal MAGs; `make_mag_summary` handles its absence. GTDB-Tk runs
here as three stages (`identify` → `align` → `classify`); see the GTDB-Tk
stages section of [Output](output.md).

---

## Kraken2 and Bracken

Kraken2 (with Bracken) does k-mer read classification, enabled by the
`profilers` list in `config.yaml`. Build or download a Kraken2 database (the
Standard database is ~200 GB) and set both keys to the **same build**:

```yaml
params:
  kraken2:
    db: /your/dbs/kraken2/Kraken2_db_Standard
    bracken-db: /your/dbs/kraken2/Kraken2_db_Standard
```

Bracken's `.kmer_distrib` files are specific to the Kraken2 build they were
generated from. Pointing `bracken-db` at a different build silently produces
wrong abundance estimates. If the `.kmer_distrib` files for your read lengths
are absent, build them once with `bracken-build` (see the
[Bracken docs](https://github.com/jenniferlu717/Bracken)); the demon Standard
database already has them for 100/150/200/250/300 bp.

---

## MetaPhlAn 4

MetaPhlAn 4 uses a marker-gene database for taxonomic profiling. Download it with:

```bash
conda activate snakemake   # or any env with metaphlan installed
metaphlan --install --bowtie2db /your/dbs/metaphlan --index mpa_vJan25_CHOCOPhlAnSGB_202503
```

Then set in `config.yaml`:

```yaml
params:
  metaphlan:
    db_path: /your/dbs/metaphlan
    db_name: mpa_vJan25_CHOCOPhlAnSGB_202503
```

Check the [MetaPhlAn wiki](https://huttenhower.sph.harvard.edu/metaphlan/) for the current database name.

---

## geNomad

geNomad classifies viral and plasmid contigs for the optional viral track.

```bash
mamba create -n db_setup -c conda-forge -c bioconda genomad -y
conda activate db_setup

genomad download-database /your/dbs/genomad/

conda deactivate
```

This creates `/your/dbs/genomad/genomad_db/`. Set `params.genomad.db_path` to
that `genomad_db` directory.

---

## CheckV

CheckV scores viral completeness and contamination for the geNomad viral
contigs.

```bash
mamba create -n db_setup -c conda-forge -c bioconda checkv -y
conda activate db_setup

checkv download_database /your/dbs/checkv/

conda deactivate
ls /your/dbs/checkv/          # note the checkv-db-v1.x directory
```

Set `params.checkv.db_path` to the `checkv-db-v1.x` directory.

---

## Host genome

The host genome FASTA is used by bowtie2 to remove host reads. Any genome FASTA will work; the bowtie2 index is built automatically on first run.

```yaml
host_filter:
  genome: /path/to/host_genome.fna
  db_dir: /path/to/index_directory/
```

The index is built once in `db_dir` using the FASTA filename stem as the index prefix. If the index already exists, the build step is skipped.

Common genomes:
- Human: [GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)
- Mouse: [GRCm39](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/)

> **WFUSM users:** All databases are pre-downloaded on DEMON. See [Running on DEMON](demon.md) for paths and configuration.

---

← [Documentation home](index.md)
