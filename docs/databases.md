# Database setup

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

## GTDB-tk

GTDB-tk classifies MAGs against the Genome Taxonomy Database. It requires a reference data package (~85 GB) that must be downloaded before the taxonomy step will run.

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

See the [GTDB-tk documentation](https://ecogenomics.github.io/GTDBTk/installing/index.html) and [GTDB releases page](https://gtdb.ecogenomics.org/) for the current database version.

> **Note:** `gtdbtk download_db` is not a valid subcommand — the database must be downloaded manually as shown above.

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

geNomad identifies and classifies viral and plasmid contigs in the [viral track](viral.md). It uses a marker and classifier database that must be downloaded before the track will run.

```bash
mamba create -n db_setup -c conda-forge -c bioconda genomad -y
conda activate db_setup

genomad download-database /your/shared/dbs/genomad/   # creates genomad_db/

conda deactivate
ls /your/shared/dbs/genomad/genomad_db/   # confirm the database directory
```

Then set `params.genomad.db_path` to the `genomad_db` directory in `config.yaml`. `params.genomad.extra` passes extra flags to `genomad end-to-end` (e.g. `--conservative` for higher-precision viral calls).

---

## CheckV

CheckV scores the completeness and contamination of the viral contigs geNomad finds. It uses a reference database that must be downloaded before the track will run.

```bash
mamba create -n db_setup -c conda-forge -c bioconda checkv -y
conda activate db_setup

checkv download_database /your/shared/dbs/checkv/   # creates checkv-db-v1.5/

conda deactivate
ls /your/shared/dbs/checkv/   # note the exact checkv-db-v* directory name
```

Then set `params.checkv.db_path` to the downloaded `checkv-db-v*` directory in `config.yaml`.

> The viral track is on by default. If the geNomad/CheckV databases are not configured, run with `--skip-virus` (see [Viral track](viral.md)).

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

← [Back to README](../README.md)
