#!/usr/bin/env bash
#SBATCH --job-name=bracken_build
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/kraken2/bracken_build_%j.log
#SBATCH --error=/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/kraken2/bracken_build_%j.err

# Build Bracken k-mer distribution files for Kraken2_db_Standard.
# Must be run once per read length. Re-run with a different -l value if your
# project uses a different read length.
#
# Usage:
#   sbatch build_bracken_db.sh           # default: 150 bp reads
#   sbatch build_bracken_db.sh 100       # 100 bp reads
#
# After completion, confirm that *.kmer_distrib files exist in the DB directory.
# These are required by the bracken rule in profile.smk.

set -euo pipefail

READ_LENGTH=${1:-150}
DB=/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/kraken2/Kraken2_db_Standard
THREADS=${SLURM_CPUS_PER_TASK:-16}

echo "Building Bracken database"
echo "  DB:          $DB"
echo "  Read length: $READ_LENGTH bp"
echo "  Threads:     $THREADS"
echo "  Started:     $(date)"

source /isilon/datalake/sprockett_lab/original/WF00SprockettLab/envs/envs/snakemake/etc/profile.d/conda.sh
conda activate /isilon/datalake/sprockett_lab/original/WF00SprockettLab/envs/envs/snakemake

bracken-build \
    -d "$DB" \
    -t "$THREADS" \
    -k 35 \
    -l "$READ_LENGTH"

echo "Done: $(date)"
echo "kmer_distrib files:"
ls "$DB"/*.kmer_distrib
