# Pipeline overview

MAGmaker runs as **two Snakemake workflows**, chained by `run_magmaker.sh`:

| Snakefile | Stage | Does |
|---|---|---|
| `Snakefile` | 1 | QC → assembly → read-level profiling → prototype selection |
| `Snakefile-bin` | 2 | mapping → binning → bin consolidation → MAG QC (+ optional viral track) |

They share the sample sheet (`resources/snakefiles/common.smk`) so the two
stages cannot disagree about it. The bridge between them is
`generate_binning_config`, which turns prototype selection into a binning plan.

---

## Stage 1 — `Snakefile`

### 1. Input

One tab-separated `metadata.txt`: `Sample`, `Sequencing_Run`, `R1_fp`,
`R2_fp`, and any covariate columns you add (carried through untouched). A run
is single-end when `R2_fp` is the literal `NA` — nothing is inferred. Plus a
host genome FASTA for read removal.

### 2. Quality filtering

FastQC (pre-trim) → **fastp** (default) or Cutadapt → FastQC (post-trim) →
`merge_seqruns` for multi-run samples → **host removal** (Bowtie2; non-host
reads kept, host reads discarded as a BAM) → FastQC (post-host) → **MultiQC**.

fastp's adapter and homopolymer trimming are made explicit
(`--detect_adapter_for_pe --adapter_fasta … --trim_poly_g --trim_poly_x`)
rather than left to its header-dependent defaults, which do less than expected
on archived/SRA data. See [Configuration](configuration.md).

The host-filtered non-host reads are the input to everything downstream:
assembly, profiling, prototype selection, mapping, and CMSeq.

### 3. Assembly

Every sample is assembled by **MEGAHIT** (default) and/or **metaSPAdes**.
Single-end samples are MEGAHIT-only (metaSPAdes rejects a single-end library).
QUAST → MultiQC for assembly stats.

### 4. Read-level profiling (optional)

- **MetaPhlAn 4** (`biobakery` list) → per-sample profiles → merged abundance
  table.
- **Kraken 2 + Bracken** (`profilers` list) → per-sample reports → Krona
  plots.

Both are toggled by their config list or by `--skip-metaphlan` /
`--skip-kraken`. Kraken2 reads a ~200 GB index per sample and is the heaviest
I/O in the pipeline — an analysis taking its taxonomy from GTDB-Tk on the MAGs
needs neither tool.

### 5. Prototype selection

sourmash sketches each sample's non-host reads → `sourmash compare` builds a
pairwise MinHash distance matrix (and a heatmap) → a greedy max-distance
heuristic picks a spread of representative samples, filtered by read depth and
honouring `params.prototypes.exclude`.

`generate_binning_config` then writes `output/config/auto_binning.txt`:
**every sample's assembly is binned**, and the `n` prototype samples
(`params.prototypes.n`) supply the reads that form the differential-coverage
basis for binning all of them.

---

## Stage 2 — `Snakefile-bin`

### 6. Mapping — differential coverage

Each assembly is indexed (minimap2 default, or Bowtie2) and the prototype
samples' non-host reads are mapped to it, producing one sorted, indexed BAM
per (prototype × assembly). This coverage profile across samples is what lets
the binners separate genomes that composition alone cannot.

### 7. Binning

Four binners run independently on each assembly + its BAMs:

| Binner | Basis |
|---|---|
| MetaBAT2 | tetranucleotide composition + coverage |
| MaxBin2 | composition + coverage + marker genes (EM) |
| CONCOCT | k-mer composition + coverage clustering |
| SemiBin2 | self-supervised neural embedding |

SemiBin2 is a different *mechanism* from the other three, so it fails
differently — the value of a fourth binner is a candidate no other produced.
A binner may decline a sample whose assembly is too sparse; that costs one
binner's contribution, not the sample.

### 8. Bin consolidation

The four per-binner bin sets are reconciled into one non-redundant MAG set per
sample by **Binette** (default) or **DAS_Tool** (`params.consolidation.tool`,
or `--binette` / `--das-tool`).

- Binette scores candidates with CheckM2 — the same estimator MAG QC uses —
  and builds extra candidates from the intersection / difference / union of
  overlapping bins.
- DAS_Tool scores with 51 single-copy genes and a `score_threshold`.

Selected FASTAs land in `DAS_Tool_Fastas/` or `Binette_Fastas/`; everything
after this reads whichever the run chose. See [Configuration](configuration.md).

### 9. MAG QC and taxonomy

Each MAG set is assessed by:

- **CheckM2** — completeness / contamination
- **GUNC** — chimerism against a reference database
- **GTDB-Tk** — run as three stages: `identify` (marker genes) → `align`
  (concatenated MSA, and `MSA_Percent`) → `classify` (pplacer placement,
  taxonomy). `--skip-gtdbtk-classify` stops after `align`; see
  [Output](output.md).
- **CMSeq** — strain heterogeneity: the sample's own reads mapped back to its
  own MAGs, counting polymorphic positions. Catches a MAG that is a blended
  consensus of co-resident strains, which CheckM2 and GUNC cannot see. Skip
  with `--skip-cmseq`.

`make_mag_summary` pools everything into `output/mag_qc/mag_summary.tsv` — one
row per MAG, globally numbered (`MAG_0001`…), sorted by taxonomy then MIMAG
quality tier, with `NONE` rows recording samples that produced no MAG.
`rename_mags` copies the FASTAs into `renamed_mags/` under the `New_Name`
column, which you can edit and re-run.

### 10. Viral / plasmid track (optional)

Branches off the assemblies: **geNomad** classifies viral and plasmid contigs
(no mapping needed), **CheckV** scores the viral ones, and
`make_virus_summary` writes `virus_summary.tsv`. Requested with the
`virus_all` target. See [Viral / plasmid track](viral.md).

---

## Running it

Every step runs in its own conda environment. On an HPC cluster each rule is a
separate SLURM job, and most rules grow their memory / runtime request on each
retry. See [Running the pipeline](running.md).

---

← [Documentation home](index.md)
