# Configuration

All configuration lives in `resources/config/`. The `config.yaml` at the repository root is a symlink to `resources/config/config.yaml`.

---

## `config.yaml`

Key settings to review before running:

### Trimmer

```yaml
trimmer: fastp   # fastp (default) or cutadapt
```

fastp auto-detects adapters. Cutadapt requires adapter sequences in `params.cutadapt`. Both options run FastQC before and after trimming and feed into MultiQC.

### Assembler(s)

```yaml
assemblers:
  - megahit       # default
# - metaspades    # uncomment to also run metaSPAdes
```

Both assemblers can run in parallel. MEGAHIT is the default — faster and more memory-efficient for most metagenomes. metaSPAdes may produce better assemblies for lower-complexity samples but requires substantially more RAM.

### Host filter

```yaml
host_filter:
  genome: /path/to/host_genome.fna
  db_dir: /path/to/bt2_index_dir/
```

The bowtie2 index is built automatically from the FASTA on first run if it doesn't exist in `db_dir`. The index name is derived from the FASTA filename stem.

### Prototype selection

```yaml
params:
  prototypes:
    n: 10            # number of representative samples to select for binning
    min_seqs: 50     # minimum reads to include a sample in sourmash sketch
```

`n` determines how many prototype samples are selected by `prototype_selection` and used by `generate_binning_config` to populate `binning.txt`. Setting `n` higher produces better binning coverage but more assembly/mapping jobs.

### Mappers and binners (binning pipeline)

```yaml
mappers:
  - minimap2      # default
# - bowtie2       # uncomment to also run bowtie2

binners:
  - concoct
  - metabat2
  - maxbin2
```

All enabled binners run independently and their results are combined by DAS_Tool.

### Taxonomy and profiling

```yaml
params:
  metaphlan:
    db_path: /path/to/metaphlan_db/
    db_name: mpa_vJan25_CHOCOPhlAnSGB_202503
  gtdbtk:
    db_path: /path/to/gtdbtk/release232/
  checkm2:
    db_path: /path/to/checkm2/uniref100.KO.1.dmnd   # leave empty to auto-download
  gunc:
    db_path: /path/to/gunc_db_progenomes2.1.dmnd     # leave empty to auto-download
```

See [Database setup](databases.md) for download instructions.

### Threads and memory

```yaml
threads:
  megahit: 16
  checkm2: 16
  gtdbtk: 16
  metaphlan: 8
  # ... one entry per rule

mem_mb:
  megahit: 256000    # ceiling; actual request auto-scales with input size
  spades: 256000
  checkm2: 32000
  gtdbtk: 128000
```

Assembly rules (`megahit`, `metaspades`) auto-scale their memory request based on input size (`max(16000, input_size_mb × 10)`) up to the configured ceiling. All other rules use their configured value directly.

---

## `samples.txt`

Tab-separated. One row per sample. The first column must be named `Sample`.

```
Sample    Subject    Timepoint
John      Beatles    1963
Paul      Beatles    1963
```

Additional columns are carried through but not used by the pipeline.

---

## `units.txt`

Tab-separated. One row per sequencing unit (e.g., per lane or per run). Columns: `Sample`, `Unit`, `R1`, `R2`.

```
Sample    Unit      R1                         R2
John      Run_1     /path/to/John_R1.fastq.gz  /path/to/John_R2.fastq.gz
Paul      Run_1     /path/to/Paul_R1.fastq.gz  /path/to/Paul_R2.fastq.gz
Paul      Run_2     /path/to/Paul_lane2_R1.gz  /path/to/Paul_lane2_R2.gz
```

Multiple units for the same sample are automatically concatenated by the `merge_units` rule before assembly. If a sample has only one unit, a symlink is created instead of copying (no I/O overhead).

Paths in `R1` and `R2` can be absolute or relative to the working directory.

---

## `binning.txt` (binning pipeline only)

Tab-separated. Defines which reads are mapped to which assemblies for binning. Columns: `Sample`, `Contigs`, `Read_Groups`, `Contig_Groups`.

```
Sample    Contigs                                        Read_Groups    Contig_Groups
John      output/assemble/megahit/John.contigs.fasta    A              A
Paul      output/assemble/megahit/Paul.contigs.fasta    A              A
George    output/assemble/megahit/George.contigs.fasta                 A
Ringo                                                   A
```

Samples that share a group label in both `Read_Groups` and `Contig_Groups` are paired: reads from all read-samples in a group are mapped to all contig-samples in that group.

- A sample with a `Contigs` path and a `Contig_Groups` label contributes an assembly to that group.
- A sample with a `Read_Groups` label contributes reads to that group.
- A sample can belong to both (contributing both reads and an assembly).
- A sample with neither label is present in `samples.txt` but skipped by the binning pipeline.

**This file is normally generated automatically** by the `generate_binning_config` rule (see [Running the pipeline](running.md)).

---

← [Back to README](../README.md)
