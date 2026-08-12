# Configuration

All configuration lives in `resources/config/`. The `config.yaml` at the repository root is a symlink to `resources/config/config.yaml`.

---

## `config.yaml`

Key settings to review before running:

### Trimmer

```yaml
trimmer: fastp   # fastp (default) or cutadapt
```

Both options run FastQC before and after trimming and feed into MultiQC.

Cutadapt requires adapter sequences in `params.cutadapt`. fastp is configured
in `params.fastp.extra`, which by default passes:

```
--detect_adapter_for_pe --adapter_fasta resources/adapters/illumina_adapters.fa --trim_poly_g --trim_poly_x
```

These are set explicitly rather than relying on fastp's defaults, because
the defaults are conditional and quietly do less than expected:

- **Adapter detection is weaker for paired-end input.** fastp only
  auto-detects adapter *sequences* for single-end reads. For paired-end it
  infers adapters from per-read overlap analysis, which leaves residual 3'
  adapter behind when the overlap is short or low quality.
  `--detect_adapter_for_pe` turns sequence detection on, and
  `--adapter_fasta` additionally trims the known TruSeq/Nextera sequences.
- **polyG trimming depends on the read header.** fastp enables it only when
  it recognises a NextSeq/NovaSeq instrument ID in the read name, and for
  archived data that is unreliable. Deflines differ across SRA datasets:

  ```
  @SRR36840144.1 A00814:609:HMHW7DSX3:1:1101:10004:4726 length=151
  @SRR36840144.1 1 length=151
  @A00814:715:HVFYYDRX2:1:1101:18114:1016 1:N:0:...          (native)
  ```

  Some runs carry the original Illumina name and some do not, and even when
  present it is not in the leading field where native output puts it. So
  whether 2-colour polyG tails get trimmed can differ between datasets that
  are meant to be compared directly. `--trim_poly_g` removes
  the dependence on header parsing, and `--trim_poly_x` additionally catches
  polyA read-through, which the adapter panel cannot cover.

### Adapter panel

`resources/adapters/illumina_adapters.fa` is a deliberately **comprehensive
guard set**, not a minimal one, so the same file can be reused across DNA and
RNA library types without editing. It covers:

| Group | Constructs |
|---|---|
| TruSeq | Read 1 / Read 2 adapters, indexed-adapter pre-index segment, post-index constant region, PE bottom adapter |
| Flowcell primers | P5, P7 (full and short forms) |
| Sequencing primers | Read 1, Read 2 |
| Nextera / Tn5 | mosaic end (both orientations), Read 1 / Read 2 adapters, i5 and i7 transposome sequences |
| Tagmentation / RNA prep | Illumina DNA PCR-Free, Illumina Stranded mRNA and Total RNA Prep |
| Small RNA | TruSeq Small RNA 5' and 3' adapters |

Sequences are taken from Illumina's
[adapter-trimming reference](https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314)
and [Illumina Adapter Sequences (doc #1000000002694)](https://support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html).

Two constructs are worth calling out because they are commonly omitted and
were the source of real misassignments: the **P5/Read 1 sequencing primer**
(`ACACTCTTTCCCTACACGACGCTCTTCCGATCT`) and the **post-index constant region**
(`ATCTCGTATGCCGTCTTCTGCTTG`, the reverse complement of the P7 primer). Reads
consisting entirely of these have been mistaken for viral genomes, because
some public assemblies contain untrimmed adapter at the matching position.

Because the panel favours coverage over minimality, a few short entries
(16-20 bp) carry a slightly higher chance of spurious trimming than the
33-61 bp constructs. That trade is intentional. Set
`params.fastp.extra: ""` to restore fastp's stock behaviour, or point
`--adapter_fasta` at a trimmed-down file for a specific library type.

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

The index name is derived from the FASTA filename stem, so `genome: /path/human.fna` expects `/path/human.*.bt2` in `db_dir`.

A **complete** existing index is detected and reused — useful for shared references, where rebuilding a human genome costs hours. The index is built only when all six parts (`.1`, `.2`, `.3`, `.4`, `.rev.1`, `.rev.2`, in either `.bt2` or `.bt2l` form) are not already present, so a partial index left by an interrupted build is rebuilt rather than silently accepted. Completion is tracked with a `<stem>.bowtie2_build.done` sentinel, since bowtie2 chooses between the small and large index layouts based on reference size and the output filenames cannot be declared up front.

### Prototype selection

```yaml
params:
  prototypes:
    n: 10                  # representative samples selected for binning
    min_seqs: 50           # depth floor; use 10000000 for real data
    max_seqs: 200000000    # depth ceiling
    exclude: []            # samples barred from the prototype pool
```

Both depth bounds are compared against fastp's `total_reads`, which counts
**both mates** — a library described as "50M reads" is 100M by this measure,
so the ceiling bites at half the depth you might expect.

`max_seqs` excludes the *deepest* samples, and those are exactly the ones
that make the best differential-coverage basis for binning. The previous
`100000000` default silently dropped the top 4% of a 341-sample gut cohort
(deepest run 138.5M reads) without warning. Sketching is cheap, so the
ceiling is now high enough to retain everything and serves only as a guard
against a pathological input. Check your own depth distribution before
lowering it.

`n` determines how many prototype samples are selected by `prototype_selection` and used by `generate_binning_config` to populate `binning.txt`. Setting `n` higher produces better binning coverage but more assembly/mapping jobs.

`exclude` lists samples by name that may never serve as a prototype:

```yaml
params:
  prototypes:
    exclude:
      - Control_Sample
```

A run control is the usual case. It should still be trimmed, screened and
profiled, because you want to know what is in it, but it is not part of the
biology and has no business forming a differential-coverage basis for
binning. The depth bounds cannot express this, since a control can be
sequenced as deeply as any real sample.

A name that matches no sample in the metadata table stops the run at load
with an error naming it. Silently ignoring a typo would leave the very
sample you meant to exclude still acting as a coverage basis, which is the
failure this option exists to prevent.

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
  gtdbtk: 128000     # first attempt
  gtdbtk_max: 320000 # ceiling; retries request gtdbtk * attempt, capped
```

Assembly rules (`megahit`, `metaspades`) auto-scale their memory request based on input size (`max(16000, input_size_mb × 10)`) up to the configured ceiling. All other rules use their configured value directly.

**GTDB-Tk is the one to watch, and it escalates.** Its memory is dominated by
pplacer during `classify_wf`, and the requirement scales with the number of
bins: on a 341-sample gut cohort, 324 samples finished inside 128 GB while 17
were OOM-killed. The failure surfaces only at the very end of the pipeline,
since MAG QC is the last stage.

Sizing every request for the worst case is the wrong fix. `mem_mb` is a
scheduler *reservation*, so a flat 300 GB request confines every sample to
whichever nodes are that large and serialises the stage. Instead `gtdbtk` is
the **first attempt** and `gtdbtk_max` the **ceiling**; the rule requests
`gtdbtk × attempt` capped at `gtdbtk_max`, with retries from
`retries.run_gtdbtk`. Most samples schedule anywhere on the first try and
only the heavy ones wait for a large node:

| attempt | request |
|---|---|
| 1 | 128 GB |
| 2 | 256 GB |
| 3 | 320 GB (capped) |

Check what your nodes can actually supply with `sinfo -o "%n %m"`. If none
reach the ceiling, pass `--scratch_dir` to gtdbtk instead: pplacer then mmaps
the reference to disk, needing far less RAM but running slower.

Runtimes come from the cluster profile rather than `config.yaml`. Note that
**`fastp` is the rule most exposed to slow shared storage**, since it reads
every raw FASTQ: on the `demon` profile it is given 360 minutes rather than
the 120-minute default, after a 2.4 GB library was cancelled four minutes
inside the limit while the filesystem was saturated. The compute is minutes;
the variance is I/O.

---

## `metadata.txt`

One tab-separated table describing every sequencing run. This replaces the
former `samples.txt` + `units.txt` pair.

**Required columns:** `Sample`, `Sequencing_Run`, `R1_fp`, `R2_fp`

```
Sample   Sequencing_Run   R1_fp                       R2_fp                       Treatment_Group   Timepoint
John     Run_1            /path/John_R1.fastq.gz      /path/John_R2.fastq.gz      Treatment         1963
Paul     Run_1            /path/Paul_R1.fastq.gz      /path/Paul_R2.fastq.gz      Treatment         1963
Paul     Run_2            /path/Paul_L2_R1.fastq.gz   /path/Paul_L2_R2.fastq.gz   Treatment         1963
George   Run_1            /path/George_R1.fastq.gz    /path/George_R2.fastq.gz    Control           1963
```

One row per sequencing run: a sample sequenced on two lanes, or resequenced,
gets one row per run with the same `Sample` value. The sample list is taken
from the unique values of `Sample`, so no separate sample sheet is needed.

The `_fp` suffix marks the two columns that hold **f**ile **p**aths. Paths may
be absolute or relative to the working directory.

**Any further columns are yours.** Study covariates -- treatment, timepoint,
subject, batch -- are carried through untouched and can be read inside a rule:

```python
metadata_table.loc[(sample, seqrun), "Treatment_Group"]
```

Multiple runs for the same sample are concatenated by the `merge_seqruns` rule
before assembly. A sample with a single run skips that rule entirely: its
trimmed FASTQ is used directly, so nothing is copied or linked.

The table is validated on load. Missing required columns, empty required
fields, duplicate `(Sample, Sequencing_Run)` pairs, and a file that is
space- rather than tab-separated each produce a specific error naming the
offending rows, rather than failing later inside the workflow.

### Single-end runs

A run is single-end when `R2_fp` holds the literal `NA`:

```
Sample   Sequencing_Run   R1_fp                      R2_fp                      Treatment_Group
John     Run_1            /path/John_R1.fastq.gz     /path/John_R2.fastq.gz     Treatment
George   Run_1            /path/George_R1.fastq.gz   NA                         Control
```

Nothing is inferred. The layout is not detected from the filesystem, from
read headers, or from a missing column, and a blank `R2_fp` is an error
rather than an implicit single-end declaration. That is deliberate: a
truncated or half-written sample sheet should fail loudly instead of
quietly assembling one mate of a paired library. The marker is exact and
case-sensitive, so `na`, `n/a` and `None` are all still errors.

Layout may vary between samples in one cohort, but not between the runs of
a single sample. `merge_seqruns` concatenates per read identifier, so a
sample mixing layouts would produce a reverse file covering only some of
its runs; this is rejected on load.

Single-end samples differ from paired-end ones in three ways:

| | Paired-end | Single-end |
|---|---|---|
| Read identifier | `R1`, `R2` | `SE` |
| Assemblers | megahit, metaspades | **megahit only** |
| Trimmer | fastp or cutadapt | **fastp only** |

metaSPAdes rejects a single-end-only library outright, so single-end
samples are dropped from the metaSPAdes target list and assembled by
MEGAHIT. With `assemblers: [metaspades, megahit]` a mixed cohort is
handled correctly: paired-end samples are assembled by both, single-end
samples by MEGAHIT alone, and `generate_binning_config` points each sample
at an assembly that exists.

With `assemblers: [metaspades]` alone, a single-end sample has no
assembler at all, and the run stops immediately with an error naming it.
The check is deliberately at load time. Without it such a sample is
trimmed, host filtered and profiled and only then dropped, so stage 1
reports success having produced no contigs for it and the first complaint
comes much later from `generate_binning_config`.

The cutadapt parameters are paired-end specific (`-U` trims the reverse
read), so `trimmer: cutadapt` is refused when the cohort contains
single-end samples.

fastp runs with `params.fastp.extra` minus `--detect_adapter_for_pe`,
which applies only to paired input; fastp detects adapter sequence for
single-end reads by default. Set `params.fastp.extra_se` to override that
derived value.

Single-end outputs never collide with paired-end ones: trimmed reads and
their fastp reports go to `output/qc/fastp/se/`, host-filtered reads to
`output/qc/host_filter/nonhost/{sample}.SE.fastq.gz`, and the host BAM to
`output/qc/host_filter/host_se/`. Paired-end paths are unchanged, so
adding single-end samples to an existing project does not invalidate work
already done.

`resources/test/metadata_single_end.txt` and
`resources/test/metadata_mixed_layout.txt` are runnable examples against
the bundled test reads. The loader behaviour is covered by
`python3 resources/test/test_metadata_layout.py`.

### Migrating from samples.txt + units.txt

`samples.txt` is no longer used -- its only role was listing sample names,
which `metadata.txt` already provides. Rename the `units.txt` header:

```bash
sed '1s/.*/Sample\tSequencing_Run\tR1_fp\tR2_fp/' units.txt > metadata.txt
```

then replace `samples:` and `units:` in your config with a single
`metadata:` key. A config still using the old keys fails immediately with
these instructions rather than a confusing error.

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
- A sample with neither label is present in `metadata.txt` but skipped by the binning pipeline.

**This file is normally generated automatically** by the `generate_binning_config` rule (see [Running the pipeline](running.md)).

---

← [Back to README](../README.md)
