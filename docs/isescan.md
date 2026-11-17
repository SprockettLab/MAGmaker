# Insertion sequences (ISEScan)

[ISEScan](https://github.com/xiezhq/ISEScan) (Xie and Tang 2017,
*Bioinformatics* 33:3340) finds insertion sequence (IS) elements, complete and
partial, in prokaryotic sequence. MAGmaker runs it in Stage 3, off by default.

```yaml
isescan: false           # true to run it on every run
params:
  isescan:
    extra: ""            # e.g. "--removeShortIS" for complete IS elements only
threads:
  isescan: 32
```

For one run only: `./run_magmaker.sh --isescan ...`. To run it on its own in
a run folder whose MAG track has finished:

```bash
snakemake --snakefile Snakefile-bin --profile resources/profiles/demon isescan_all \
  --config binning=output/config/auto_binning.txt isescan=True
```

ISEScan is Linux only.

---

## What it runs on, and why

**Every contig of each binning row's own assembly, not the MAGs.** IS
elements are repeats present in several copies, so assemblers often break at
them and binners often leave the resulting contigs out of bins (their
coverage and composition do not match the genome's). Running on the bins
would undercount IS elements, worst in exactly the IS-rich strains.

**Each call is then labelled with the MAG its contig was binned into**, if
any, with that MAG's taxonomy and MIMAG quality from `mag_summary.tsv`. Calls
on unbinned contigs stay in the table with `NA`.

**A per-species summary is built downstream from this table**, by mapping
MAGs to their species-level genome bins, without running ISEScan again. A
dereplicated genome set alone cannot show IS presence or absence differing
between strains or samples, which is usually the question.

**Partial elements.** Short-read assemblies often break at IS copies, so many
calls there are partial elements at contig ends. ISEScan's `type` column says
which (`c` complete, `p` partial); report the two separately. `extra:
"--removeShortIS"` keeps only complete elements.

---

## Output

```
output/isescan/
├── {assembler}/{row}/
│   ├── is.tsv           ISEScan's own table for that assembly
│   ├── run/             ISEScan's other outputs (.sum, .gff, .is.fna, .orf.faa, ...)
│   ├── params.txt       ISEScan version, options, input, contig and IS counts
│   └── NO_IS_FOUND      only when its HMM search found no transposase at all
└── isescan_summary.tsv  every row's calls, labelled with sample and MAG
```

`isescan_summary.tsv` has one row per IS element copy:

| Column | Meaning |
|---|---|
| `Sample_ID`, `Assembler`, `Binning_Row` | where the assembly came from |
| `Contig` | first word of ISEScan's `seqID` (ISEScan keeps the whole FASTA header, e.g. MEGAHIT's `flag=… multi=… len=…`) |
| ISEScan's columns | `seqID`, `family`, `cluster`, `isBegin`, `isEnd`, `isLen`, `ncopy4is`, TIR and ORF columns, `E-value`, `E-value4copy`, `type`, `ov`, `tir`, as ISEScan writes them |
| `MAG_ID`, `MAG_New_Name`, `MAG_GTDB_Classification`, `MAG_MIMAG_Quality` | the MAG holding the contig; `NA` if none; `;`-joined if MAGs from more than one mapper hold it |

Contig names repeat between assemblies (MEGAHIT numbers every assembly from
`k141_0`), so a MAG's contigs are matched only within the binning row its
FASTA came from.

The summary's log (`output/logs/isescan/make_isescan_summary.log`) gives the
IS copies per row and how many sit on binned contigs, and names any MAG FASTA
it could not read.

---

## Runtime

Measured on WF18 on demon, one assembly per job, 32 threads, ISEScan 1.7.3:

| Assembly | Size | Contigs | Time | Peak memory |
|---|---|---|---|---|
| flye (long reads) | 82.1 Mbp | 1,526 | 19.6 min | 607 MB |
| MEGAHIT, median of 88 by size | 214.8 Mbp | 241,571 | 82.8 min | 753 MB |
| myloasm, largest of 96 | 2,261 Mbp | 60,457 | 6.46 h | 5.4 GB |

About 10–23 s per Mbp, slowest on short-read assemblies with many small
contigs. The rule reserves 16 GB and 24 h on the first attempt, doubling
memory and adding runtime on each retry (`mem_mb.run_isescan`,
`runtime.run_isescan`).

---

## Traps

- **A retry starts from an empty folder.** ISEScan skips any HMM search whose
  result file already exists, so resuming after a killed attempt would build
  calls from a half-written search.
- **ISEScan writes a file list into its working directory**, named after the
  input's folder and file. The rule runs it from inside its own output folder
  so parallel jobs cannot collide.
- **"Found nothing" and "wrote nothing" are kept apart.** When the HMM search
  finds no transposase at all, ISEScan writes no table; the rule writes an
  empty `is.tsv` and `NO_IS_FOUND`. Any other missing table fails the job.
