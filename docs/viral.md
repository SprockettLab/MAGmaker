# Viral / plasmid track

A per-sample viral and plasmid track that runs parallel to the MAG track and
does not touch it. `run_magmaker.sh` runs it in Stage 3 as a second target
next to `rename_mags`; `--skip-virus` leaves it out for one run. To run it on
its own (for example in a run folder whose MAG track already finished):

```bash
snakemake --snakefile Snakefile-bin --profile resources/profiles/demon virus_all \
  --config binning=output/config/auto_binning.txt
```

---

## What it does

```
per-sample contigs ─┬─ geNomad ──────────┬─ viral contigs ─ CheckV ─┐
                    ├─ VIBRANT ──────────┤                          │
                    ├─ Cenote-Taker3 ────┤   (CheckV per predictor)  ├─ virus_summary.tsv
                    └─ VirSorter2 ───────┘                          │  virus_name_mapping.tsv
                       geNomad only ────── plasmid contigs ─────────┘

geNomad's viral contigs + MetaBAT2 coverage table ─ vRhyme ─ vMAGs
                                                             vrhyme_summary.tsv
                                                             vrhyme_membership.tsv
```

**Predictors.** Four, each with its own switch. geNomad is on by default; the
other three are off until a benchmark decides, because they disagree with each
other far more than they agree (see *Agreement*, below).

**CheckV runs after each predictor, on that predictor's own calls.** It
corrects predicted boundaries and the tools differ in exactly that, so running
it once on a pooled set would make quality identical by construction and
discard what is being compared.

**vRhyme is a binner, not a predictor.** It joins geNomad's viral contigs into
viral bins (vMAGs), so a phage that assembled in several pieces is reported as
one genome instead of several partial ones. It computes nothing new: it reads
geNomad's viral FASTA and the MetaBAT2 coverage table the MAG track already
wrote. It is on by default and needs `genomad: true` and `metabat2` among the
`binners`; without either, the run says once that vRhyme will not run, and why.

**Plasmids are geNomad's.** It is the only one of the four that calls them,
and they come out of the same single run, so geNomad is one switch for both:
with it on you get viruses and plasmids, with it off you get neither.

---

## Configuration

```yaml
virus:
  genomad: true          # viruses AND plasmids; also what vRhyme bins
  vrhyme: true           # viral binning; needs genomad + metabat2
  vibrant: false
  cenotetaker3: false
  virsorter2: false

params:
  genomad:
    db_path: /path/to/genomad_db        # genomad download-database <dir>
    extra: ''                           # e.g. --conservative
    score_calibration: true             # FDR column + geNomad's --max-fdr filter
  checkv:
    db_path: /path/to/checkv-db-v1.5    # checkv download_database <dir>
  vibrant:
    env: /path/to/envs/vibrant-1.2.1
    db_dir: /path/to/dbs/vibrant/databases
    files_dir: /path/to/dbs/vibrant/files
  cenotetaker3:
    env: /path/to/envs/cenotetaker3-3.4.4
    db_dir: /path/to/dbs/cenotetaker3/ct3_DBs
    prune: 'T'
  virsorter2:
    env: /path/to/envs/virsorter2-2.2.4-single
    db_dir: /path/to/dbs/virsorter2/db
  vrhyme:
    env: /path/to/envs/vrhyme-1.1.0
    min_length: 2000     # vRhyme's own floor; lower is refused
    mapper: ''           # which mapper's coverage table; blank = first
```

Per run, without editing a config:

```bash
./run_magmaker.sh --profile resources/profiles/demon \
  --virus-predictors=genomad,vibrant      # replaces the predictor switches
./run_magmaker.sh --profile resources/profiles/demon --skip-vrhyme
```

### geNomad score calibration

`params.genomad.score_calibration` (on by default) passes
`--enable-score-calibration`. geNomad computes FDR **only** when scores are
calibrated, and calibration is off in the tool's own default, so without this
the `FDR` column of `virus_summary.tsv` is always empty.

It is not purely additive. Calibration also activates geNomad's `--max-fdr`
filter — 0.10 in its default preset, 0.05 under `--conservative`, 1.00 under
`--relaxed` — which is ignored when scores are not calibrated. Turning it on
therefore both fills in `FDR` and drops the calls that fail that threshold,
which is the behavior geNomad's defaults are built around. Set it to `false`
to get the unfiltered, FDR-less call set back.

geNomad disables calibration itself on an input of fewer than 1,000
sequences, because calibration estimates the sample's composition. `FDR` is
then empty for that assembly and no FDR filtering happens. The rule's
`params.txt` records whether the flag was passed, and geNomad's own log says
whether calibration ran, so an empty `FDR` is always explainable.

### Why four of these name an environment instead of a conda yaml

geNomad and CheckV install from `resources/env/*.yaml` like everything else.
The other four cannot: their dependencies are too fragile to re-solve, so they
run in an environment that already exists.

- **VIBRANT** pins `scikit-learn==0.21.3`, which uses `np.float`, removed in
  numpy 1.24. The environment works only because its numpy still deprecates
  rather than removes it. Never update numpy in it.
- **VirSorter2** needs one environment satisfying its own per-rule spec
  (`virsorter=2.2.4`, `python=3.8`, `scikit-learn=0.22.1`, `screed=1`,
  `pandas=1.2`, `numpy<1.24`, `prodigal=2.6`, `hmmer!=3.3.1`), and runs with
  `--use-conda-off`. Its per-rule environments never activate under a
  non-interactive shell, because conda is an interactive shell function and
  rules run in `bash -c`: the rules then fall through to whatever python is on
  PATH and die on `import screed`. The pins are kept and the activation
  removed, not the other way round. `--use-conda-off` against any other
  environment would classify under a scikit-learn the model was not trained
  with, turning a crash into a silent wrong answer, so the rule checks the
  version and fails if it is not 0.22.x.
- **vRhyme** ships a pickled ExtraTrees model built with an old scikit-learn.
  Version 1.3 added `missing_go_to_left` to tree nodes and cannot load it, and
  the README's floor of `>=0.23.0` has no ceiling, so conda takes the latest
  and the tool breaks. Pin `scikit-learn=1.2.2` and check the model unpickles
  at install time: on real data vRhyme never reaches the model load when it
  has no candidate pairs, so it exits 0 reporting success with 0 bins.
- **Cenote-Taker3** is less fragile, but is kept in the same place for
  consistency with the other three.

On the demon cluster these are already built. Elsewhere, build them the same
way — one conda prefix per tool, exact pins, `conda list --explicit` saved
beside each — and point `params.<tool>.env` at them.

---

## Output

`output/virus/`:

```
genomad/{assembler}/{row}/           geNomad output (summary kept, intermediates removed)
vibrant/{assembler}/{row}/           VIBRANT output, with params.txt
cenotetaker3/{assembler}/{row}/      Cenote-Taker3 output, with params.txt
virsorter2/{assembler}/{row}/        VirSorter2 output, with params.txt
checkv/{assembler}/{row}/            CheckV for geNomad's calls
checkv_{predictor}/{assembler}/{row}/   CheckV for each other predictor's calls
vrhyme/{assembler}/{row}/            vRhyme inputs, prep_params.txt and results
virus_summary.tsv                    one row per call per predictor
virus_name_mapping.tsv               per predictor and assembly: how names mapped
vrhyme_summary.tsv                   one row per vMAG
vrhyme_membership.tsv                one row per contig, naming its vMAG
```

A predictor that finishes without writing its viral FASTA gets an empty one
plus a `NO_VIRAL_FASTA` note, and `virus_name_mapping.tsv` records
`Output_Present = no`. "Wrote nothing" and "called nothing" stay apart.

Every predictor rule clears its output folder before running, as geNomad's
does: all of them resume from whatever they find, so a retry after a killed
attempt would otherwise mix new results with half-written ones.

### Sequence names

Every tool renames sequences. **Nothing downstream should parse a
tool-specific name**: `virus_summary.tsv` carries `Predictor` and
`Source_Contig`, which is the contig in the assembly the call came from.

| Predictor | File read | Rule |
|---|---|---|
| geNomad | `<stem>_virus.fna` | first whitespace token, then cut at `\|` (which strips `\|provirus_<start>_<end>`) |
| VIBRANT | `<stem>.phages_combined.fna` | first whitespace token |
| Cenote-Taker3 | `ct3_virus_summary.tsv` | first token of `input_name` |
| VirSorter2 | `final-viral-combined.fa` | cut at `\|\|`, then first whitespace token |

The traps each rule exists for:

- **VIBRANT** appends `_fragment_1` to the END of the header, after `len=`,
  not to the contig name:
  `>k141_10643 flag=0 multi=22.9958 len=231935_fragment_1`. Splitting on `_`
  mangles it.
- **Cenote-Taker3** renames beyond recovery (`>ct3_99@C3`), so its FASTA
  headers are never parsed; the join goes through its own summary table,
  which carries both its name and the full original header.
- **geNomad** puts provirus coordinates on the name itself, so without the
  `|` cut no provirus ever joins its contig.
- **VirSorter2** appends `||full`, `||{i}_partial` and `||lt2gene`.

A wrong rule does not fail; it produces a small overlap that reads like
genuine disagreement between tools. So `virus_name_mapping.tsv` reports, per
predictor and assembly, the sequences seen, how many mapped to a source
contig, how many did not, and how many mapped to a name that is not in the
assembly, with examples. A non-zero count in either of the last two means the
rule is wrong.

CheckV refuses a FASTA with a repeated sequence ID, which VIBRANT produces
whenever one contig yields several phage regions. Repeats are renamed
`<id>__dup<N>` rather than dropped, since each is a distinct predicted region,
and the suffix is stripped again when names are mapped.

### Provirus columns

Five, not one boolean, because the provenance is what prophage work needs.

| Column | Meaning |
|---|---|
| `Provirus_Predictor` | the predicted name carried provirus coordinates |
| `Provirus_Predictor_Coords` | that start and end, in the ORIGINAL contig — what locates the integration site and its flanking host sequence |
| `provirus` | CheckV's own call (kept under its old name) |
| `Provirus_CheckV_Coords` | CheckV's trimmed region — what was actually scored |
| `Provirus_Any` | the OR of the two, for filtering |

CheckV reports `No` for a sequence geNomad already excised, because no host DNA
is left to detect. So `Provirus_Predictor = Yes` with `provirus = No` is the
**expected** combination for an excised provirus, not a disagreement. With one
boolean that distinction is unrecoverable: on an earlier vervet catalog CheckV
alone marked 339 of 3,910 medium-or-better genomes (8.7%) as proviruses, while
geNomad OR CheckV marked 1,157 of 3,910 (29.6%).

### vRhyme

`vrhyme_prep` selects which of geNomad's viral contigs to bin:

- `Element_Type == virus` — plasmids are a different question.
- **not a provirus**, by geNomad's name OR CheckV's call. A dormant prophage's
  coverage tracks its host, so vRhyme would bin it with whatever host contigs
  are nearby. Because CheckV says `No` for the ones geNomad excised, the CheckV
  column alone misses most of them.
- at least `params.vrhyme.min_length` (2000, vRhyme's floor; lower is refused).
- **every CheckV quality tier**, `Low-quality` and `Not-determined` included.
  The pieces of a fragmented phage each score badly alone, and joining them is
  the point. Quality is judged after binning, on the vMAGs.

`prep_params.txt` records every count with its denominator, including how many
selected names were found in the coverage table. A poor join invalidates the
binning, so it is reported rather than absorbed.

**Headers are truncated to their first token, and this is essential.** vRhyme
keys scaffolds on the ENTIRE header line, while the coverage table carries
only the contig name. With full headers not one coverage row matches any
scaffold: vRhyme finds no candidate pairs, writes a log saying "vRhyme binning
complete", reports 0 bins and **exits 0**. That cost a day of 13 samples
across six `--cov` and four `--penalty_n` settings, all returning exactly 0
bins, before truncation took one sample from 0 bins to 473. Upstream's example
data hides it, because its headers have no description.

Two more things about vRhyme:

- It is **not deterministic**: it keeps the best of 20 iterations and has no
  seed option, so an identical input can give a slightly different answer.
  Record the run; do not assume a rerun reproduces it.
- It sets **circular scaffolds** aside before binning (DTR/ITR detection).
  Those are complete genomes, so `vrhyme_summary.tsv` lists them as their own
  single-contig vMAGs marked `Circular`, rather than letting them go missing.

The number that decides whether vRhyme earned its place is how many contigs
ended up in a **multi-contig** bin, since a bin of one contig is just that
contig renamed. `vrhyme_summary.tsv` carries `Multi_Contig` and `N_Contigs`
per vMAG, and the rule's log states both counts with their denominators.

---

## Agreement between the predictors

Measured on 4 of 88 vervet samples, MEGAHIT only, on assemblies host-filtered
against a host reference since found to be contaminated. VirSorter2 ran on 1
of those 4 (the smallest). **These describe tool behavior, not biology, and
are not results.** They are here because they are why only geNomad is on by
default.

All four name-mapping rules were clean on that run: zero unmapped sequences
and zero names absent from the assembly, on every sample.

Contigs called, as a share of assembly contigs:

| Sample | Assembly contigs | geNomad | VIBRANT | Cenote-Taker3 | VirSorter2 |
|---|---|---|---|---|---|
| A | 10,575 | 108 (1.02%) | 31 (0.29%) | 29 (0.27%) | 193 (1.83%) |
| B | 243,698 | 2,058 (0.84%) | 674 (0.28%) | 731 (0.30%) | not run |
| C | 608,114 | 5,535 (0.91%) | 1,743 (0.29%) | 1,972 (0.32%) | not run |
| D | 642,399 | 7,920 (1.23%) | 2,759 (0.43%) | 2,857 (0.44%) | not run |

Agreement with geNomad, as a share of each tool's own calls: VIBRANT 68.2 to
87.1%, Cenote-Taker3 64.6 to 82.8%, VirSorter2 40.4% on its one sample.
Across the three large samples, 1,563 of VIBRANT's calls (30.2% of its own)
and 1,920 of Cenote-Taker3's (34.5%) were contigs geNomad did not call.

About 70% of all called contigs were called by exactly one tool and about 10%
by all three, with pairwise Jaccard 20 to 25%. Only 1.0 to 2.1% of assembly
contigs were called by any tool.

**What this does not answer:** whether a tool's unique calls are genuine
genomes or low-quality fragments. That is what per-predictor CheckV is for,
and it is the number an adoption decision rests on: per tool, how many contigs
it calls that geNomad did not reach medium quality or better.

### Runtimes

Same 4 samples, 16 CPUs and 128 GB each, 16 jobs running at once, so these
are contended rather than clean timings. Runtime scales with assembly size,
not with how many viruses are found.

| Assembly | geNomad | VIBRANT | Cenote-Taker3 | VirSorter2 |
|---|---|---|---|---|
| 10.3 Mb | 2:49 | 1:51 | 2:09 | 27:02 |
| 297 Mb | 37:55 | 22:24 | 13:18 | >81 min |
| 592 Mb | >81 min | 34:12 | 23:30 | >81 min |
| 620 Mb | >81 min | 36:41 | 26:07 | >81 min |

VirSorter2 is roughly ten times the others on the smallest assembly and had
not finished any of the three larger ones at 81 minutes, which is why its rule
asks for 48 hours. Peak memory was not captured reliably on that run; size
later jobs from `sacct MaxRSS` on the `.batch` step.

---

## Scope

This is a **per-sample** track. Two steps are deliberately downstream of it:

- **vOTU dereplication** across samples (95% ANI / 85% AF, MIUViG), kept out
  for the same reason dRep is. A catalog pools genomes from every run folder
  at once, so it belongs to no single run, and adding a study can change every
  cluster and representative. It is an analysis product versioned with the
  paper, not an artifact versioned with a pipeline run.
- **Gene-sharing network taxonomy** (vConTACT3), which runs on that catalog's
  representatives rather than on every member genome. It is wanted because
  geNomad leaves most calls at class level, but it has nothing to work on until
  the catalog exists.

The contract between this track and whatever does those steps is
`virus_summary.tsv`: one row per call per predictor, with `Predictor` and
`Source_Contig`.

---

← [Documentation home](index.md)
