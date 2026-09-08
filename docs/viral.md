# Viral / plasmid track

A per-sample viral and plasmid discovery track that runs parallel to the MAG
track and does not touch it. It is **opt-in**: `run_magmaker.sh` does not run
it, and it is requested as its own Stage 3 target.

```bash
snakemake --snakefile Snakefile-bin --profile resources/profiles/demon virus_all \
  --config binning=output/config/auto_binning.txt
```

---

## What it does

```
per-sample contigs → geNomad end-to-end → viral contigs + plasmid contigs
viral contigs → CheckV → completeness / contamination
→ output/virus/virus_summary.tsv
```

- **geNomad** (`run_genomad`) classifies viral and plasmid contigs straight
  from each assembly. No read mapping is needed, so this branch depends only
  on the assemblies from Stage 1 — it can run before, during or after the MAG
  track.
- **CheckV** (`run_checkv`) scores completeness and contamination for the
  viral contigs geNomad called. It skips cleanly when a sample has none.
- **`make_virus_summary`** merges every sample's geNomad and CheckV output
  into one table.

geNomad's own `--cleanup` is deliberately **not** used — it races NFS
silly-rename on shared storage and exits non-zero after the real outputs are
already written. The rule lets geNomad finish and then removes the heavy
module intermediates itself, keeping the `*_summary` directory downstream
rules read.

---

## Configuration

```yaml
params:
  genomad:
    db_path: /path/to/genomad_db      # genomad download-database <dir>
    extra: ''                          # e.g. --conservative for higher-precision calls
  checkv:
    db_path: /path/to/checkv-db-v1.5  # checkv download_database <dir>
```

See [Database setup](databases.md) for how to download both.

---

## Output

`output/virus/`:

```
genomad/{assembler}/{sample}/    geNomad per-sample output (summary dir kept, intermediates removed)
checkv/{assembler}/{sample}/     CheckV quality for the viral contigs
virus_summary.tsv                merged table, one row per viral or plasmid contig
```

`virus_summary.tsv` columns — see the *Viral summary table* section of
[Output](output.md). In short: geNomad contig stats and score for every row,
geNomad viral taxonomy on `virus` rows, geNomad plasmid annotations
(`Conjugation_Genes`, `AMR_Genes`) on `plasmid` rows, and CheckV quality /
completeness / contamination joined onto the `virus` rows. CheckV does not
score plasmids, so those fields are blank there.

---

## Scope

This is a **per-sample** track. Dereplicating viruses into vOTUs across
samples (95% ANI / 85% AF, MIUViG) is a downstream cross-sample analysis
choice, kept out of MAGmaker for the same reason dRep is — it is not part of
per-sample genome recovery.

---

← [Documentation home](index.md)
