# Viral and plasmid track

MAGmaker includes a viral and plasmid discovery track that runs [geNomad](https://portal.nersc.gov/genomad/) to identify viral and plasmid contigs directly from the assemblies, then scores viral completeness and contamination with [CheckV](https://bitbucket.org/berkeleylab/CheckV/). It runs alongside the MAG track and does not touch it — geNomad classifies contigs straight from the assembly, so no read mapping or binning is required.

The track is **on by default** in `run_magmaker.sh`. Disable it with `--skip-virus` (for example, if the geNomad/CheckV databases are not configured).

---

## What it produces

For every sample and assembler, geNomad extracts viral and plasmid contigs and classifies them, and CheckV then assesses the viral contigs. The per-sample results are merged into a single table, `output/virus/virus_summary.tsv`, with one row per viral or plasmid contig.

Dereplicating viruses into cross-sample vOTUs (95% ANI / 85% AF, MIUViG) is a downstream, cross-sample analysis choice and is intentionally kept out of MAGmaker, exactly as MAG dereplication (dRep) is.

---

## Databases

The track needs the geNomad and CheckV databases. Download each once to a shared location and point `config.yaml` at it. See [Database setup](databases.md#genomad) for full instructions; in brief:

```bash
mamba create -n db_setup -c conda-forge -c bioconda genomad checkv -y
conda activate db_setup

genomad download-database /your/shared/dbs/genomad/    # creates genomad_db/
checkv  download_database /your/shared/dbs/checkv/     # creates checkv-db-v1.5/

conda deactivate
```

Then set in `config.yaml`:

```yaml
params:
  genomad:
    db_path: /your/shared/dbs/genomad/genomad_db
    extra: ''            # e.g. '--conservative' for higher-precision viral calls
  checkv:
    db_path: /your/shared/dbs/checkv/checkv-db-v1.5
```

If these databases are not configured, run with `--skip-virus`.

---

## Running

### With the wrapper (default)

The viral track is stage 4 of `run_magmaker.sh` and runs automatically:

```bash
# Full pipeline including the viral track
./run_magmaker.sh --profile resources/profiles/demon

# Skip the viral track
./run_magmaker.sh --profile resources/profiles/demon --skip-virus
```

`--skip-virus` is consumed by the wrapper and is not forwarded to Snakemake.

### Standalone

`virus_all` lives in the binning Snakefile (`Snakefile-bin`), so it needs the binning config override to parse, like the other stage-3 targets:

```bash
snakemake --snakefile Snakefile-bin --profile resources/profiles/demon virus_all \
  --config binning=output/config/auto_binning.txt
```

The viral rules are per-sample and do not use the binning groups; the override is only needed so `Snakefile-bin` can be parsed. `output/config/auto_binning.txt` must already exist (it is written by `generate_binning_config` in stage 2).

---

## Output

```
output/virus/
├── genomad/{assembler}/{sample}/{sample}.contigs_summary/
│   ├── {sample}.contigs_virus.fna              viral contig sequences
│   ├── {sample}.contigs_plasmid.fna            plasmid contig sequences
│   ├── {sample}.contigs_virus_summary.tsv      geNomad viral calls
│   └── {sample}.contigs_plasmid_summary.tsv    geNomad plasmid calls
├── checkv/{assembler}/{sample}/
│   ├── quality_summary.tsv                     CheckV completeness/contamination
│   ├── viruses.fna                             CheckV-assessed viral sequences
│   └── proviruses.fna                          CheckV-trimmed proviral sequences
└── virus_summary.tsv                           merged table (all samples)
```

### `virus_summary.tsv`

One row per viral or plasmid contig, across all samples and assemblers. CheckV quality fields are populated for viral rows and blank for plasmid rows (CheckV does not score plasmids). geNomad's conjugation and AMR gene calls are populated for plasmid rows.

| Column | Description |
|---|---|
| `Sample_ID` | Sample the contig came from |
| `Assembler` | Assembler that produced the contig (megahit / metaspades) |
| `Element_Type` | `virus` or `plasmid` |
| `Seq_Name` | Contig name (geNomad) |
| `Length_BP` | Contig length (bp) |
| `Topology` | geNomad topology call (e.g. linear, DTR, provirus) |
| `N_Genes` | Number of predicted genes |
| `N_Hallmarks` | Number of viral hallmark genes (capsid, terminase, portal, etc.) |
| `Score` | geNomad virus or plasmid score |
| `FDR` | geNomad false discovery rate estimate |
| `Taxonomy` | geNomad viral taxonomy (viral rows) |
| `Conjugation_Genes` | Plasmid mobilization/conjugation gene families (plasmid rows) |
| `AMR_Genes` | Antimicrobial resistance gene identifiers (plasmid rows) |
| `checkv_quality` | CheckV tier: Complete / High-quality / Medium-quality / Low-quality / Not-determined |
| `miuvig_quality` | MIUViG quality tier |
| `checkv_completeness` | CheckV completeness estimate (%) |
| `checkv_contamination` | CheckV contamination estimate (%) |
| `provirus` | Whether CheckV flagged the contig as a provirus |
| `warnings` | CheckV warnings |

A common downstream filter for a working viral catalog is `Element_Type == virus` with `checkv_quality` in {Complete, High-quality, Medium-quality} — no length or hallmark cut is needed, since CheckV completeness already does the confidence work and a length cut would drop small complete viruses.

---

## Implementation notes

- geNomad's own `--cleanup` is deliberately not used. On NFS filesystems it can race silly-rename (files deleted while a handle is still open become `.nfs*` entries), so `shutil.rmtree` fails with `OSError ENOTEMPTY` after the real outputs are already written. MAGmaker instead lets geNomad finish and then removes the heavy module intermediates itself, keeping the summary directory and the virus/plasmid FASTAs that downstream steps read.
- CheckV is skipped cleanly for any sample with no viral contigs.

---

← [Back to README](../README.md)
