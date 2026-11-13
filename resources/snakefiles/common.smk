# =============================================================================
# common.smk — sample metadata loading, shared by Snakefile and Snakefile-bin
#
# Both entry points need the same metadata, so the loader lives here rather
# than being duplicated. This file defines no rules; it only exposes:
#
#   metadata_table   pandas DataFrame indexed by (Sample, Sequencing_Run);
#                    short-read (Platform illumina) rows only
#   samples          unique sample names, in file order
#   long_read_table  the Platform ont rows, same index
#   long_read_samples   samples with long reads, in file order
#   long_reads_for()    sample -> list of its long-read FASTQ paths
#   seqruns          the (Sample, Sequencing_Run) MultiIndex
#   reads            read identifiers from config, normally ['R1', 'R2']
#   get_read()       (sample, seqrun, read) -> FASTQ path
#   seqruns_for()    sample -> list of its sequencing runs
#   is_single_end()  sample -> True for a single-end sample
#   reads_for()      sample -> ['R1', 'R2'] or ['SE']
#   nonhost_reads()  sample -> host-filtered FASTQ path(s)
#   samples_for_assembler()  assembler -> samples it can assemble
#
# Columns beyond the four required ones are carried through untouched and can
# be read from metadata_table inside a rule, so study covariates live
# alongside the file paths instead of in a separate sheet.
#
# PLATFORM
#
# Every row states its sequencing platform in the required Platform column:
# illumina for short reads, ont for Oxford Nanopore long reads. A long-read
# run is one more row of the same sample, with its FASTQ in R1_fp and the
# literal NA in R2_fp. Nothing is detected from the reads. Short-read steps
# see only the illumina rows, so a sheet without ont rows runs exactly as
# before; ont rows go to long_read_table.
#
# SINGLE-END INPUT
#
# A row declares single-end sequencing by putting the literal NA in R2_fp.
# Nothing is inferred: not from the filesystem, not from read headers, not
# from a missing column. An empty, absent, "nan" or "None" R2_fp stays a
# hard error exactly as before, so an incomplete or half-written sample
# sheet still fails loudly instead of quietly assembling single-end.
#
# Single-end samples carry the read identifier 'SE' instead of 'R1'/'R2'.
# Because every QC rule is already parameterized by read identifier, that
# one substitution carries most of the way through the pipeline, and it
# keeps single-end outputs from ever colliding with paired-end ones.
# =============================================================================

import os
import sys
import pandas as pd

REQUIRED_COLUMNS = ["Sample", "Sequencing_Run", "Platform", "R1_fp", "R2_fp"]

# Platform values, exact and lowercase.
SHORT_READ_PLATFORM = "illumina"
LONG_READ_PLATFORM = "ont"
PLATFORMS = [SHORT_READ_PLATFORM, LONG_READ_PLATFORM]

# Separates a sample from an assembler in long-read binning rows
# (<sample>__flye), so it may not occur in the name of a long-read sample.
ROW_SEPARATOR = "__"

# The literal that declares single-end input. Case-sensitive and exact:
# "na", "n/a" and "None" are all still errors.
SINGLE_END_MARKER = "NA"

# Read identifier used for single-end samples throughout the pipeline.
SINGLE_END_READ = "SE"

# Values that mean "this field was left blank", and are always an error.
BLANK_VALUES = ["", "nan", "NaN", "None"]


# The sample and seqrun wildcard constraints live in Snakefile and
# Snakefile-bin, right after the include of this file: a directive here
# would stop this file being plain Python, and the metadata tests exec it.


def _fail(message):
    """Abort with a message readable without a Python traceback."""
    sys.exit("\nMAGmaker metadata error\n" + "-" * 60 + "\n" + message + "\n")


def load_metadata(config):
    """in : the snakemake config dict
       out: (metadata_table, samples, long_read_table)

    Validation happens up front so a malformed sheet fails immediately with a
    specific message, instead of surfacing later as a missing-input or
    KeyError somewhere deep in the DAG."""

    if "metadata" not in config:
        # Point at the migration rather than just reporting a missing key.
        if "samples" in config or "units" in config:
            _fail(
                "config.yaml still uses the old 'samples:' and 'units:' keys.\n"
                "They are replaced by one metadata table:\n\n"
                "    metadata: resources/config/metadata.txt\n\n"
                "The metadata file is units.txt with renamed columns:\n"
                "    Sample  Sequencing_Run  R1_fp  R2_fp  [covariates...]\n\n"
                "samples.txt is no longer needed; the sample list comes from\n"
                "the Sample column. To convert an existing units.txt:\n\n"
                "    sed '1s/.*/Sample\\tSequencing_Run\\tR1_fp\\tR2_fp/' \\\n"
                "        units.txt > metadata.txt\n\n"
                "See docs/configuration.md."
            )
        _fail("config.yaml must define 'metadata' (path to the metadata table).")

    metadata_fp = config["metadata"]

    try:
        # keep_default_na=False is load-bearing: pandas otherwise converts the
        # literal NA to NaN, which is exactly the value used to declare
        # single-end input, and it would arrive here indistinguishable from a
        # genuinely empty cell. It also stops NA in any carried-through
        # covariate column from being silently rewritten.
        table = pd.read_csv(metadata_fp, sep="\t", header=0, dtype=str,
                            keep_default_na=False)
    except FileNotFoundError:
        _fail("metadata file not found: %s" % metadata_fp)
    except Exception as exc:
        _fail("could not read metadata file %s\n  %s" % (metadata_fp, exc))

    # Tabs are easily lost when a sheet is pasted through a terminal or
    # reformatted by an editor; a single-column frame is the signature.
    if table.shape[1] == 1:
        _fail(
            "metadata file %s parsed as a single column.\n"
            "It must be TAB-separated. Check with:  cat -A %s | head -2\n"
            "(real tabs appear as ^I)" % (metadata_fp, metadata_fp)
        )

    table.columns = [str(c).strip() for c in table.columns]

    if "Platform" not in table.columns and all(
            c in table.columns for c in REQUIRED_COLUMNS if c != "Platform"):
        _fail(
            "metadata file %s has no Platform column.\n\n"
            "Every row must state its sequencing platform: illumina for short\n"
            "reads, ont for Oxford Nanopore long reads. For a sheet of short\n"
            "reads only, add the column with illumina on every row:\n\n"
            "    awk 'BEGIN{FS=OFS=\"\\t\"} NR==1{print $0, \"Platform\"; next}"
            " {print $0, \"illumina\"}' \\\n"
            "        %s > metadata_with_platform.txt\n\n"
            "See docs/configuration.md." % (metadata_fp, metadata_fp)
        )

    missing = [c for c in REQUIRED_COLUMNS if c not in table.columns]
    if missing:
        _fail(
            "metadata file %s is missing required column(s): %s\n"
            "Required: %s\n"
            "Found   : %s"
            % (metadata_fp, ", ".join(missing),
               ", ".join(REQUIRED_COLUMNS), ", ".join(table.columns))
        )

    # fillna before astype: in pandas >=3 a missing value survives .astype(str)
    # as a float NaN rather than the string "nan", so comparing to "" alone
    # would silently miss empty fields.
    for column in REQUIRED_COLUMNS:
        table[column] = table[column].fillna("").astype(str).str.strip()

    # R2_fp is checked separately below, because the literal NA is meaningful
    # there and an error everywhere else.
    always_required = [c for c in REQUIRED_COLUMNS if c != "R2_fp"]

    blank = table[
        table[always_required].isin(BLANK_VALUES + [SINGLE_END_MARKER]).any(axis=1)
    ]
    if not blank.empty:
        # +2 converts a 0-based row index to the line number in the file.
        lines = ", ".join(str(i + 2) for i in blank.index)
        _fail(
            "metadata file %s has empty required fields on line(s): %s\n"
            "Every row needs Sample, Sequencing_Run, Platform and R1_fp."
            % (metadata_fp, lines)
        )

    # A blank R2_fp is still an error. Only the exact literal NA declares
    # single-end, so a truncated sheet cannot silently downgrade a paired
    # library to single-end and assemble half the data.
    blank_r2 = table[table["R2_fp"].isin(BLANK_VALUES)]
    if not blank_r2.empty:
        lines = ", ".join(str(i + 2) for i in blank_r2.index)
        _fail(
            "metadata file %s has a blank R2_fp on line(s): %s\n\n"
            "R2_fp is never optional. Give the path to the reverse reads,\n"
            "or write the literal\n\n"
            "    NA\n\n"
            "to declare that the run is single-end. Declaring it explicitly\n"
            "is required so that an incomplete sample sheet cannot be\n"
            "mistaken for single-end data." % (metadata_fp, lines)
        )

    bad_platform = table[~table["Platform"].isin(PLATFORMS)]
    if not bad_platform.empty:
        lines = ", ".join(str(i + 2) for i in bad_platform.index)
        _fail(
            "metadata file %s has an unknown Platform on line(s): %s\n"
            "Platform must be one of: %s (exact, lowercase)."
            % (metadata_fp, lines, ", ".join(PLATFORMS))
        )

    # A repeated (Sample, Sequencing_Run) makes the index ambiguous, across
    # platforms too, so it is checked before the rows are split.
    duplicated = table.duplicated(subset=["Sample", "Sequencing_Run"], keep=False)
    if duplicated.any():
        pairs = (table.loc[duplicated, ["Sample", "Sequencing_Run"]]
                 .drop_duplicates().itertuples(index=False, name=None))
        _fail(
            "metadata file %s has duplicate (Sample, Sequencing_Run) pairs:\n  %s\n"
            "Each row must describe one sequencing run of one sample."
            % (metadata_fp, "\n  ".join("%s / %s" % p for p in pairs))
        )

    long_rows = table[table["Platform"] == LONG_READ_PLATFORM].copy()
    table = table[table["Platform"] == SHORT_READ_PLATFORM].copy()

    paired_long = long_rows[long_rows["R2_fp"] != SINGLE_END_MARKER]
    if not paired_long.empty:
        lines = ", ".join(str(i + 2) for i in paired_long.index)
        _fail(
            "metadata file %s has an ont row with an R2_fp on line(s): %s\n"
            "A long-read run has one FASTQ, in R1_fp; R2_fp must be NA."
            % (metadata_fp, lines)
        )
    no_short = sorted(set(long_rows["Sample"]) - set(table["Sample"]))
    if no_short:
        _fail(
            "metadata file %s has long reads for sample(s) with no illumina\n"
            "rows: %s\n\n"
            "Long-read assemblies are polished with the sample's short reads\n"
            "and binned with short-read coverage, so a sample with long reads\n"
            "needs short reads too. Leave samples without short reads (e.g.\n"
            "negative controls) out of the sheet." % (metadata_fp, ", ".join(no_short))
        )
    separator = sorted(s for s in set(long_rows["Sample"]) if ROW_SEPARATOR in s)
    if separator:
        _fail(
            "metadata file %s: long-read sample name(s) contain %r: %s\n"
            "%r separates a sample from its assembler in binning rows\n"
            "(e.g. v1928__flye), so it may not occur in these names."
            % (metadata_fp, ROW_SEPARATOR, ", ".join(separator), ROW_SEPARATOR)
        )

    table["Layout"] = table["R2_fp"].apply(
        lambda v: "SINGLE" if v == SINGLE_END_MARKER else "PAIRED"
    )

    # merge_seqruns concatenates per read identifier, so a sample whose runs
    # disagree about layout would produce a reverse file covering only some
    # of its runs. Caught here rather than as a short R2 much later.
    mixed = (table.groupby("Sample")["Layout"].nunique() > 1)
    if mixed.any():
        names = ", ".join(sorted(mixed[mixed].index))
        _fail(
            "metadata file %s mixes single-end and paired-end runs within\n"
            "the same sample: %s\n\n"
            "Every sequencing run of one sample must have the same layout."
            % (metadata_fp, names)
        )

    metadata_table = table.set_index(["Sample", "Sequencing_Run"])
    metadata_table.sort_index(inplace=True)
    long_read_table = long_rows.set_index(["Sample", "Sequencing_Run"])
    long_read_table.sort_index(inplace=True)

    # Order-preserving unique, so run order follows the file rather than sort
    # order -- keeps logs and reports in the order the user wrote them.
    samples = list(dict.fromkeys(table["Sample"]))

    return metadata_table, samples, long_read_table


metadata_table, samples, long_read_table = load_metadata(config)

# Samples with long reads, in file order (a subset of samples).
long_read_samples = [s for s in samples
                     if s in set(long_read_table.index.get_level_values("Sample"))]


def long_reads_for(sample):
    """in : sample name
       out: its long-read FASTQ paths, one per ont row, in index order"""
    return list(long_read_table.loc[sample, "R1_fp"])


def _validate_prototype_exclusions(config, samples):
    """Check params.prototypes.exclude against the sample list, at load.

    Validated here rather than inside the prototype_selection rule so a
    typo stops the run immediately instead of after sketching and
    comparing every sample, and so the failure names a metadata problem
    rather than surfacing as a Python traceback mid-DAG."""
    excluded = (config.get("params", {}).get("prototypes", {}).get("exclude")
                or [])
    if isinstance(excluded, str):
        _fail("params.prototypes.exclude must be a list of sample names, "
              "not a single string. Write it as:\n"
              "    exclude:\n      - %s" % excluded)
    unknown = [e for e in excluded if e not in samples]
    if unknown:
        _fail(
            "params.prototypes.exclude names sample(s) that are not in the "
            "metadata table:\n"
            "  %s\n\n"
            "Known samples:\n  %s"
            % (", ".join(unknown), ", ".join(samples)))
    return set(excluded)


prototype_exclusions = _validate_prototype_exclusions(config, samples)


def cmseq_enabled():
    """Whether strain heterogeneity is computed for this run.

    Defined here rather than in cmseq.smk because mag_qc.smk calls it from
    an input function, and relying on one snakefile's functions being
    visible to another makes correctness depend on include order. common.smk
    is included first by both Snakefile and Snakefile-bin, so anything
    defined here is available everywhere."""
    return str(config.get('skip_cmseq', '')).strip().lower() not in (
        'true', '1', 'yes')


def vsc_enabled():
    """Whether MetaPhlAn also profiles Viral Sequence Clusters.

    On by default. VSC profiling costs one extra bowtie2 pass of the same
    reads inside the metaphlan job, and needs no database beyond the
    {index}_VSG.fna and {index}_VINFO.csv that ship with the CHOCOPhlAn
    index already in use.

    It is a toggle rather than always-on because of one failure mode:
    metaphlan exits non-zero when a sample yields no viral mappings at all
    (metaphlan.py, "No viral mappings in the viral sam file"). That would
    take down the whole metaphlan job for that sample, taxonomic profile
    included, not just the viral part. A host whose viruses are poorly
    represented in the VSC database is the case where this bites. Three
    vervet samples of 88 (the smallest, median and largest by read volume)
    were tested on 2026-11-10 and all three produced 6 to 16 clusters, so
    the abort was never provoked -- which is not the same as showing it
    cannot happen on the other 85.

    Set params.metaphlan.profile_vsc: false in the config, or pass
    --skip-vsc for one run."""
    if str(config.get('skip_vsc', '')).strip().lower() in ('true', '1', 'yes'):
        return False
    return str(
        config.get('params', {}).get('metaphlan', {}).get('profile_vsc', True)
    ).strip().lower() not in ('false', '0', 'no')


# Which tool reconciles the three per-binner bin sets into one set of MAGs,
# and the directory its selected FASTAs land in. Everything downstream of
# selection -- checkm2, gunc, gtdbtk, cmseq, mag_summary -- reads that one
# directory, so switching tools moves a single path rather than editing
# every consumer.
#
# DAS_Tool's directory name is unchanged, so an existing run keeps every
# path it already has and nothing reruns on upgrade.
CONSOLIDATION_TOOLS = {
    'das_tool': 'DAS_Tool_Fastas',
    'binette': 'Binette_Fastas',
}


def consolidation_tool():
    """in : nothing; reads config
       out: 'das_tool' or 'binette'

    Settable two ways, because the two have different lifetimes: a config
    file records what an arm is for, while --config consolidation_tool=...
    overrides it for one invocation without editing a tracked file, which
    is what comparing the two tools on one arm needs.

    Validated here rather than where it is consumed so a typo stops the run
    at DAG construction, naming the valid values, instead of surfacing much
    later as an empty bins directory that looks like a sample with no MAGs.
    """
    tool = str(
        config.get('consolidation_tool')
        or config.get('params', {}).get('consolidation', {}).get('tool')
        or 'das_tool'
    ).strip().lower()
    if tool not in CONSOLIDATION_TOOLS:
        _fail(
            "unknown consolidation tool: %r\n\n"
            "Valid values are: %s\n\n"
            "Set it in config.yaml as\n\n"
            "    params:\n      consolidation:\n        tool: binette\n\n"
            "or for a single run as\n\n"
            "    --config consolidation_tool=binette"
            % (tool, ", ".join(sorted(CONSOLIDATION_TOOLS)))
        )
    return tool


SELECTED_FASTAS = CONSOLIDATION_TOOLS[consolidation_tool()]

# Binners whose bins can be handed to a consolidation tool, in a fixed
# order so labels are stable between runs. Which of them actually run is
# `binners:` in the config; this is the set the pipeline knows how to wire.
KNOWN_BINNERS = ['metabat2', 'maxbin2', 'concoct', 'semibin2']


def enabled_binners():
    """The configured binners, in KNOWN_BINNERS order.

    Ordered here rather than taken as the config lists them, so that two
    runs whose config differs only in the order of `binners:` produce the
    same labels in the same order.
    """
    configured = config.get('binners', []) or []
    if isinstance(configured, str):
        configured = [configured]
    # A name this pipeline does not know is a hard error, not a binner
    # quietly left out. Silently dropping `maxbinn2` would run three binners
    # where four were asked for, produce a smaller MAG set, and say nothing
    # -- and the result would look like a real answer.
    unknown = [b for b in configured if b not in KNOWN_BINNERS]
    if unknown:
        _fail(
            "config `binners:` names binner(s) this pipeline does not know:\n"
            "  %s\n\n"
            "Known binners: %s\n\n"
            "Adding a new one means a rule producing\n"
            "  output/selected_bins/<binner>/<mapper>/scaffolds2bin/"
            "<sample>_scaffolds2bin.tsv\n"
            "and an entry in KNOWN_BINNERS and BINNER_BIN_DIRS in common.smk."
            % (", ".join(unknown), ", ".join(KNOWN_BINNERS)))
    if not configured:
        _fail("config `binners:` is empty; at least one binner is required.")
    return [b for b in KNOWN_BINNERS if b in configured]


def binner_of_table(path):
    """in : a .../selected_bins/<binner>/<mapper>/scaffolds2bin/... path
       out: the binner name

    Read back out of the path rather than carried alongside it. Both
    consolidation tools need bin sets paired with their labels, and pairing
    two lists by position breaks silently the moment one binner declines a
    sample and the lists stop lining up.
    """
    parts = os.path.normpath(path).split(os.sep)
    try:
        return parts[parts.index('selected_bins') + 1]
    except (ValueError, IndexError):
        return 'unknown'


# Long-read assemblers (long_reads.smk); their binning rows also get the
# circular-contig candidate set.
LONG_READ_ASSEMBLER_NAMES = ("flye", "myloasm")


def is_long_read_row(row):
    """Overridden in Snakefile-bin, which knows the binning table."""
    return False


def binner_tables(wildcards):
    """Every enabled binner's contig-to-bin table for one binning row, plus,
    for a long-read row, the circular-contig candidates (one bin per closed
    contig; long_read_binning.smk)."""
    tables = expand(
        "output/selected_bins/{binner}/{mapper}/scaffolds2bin/"
        "{contig_sample}_scaffolds2bin.tsv",
        binner=enabled_binners(),
        mapper=config['mappers'],
        contig_sample=wildcards.contig_sample)
    if is_long_read_row(wildcards.contig_sample):
        tables += expand(
            "output/selected_bins/circular/{mapper}/scaffolds2bin/"
            "{contig_sample}_scaffolds2bin.tsv",
            mapper=config['mappers'], contig_sample=wildcards.contig_sample)
    return tables


# Where each binner leaves its bins. CONCOCT's differs because its clusters
# are written as a table and turned into FASTAs by a later rule, so the
# directory that holds bins is not the one named after the binner.
BINNER_BIN_DIRS = {
    'metabat2': "output/binning/metabat2/{mapper}/run_metabat2/{contig_sample}/",
    'maxbin2':  "output/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}/",
    'concoct':  "output/binning/concoct/{mapper}/extract_fasta_bins/{contig_sample}_bins/",
    'semibin2': "output/binning/semibin2/{mapper}/run_semibin2/{contig_sample}/",
}


def binner_bin_dirs(samples):
    """Bin directories for every enabled binner, across samples."""
    out = []
    for b in enabled_binners():
        out += expand(BINNER_BIN_DIRS[b],
                      mapper=config['mappers'], contig_sample=samples)
    return out



def config_flag(name):
    """A --config boolean, tolerant of True/true/1/yes.

    Snakemake hands --config values over either parsed or as bare strings
    depending on how they were written, so both are accepted."""
    return str(config.get(name, "")).strip().lower() in ("true", "1", "yes")


def gtdbtk_mode():
    """in : nothing; reads config
       out: 'full', 'align_only' or 'none'

    GTDB-Tk runs in three stages. `identify` calls genes and searches the
    marker HMMs, `align` builds the concatenated alignment, and `classify`
    places each genome in the reference tree with pplacer.

    Nearly all of the cost is the last stage. GTDB-Tk's documentation puts
    the bacterial memory requirement at about 140 GB and attributes it to
    pplacer; this pipeline reserves 128 GB on the first attempt and
    escalates to 320 GB. identify and align do not run pplacer at all, so
    stopping before classify is a change of memory class rather than a
    tuning change: a rule only the largest nodes can schedule becomes one
    that runs anywhere.

    What align produces is the alignment itself, and MSA_Percent -- the
    share of the concatenated alignment a genome fills with residues -- is
    recoverable from it exactly. That is the quantity that decides whether
    a genome can be placed in a tree, so an analysis that filters on
    MSA_Percent but takes its taxonomy from elsewhere never needs classify.

      full        identify + align + classify. The default, unchanged.
      align_only  --skip-gtdbtk-classify. MSA_Percent is populated,
                  taxonomy columns are NA, pplacer never runs.
      none        --skip-gtdbtk. Both are NA.
    """
    if config_flag('skip_gtdbtk'):
        return 'none'
    if config_flag('skip_gtdbtk_classify'):
        return 'align_only'
    return 'full'
seqruns = metadata_table.index
reads = config["reads"]

# sample -> "SINGLE" or "PAIRED". Validated in load_metadata to be constant
# within a sample, so the first run of each is representative.
layout_index = (metadata_table.reset_index()
                .groupby("Sample")["Layout"].first().to_dict())


def get_read(sample, seqrun, read):
    """in : sample name, sequencing run name, read id ('R1', 'R2' or 'SE')
       out: path to that FASTQ

    `read` is the bare identifier so it can double as a filename wildcard;
    the metadata column it maps to is <read>_fp. 'SE' has no column of its
    own -- a single-end run keeps its only FASTQ in R1_fp -- so it reads
    R1_fp while carrying a distinct identifier through output filenames."""
    column = "R1_fp" if read == SINGLE_END_READ else "%s_fp" % read
    return metadata_table.loc[(sample, seqrun), column]


def is_single_end(sample):
    """in : sample name
       out: True if the sample was sequenced single-end

    Layout is validated to be consistent across a sample's runs, so the
    first run answers for all of them."""
    return layout_index[sample] == "SINGLE"


def reads_for(sample):
    """in : sample name
       out: ['SE'] for single-end, otherwise the configured read ids

    The list every per-read expansion should iterate, in place of the
    global `reads`, so single-end samples never have an R2 target."""
    return [SINGLE_END_READ] if is_single_end(sample) else reads


def nonhost_reads(sample):
    """in : sample name
       out: list of host-filtered FASTQ paths, one entry or two

    Rules take this instead of referring to host_filter's outputs directly.
    The paired and single-end host filters are separate rules with separate
    outputs, so a fixed `rules.host_filter.output.nonhost_R2` cannot serve
    both."""
    return expand("output/qc/host_filter/nonhost/{sample}.{read}.fastq.gz",
                  sample=sample, read=reads_for(sample))


def samples_for_assembler(assembler):
    """in : assembler name
       out: the samples that assembler can run on

    metaSPAdes rejects a single-end-only library outright, so single-end
    samples are assembled by MEGAHIT alone. Filtering the target list keeps
    that a planning decision rather than a run-time crash partway through a
    long job."""
    if assembler == "metaspades":
        return [s for s in samples if not is_single_end(s)]
    return list(samples)


def mem_escalate(key, base_default=8000, cap_multiple=4):
    """in : a config['mem_mb'] key, plus fallbacks
       out: a snakemake resource callable of (wildcards, attempt)

    Returns config['mem_mb'][key] on the first attempt and doubles it on each
    retry, capped by config['mem_mb'][key + '_max'].

    mem_mb is a scheduler RESERVATION, not a ceiling: a job requesting 300 GB
    will only run on a node with 300 GB free, however much it actually uses.
    Sizing every request for the worst-observed case therefore confines the
    whole rule to the largest nodes and serialises it, to benefit the few
    samples that need the headroom. Escalating instead keeps the common case
    schedulable anywhere and pays the queueing cost only where it is earned.

    Rules using this need a `retries:` directive, or `attempt` never exceeds 1
    and the escalation is inert."""
    def _resolve(wildcards, attempt):
        base = config["mem_mb"].get(key, base_default)
        cap = config["mem_mb"].get("%s_max" % key, base * cap_multiple)
        return min(cap, base * attempt)
    return _resolve


def runtime_escalate(key, base_default, cap_multiple=3):
    """in : a config['runtime'] key, plus a fallback base value (minutes)
       out: a snakemake resource callable of (wildcards, attempt)

    Mirrors mem_escalate above: returns base minutes on the first attempt
    and base * attempt on each retry, capped by base * cap_multiple
    (override per rule with config['runtime'][key + '_max']).

    Added 2026-08-19 after a run of TIMEOUT failures across unrelated
    rules and arms -- Chlorocebus_sabaeus__Jasinska_2013_tentative/
    run_semibin2, Callithrix_jacchus/run_semibin2, Pan_troglodytes_
    troglodytes/run_concoct, Rhinopithecus_roxellana/megahit and
    sourmash_sketch_reads -- turned out to split into two causes a single
    flat runtime can't tell apart in advance: genuine per-sample
    underprovisioning, and contention. `sacct` on 2026-08-19 confirmed two
    Rhinopithecus_roxellana megahit jobs hit their exact 480 min TIMEOUT,
    but each lost ~2h before megahit's own log printed a single line --
    conda activation / isilon I/O under concurrent-arm load, not assembly.
    A retry can't tell which cause it hit either, but does not need to: a
    fresh attempt lands in a different, often less contended window (this
    already resolved two of the cases above -- Callithrix_jacchus and
    Jasinska_2013_tentative's run_semibin2 -- via the existing retries: 2
    alone, no runtime change), and escalating the ceiling on top means a
    retry that hits genuine data-driven slowness gets more room too,
    instead of failing the same way twice. Every static profile runtime:
    override that existed before this was added is now the base_default
    passed in here, so first-attempt behaviour is unchanged; only retries
    (now added wherever a rule lacked one) get more room.

    Rules using this need a `retries:` directive, or `attempt` never
    exceeds 1 and the escalation is inert -- same caveat as mem_escalate."""
    def _resolve(wildcards, attempt):
        base = config["runtime"].get(key, base_default)
        cap = config["runtime"].get("%s_max" % key, base * cap_multiple)
        return min(cap, base * attempt)
    return _resolve


def disable_inner_retries():
    """Keep retries in the main Snakemake process only.

    Under the SLURM executor each job runs a second Snakemake inside its
    allocation (exec mode "remote"), launched with --retries 0. A rule's own
    `retries:` directive overrides that flag, so the inner process re-ran a
    failed command in the same allocation, at the same memory and time,
    before the main process retried with mem_escalate/runtime_escalate.
    WF18_FullTest_2026 (2026-11-13): every run_cenotetaker3 (retries: 2)
    first attempt shows three srun steps, all OUT_OF_MEMORY at 32 GB, ~3 h
    per job; polish_long_read_assembly (retries: 1) showed two. The main
    process then resubmitted correctly at 64 GB.

    Call at the end of each Snakefile, after every rule is defined. Only
    the remote (in-allocation) process is changed; local runs and the main
    process keep their retries."""
    if workflow.remote_exec:
        for r in workflow.rules:
            r.restart_times = 0


def megahit_runtime(wildcards, input, attempt):
    """Snakemake resource callable for rule megahit's runtime.

    Size-scaled floor from real throughput data (2026-08-19, Pan_
    troglodytes_schweinfurthii completed benchmarks, pure compute time
    from the benchmark directive, not including cluster overhead):
    12.70 GB/258.7 min, 16.92 GB/351.5 min, 21.94 GB/390.8 min -- a linear
    fit of ~77 + 14.3 min/GB, padded to 90 + 15 min/GB for safety margin.

    Never below config['runtime'].get('megahit', 480): that flat value
    already covers every sample seen so far up to ~22 GB, and a pure
    size formula would UNDER-predict a small-but-complex sample --
    Rhinopithecus_roxellana/summer-Wild06, only 2.64 GB, genuinely needed
    ~8h because of k-mer graph complexity at k=21, not data volume (its
    own log spent 2h18m on just the k=21 build-to-assemble transition).
    Size alone does not capture that; the flat floor is what keeps a
    small-but-complex sample from being under-provisioned by a formula
    that only knows about GB.

    Escalates by attempt on top of whichever floor wins, same as
    runtime_escalate, so a retry gets more room whether the sample was
    underprovisioned by size or by complexity."""
    size_gb = getattr(input, 'size_mb', 0) / 1024
    base = max(config['runtime'].get('megahit', 480), 90 + 15 * size_gb)
    cap = config['runtime'].get('megahit_max', base * 3)
    return min(cap, base * attempt)


def host_filter_runtime(wildcards, input, attempt):
    """Snakemake resource callable for host_filter/host_filter_se's runtime.

    Was a flat runtime_escalate('host_filter', base_default=480) until
    2026-08-21, when two Pan_troglodytes_schweinfurthii samples
    (SAMN28679412/SRR19415659, SAMN28679423/SRR19415682) timed out at
    that flat 480 min. ENA's own read/base counts confirmed these are
    genuinely deep at the source (553.8M and 507.5M read pairs, 167.3 and
    153.3 Gbp) -- not corrupted or duplicated data -- roughly 3.5-4x
    deeper than a typical sample in the same arm (143.1M reads, 43.2 Gbp)
    that passes within the flat limit. A flat value escalated only by
    doubling on retry (960 min) risked still being insufficient for an
    outlier this deep, and would waste two full retry cycles reaching
    that point regardless.

    Floor stays at config['runtime'].get('host_filter', 480) -- unchanged
    for the great majority of samples, which fall well inside it. The
    size-scaled term only matters for volume well beyond that.

    Rate (20 min per combined GB of R1+R2, or the single SE file) is a
    deliberately generous first-pass estimate, NOT fit from a completed
    benchmark the way megahit_runtime's 90 + 15/GB was fit from three
    real runs -- there is only the one confirmed-insufficient data point
    here (480 min failed at ~76 GB), no confirmed completion time to
    anchor a precise rate. Revisit once these two samples' own
    benchmark files exist and give something real to fit against.

    Escalates by attempt on top of whichever floor wins, same as
    runtime_escalate."""
    size_gb = getattr(input, 'size_mb', 0) / 1024
    base = max(config['runtime'].get('host_filter', 480), 20 * size_gb)
    cap = config['runtime'].get('host_filter_max', base * 3)
    return min(cap, base * attempt)


def seqruns_for(sample):
    """in : sample name
       out: list of that sample's sequencing run names, in index order"""
    return list(metadata_table.loc[sample].index)
