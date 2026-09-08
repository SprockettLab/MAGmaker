import pandas as pd
from os.path import join

configfile: "config.yaml"

# Loads the metadata table and defines metadata_table, samples, seqruns,
# reads, get_read() and seqruns_for(). Shared with Snakefile-bin so the
# two stages cannot disagree about the sample sheet.
include: "resources/snakefiles/common.smk"

include: "resources/snakefiles/qc.smk"
include: "resources/snakefiles/assemble.smk"
include: "resources/snakefiles/prototype_selection.smk"
include: "resources/snakefiles/profile.smk"

def _flag(name):
    """A --config boolean, tolerant of True/true/1/yes as Snakemake may
    hand it over either parsed or as a bare string."""
    return str(config.get(name, "")).strip().lower() in ("true", "1", "yes")


# Read-level profiling is the most I/O-expensive part of the run for what
# it contributes: kraken2 loads a ~200 GB index for every sample, and an
# analysis that takes its taxonomy from GTDB-Tk on the MAGs never reads
# these profiles. profilers and biobakery already declare which tools run;
# the flags let a single run opt out without editing the config.
_profilers = [] if _flag('skip_kraken') else config.get('profilers', [])
_biobakery = [] if _flag('skip_metaphlan') else config.get('biobakery', [])

kraken_targets = (
    expand("output/profile/kraken2/{sample}.report.txt", sample=samples)
    if 'kraken2' in _profilers else []
)

# metaphlan was required unconditionally even though `biobakery` was
# already there to declare it, so removing it from the list had no effect.
metaphlan_targets = (
    ["output/profile/metaphlan/merged_abundance_table.txt"]
    if 'metaphlan' in _biobakery else []
)

rule all:
    input:
        "output/qc/multiqc/multiqc.html",
        "output/assemble/multiqc_assemble/multiqc.html",
        "output/prototype_selection/sourmash_plot",
        "output/prototype_selection/prototype_selection/selected_prototypes.yaml",
        *metaphlan_targets,
        *kraken_targets
