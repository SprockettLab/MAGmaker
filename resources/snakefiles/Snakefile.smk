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

kraken_targets = (
    expand("output/profile/kraken2/{sample}.report.txt", sample=samples)
    if 'kraken2' in config.get('profilers', []) else []
)

rule all:
    input:
        "output/qc/multiqc/multiqc.html",
        "output/assemble/multiqc_assemble/multiqc.html",
        "output/prototype_selection/sourmash_plot",
        "output/prototype_selection/prototype_selection/selected_prototypes.yaml",
        "output/profile/metaphlan/merged_abundance_table.txt",
        *kraken_targets
