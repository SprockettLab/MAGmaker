# Host index, BoarderCollie on: BoarderCollie's shared bowtie2 index.
# Included by host_reference.smk, which sets host_reference, host_base,
# BOARDERCOLLIE_SNAKEFILE and BOARDERCOLLIE_CONFIG (with reference_out and
# bowtie2_out pointing into output/qc/host_reference/).
#
# BoarderCollie builds the index once per host key under its index_dir
# (bt2/<host>/<host key>/) and hard-links it to host_base.*.bt2l. Any other run
# folder with the same host links to the same files; the build runs only when
# that host key has no index yet. A changed host (other genomes, other
# settings) has another key, so it gets its own index, never a stale one.

module boardercollie:
    snakefile: BOARDERCOLLIE_SNAKEFILE
    config: BOARDERCOLLIE_CONFIG

use rule * from boardercollie as boardercollie_*

localrules: host_bowtie2_build


rule host_bowtie2_build:
    # The name host_filter asks for. Nothing is built here: the index is the
    # shared one BoarderCollie linked to host_base; this writes the record of
    # which reference and which shared index host_filter uses.
    input:
        index=host_base + ".index_done.txt",
        # ancient: a hard link of BoarderCollie's shared library FASTA, which
        # other run folders with this host touch (see host_index_boardercollie.smk).
        reference=ancient(host_reference),
        members=host_reference + ".members.tsv"
    output:
        host_base + ".bowtie2_build.done"
    params:
        indexbase=host_base
    log:
        "output/logs/qc/host_bowtie2_build/host_bowtie2_build.log"
    shell:
        """
        cp {input.members} {params.indexbase}.members.tsv
        {{ printf 'reference\\t%s\\nindex\\t%s\\n' {input.reference} {params.indexbase}; cat {input.index}; }} > {output}
        echo "Shared BoarderCollie bowtie2 index linked at {params.indexbase}; nothing built in this rule." > {log}
        """
