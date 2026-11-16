# Host reads split by the genome they came from.
#
# A BoarderCollie host can join several genomes (a primate and a parasite,
# say) into the one reference host_filter maps against. host_filter removes
# reads that hit any of them; these rules split that count by member genome,
# from the host BAM host_filter already keeps. Nothing is mapped again.
#
#   host_contig_map       which member genome each reference sequence is in
#   host_genome_counts    one sample: removed pairs (or reads) per member
#   merge_host_genome_counts
#                         output/qc/host_filter/host_genome_counts.tsv and
#                         host_genome_counts_params.txt
#
# Each sample's removed total is checked against bowtie2's own summary in the
# host_filter log, and the BAM's sequence names against the contig map; a
# mismatch fails the job. Column definitions: docs/output.md.
#
# Included by Snakefile.smk after qc.smk, which defines host_reference,
# host_bowtie2_build and USE_BOARDERCOLLIE.

HOST_GENOME_COUNTS_SCRIPT = os.path.normpath(os.path.join(
    workflow.current_basedir, "..", "scripts", "host_genome_counts.py"))
HOST_GENOME_COUNTS_DIR = "output/qc/host_filter/genome_counts"
HOST_GENOME_COUNTS_MIN_MAPQ = int(
    config['params'].get('host_genome_counts', {}).get('min_mapq', 10))

if USE_BOARDERCOLLIE:
    HOST_NAME = BOARDERCOLLIE_CONFIG['hosts'][0]['name']
    HOST_MEMBERS = host_reference + ".members.tsv"
    HOST_CONTIG_MAP_SOURCE = "--members " + HOST_MEMBERS
else:
    # A host FASTA not made by BoarderCollie: one member, named after it.
    HOST_NAME = os.path.splitext(os.path.basename(host_reference))[0]
    HOST_MEMBERS = None
    HOST_CONTIG_MAP_SOURCE = ("--reference " + host_reference
                              + " --single-name " + HOST_NAME)


def host_bam(sample):
    """in : sample name
       out: the host BAM host_filter (or host_filter_se) wrote for it"""
    sub = "host_se" if is_single_end(sample) else "host"
    return "output/qc/host_filter/%s/%s.bam" % (sub, sample)


def host_filter_log(sample):
    """in : sample name
       out: its host_filter log, which holds bowtie2's alignment summary"""
    suffix = ".se.log" if is_single_end(sample) else ".log"
    return "output/logs/qc/host_filter/" + sample + suffix


localrules: host_contig_map, merge_host_genome_counts


rule host_contig_map:
    """
    Contig -> member genome, from members.tsv and each member's cleaned
    FASTA in BoarderCollie's library. Reads only FASTA headers.
    """
    input:
        # The record host_bowtie2_build writes once the reference and its
        # members.tsv are in place.
        done=rules.host_bowtie2_build.output,
        reference=ancient(host_reference)
    output:
        HOST_GENOME_COUNTS_DIR + "/host_contig_map.tsv"
    params:
        script=HOST_GENOME_COUNTS_SCRIPT,
        source=HOST_CONTIG_MAP_SOURCE
    conda:
        "../env/qc.yaml"
    log:
        "output/logs/qc/host_genome_counts/host_contig_map.log"
    shell:
        """
        python {params.script} contig-map {params.source} \
          --out {output}.partial 2> {log}
        mv {output}.partial {output}
        """


rule host_genome_counts:
    """
    One sample's host reads, split by member genome. Every record of the
    host BAM is read (samtools view -h, unsorted bowtie2 order).
    """
    input:
        bam=lambda wildcards: host_bam(wildcards.sample),
        log=lambda wildcards: host_filter_log(wildcards.sample),
        contig_map=rules.host_contig_map.output
    output:
        HOST_GENOME_COUNTS_DIR + "/{sample}.tsv"
    params:
        script=HOST_GENOME_COUNTS_SCRIPT,
        layout=lambda wildcards: "SE" if is_single_end(wildcards.sample) else "PE",
        min_mapq=HOST_GENOME_COUNTS_MIN_MAPQ
    conda:
        "../env/qc.yaml"
    threads: 2
    retries:
        config['retries'].get('host_genome_counts', 2)
    resources:
        mem_mb=mem_escalate('host_genome_counts', base_default=4000),
        runtime=runtime_escalate('host_genome_counts', base_default=240)
    log:
        "output/logs/qc/host_genome_counts/{sample}.log"
    benchmark:
        "output/benchmarks/qc/host_genome_counts/{sample}_benchmark.txt"
    shell:
        """
        set -o pipefail
        samtools view -h -@ 1 {input.bam} \
          | python {params.script} count --sample {wildcards.sample} \
              --layout {params.layout} --contig-map {input.contig_map} \
              --bowtie2-log {input.log} --min-mapq {params.min_mapq} \
              --out {output}.partial 2> {log}
        mv {output}.partial {output}
        """


rule merge_host_genome_counts:
    input:
        counts=expand(HOST_GENOME_COUNTS_DIR + "/{sample}.tsv", sample=samples),
        contig_map=rules.host_contig_map.output
    output:
        table="output/qc/host_filter/host_genome_counts.tsv",
        params="output/qc/host_filter/host_genome_counts_params.txt"
    params:
        script=HOST_GENOME_COUNTS_SCRIPT,
        host=HOST_NAME,
        reference=host_reference,
        members=("--members " + HOST_MEMBERS) if HOST_MEMBERS else "",
        min_mapq=HOST_GENOME_COUNTS_MIN_MAPQ
    conda:
        "../env/qc.yaml"
    log:
        "output/logs/qc/host_genome_counts/merge.log"
    shell:
        """
        python {params.script} merge --out {output.table} \
          --params {output.params} --host {params.host} \
          --reference {params.reference} {params.members} \
          --min-mapq {params.min_mapq} {input.counts} 2> {log}
        """
