# =============================================================================
# long_reads.smk -- Oxford Nanopore long reads (Platform ont rows)
#
# A sample's long reads are optional extra input. Everything here is planned
# only for long_read_samples, so a sheet without ont rows plans exactly the
# same jobs as before. Steps, in order:
#
#   host_minimap2_index   the run's host reference, indexed for map-ont
#   host_filter_ont       host reads removed (union of aligned bases >= 50%)
#
# The host reference is the one the run's short reads are filtered against
# (host_reference, from host_reference.smk), so both read types see one host.
# =============================================================================

import os

LONG_READ_SCRIPTS = os.path.normpath(os.path.join(workflow.current_basedir, "..", "scripts"))
LONG_READ_PARAMS = config.get("params", {}).get("long_reads", {}) or {}


def long_read_param(key, default):
    return LONG_READ_PARAMS.get(key, default)


rule host_minimap2_index:
    """
    The run's host reference indexed for minimap2 map-ont. Built as ONE index
    part (-I larger than the reference): minimap2 splits a reference larger
    than -I (default 8 Gbp; Chlorocebus_spp is about 8.8 Gbp) into parts and,
    without --split-prefix, aligns each part independently, so a read could
    receive one "primary" alignment per part and the host rule would count
    them as if they were one alignment set.
    """
    input:
        ref=host_reference
    output:
        "output/qc/host_filter_ont/index/host.map-ont.mmi"
    params:
        split=long_read_param("host_index_max_bases", "16G")
    conda:
        "../env/long_read_qc.yaml"
    threads:
        config['threads'].get('host_minimap2_index', 8)
    retries:
        config['retries'].get('host_minimap2_index', 1)
    resources:
        mem_mb=mem_escalate('host_minimap2_index', base_default=64000),
        runtime=runtime_escalate('host_minimap2_index', base_default=240)
    log:
        "output/logs/qc/host_filter_ont/host_minimap2_index.log"
    benchmark:
        "output/benchmarks/qc/host_filter_ont/host_minimap2_index_benchmark.txt"
    shell:
        """
        minimap2 -x map-ont -I {params.split} -t {threads} -d {output}.partial {input.ref} 2> {log}
        mv {output}.partial {output}
        """


rule host_filter_ont:
    """
    Long-read host removal. A read is host when the union of the read bases
    its primary and supplementary alignments cover (--secondary=no) is at
    least params.long_reads.host_min_aligned_fraction (0.5) of its length.
    Reads with some host alignment below that are kept and counted.

    Outputs: the nonhost reads; the host alignments as a BAM (as host_filter
    keeps the short reads'), and the host read IDs; counts.tsv with every
    parameter. The nonhost read count is checked against counts.tsv. Several
    ont rows of one sample are read as one stream.
    """
    input:
        reads=lambda wildcards: long_reads_for(wildcards.sample),
        index=rules.host_minimap2_index.output
    output:
        nonhost="output/qc/host_filter_ont/nonhost/{sample}.ont.fastq.gz",
        bam="output/qc/host_filter_ont/host/{sample}.bam",
        ids="output/qc/host_filter_ont/host/{sample}.host_ids.txt.gz",
        counts="output/qc/host_filter_ont/{sample}.counts.tsv"
    params:
        threshold=long_read_param("host_min_aligned_fraction", 0.5),
        reference=host_reference,
        script=os.path.join(LONG_READ_SCRIPTS, "ont_host_filter.py")
    conda:
        "../env/long_read_qc.yaml"
    threads:
        config['threads'].get('host_filter_ont', 32)
    retries:
        config['retries'].get('host_filter_ont', 2)
    resources:
        mem_mb=mem_escalate('host_filter_ont', base_default=64000),
        runtime=runtime_escalate('host_filter_ont', base_default=1440)
    log:
        "output/logs/qc/host_filter_ont/{sample}.log"
    benchmark:
        "output/benchmarks/qc/host_filter_ont/{sample}_benchmark.txt"
    shell:
        """
        mkdir -p output/qc/host_filter_ont/tmp
        TMP=$(mktemp -d output/qc/host_filter_ont/tmp/{wildcards.sample}.XXXXXX)
        trap 'rm -rf "$TMP"' EXIT
        MM2=$(minimap2 --version)

        cat {input.reads} \
          | minimap2 -a -x map-ont --secondary=no -t {threads} {input.index} - 2> {log} \
          | python {params.script} --threshold {params.threshold} \
              --host-ids "$TMP/host_ids.txt.gz" --counts "$TMP/counts.tsv" \
              --sample {wildcards.sample} --reference {params.reference} \
              --minimap2-version "$MM2" \
          | samtools view -b -o "$TMP/host.bam" -

        cat {input.reads} \
          | seqkit grep -v -f <(gzip -dc "$TMP/host_ids.txt.gz") -j {threads} \
              -o "$TMP/nonhost.fastq.gz" 2>> {log}

        # The kept reads must be exactly the reads counted as kept.
        got=$(seqkit stats -T "$TMP/nonhost.fastq.gz" | awk -F'\\t' 'NR == 2 {{print $4}}')
        want=$(awk -F'\\t' '$1 == "reads_kept" {{print $2}}' "$TMP/counts.tsv")
        if [ "$got" != "$want" ]; then
            echo "CHECK FAILED: $got nonhost reads written, $want counted as kept" >> {log}
            exit 1
        fi

        mv "$TMP/nonhost.fastq.gz" {output.nonhost}
        mv "$TMP/host.bam" {output.bam}
        mv "$TMP/host_ids.txt.gz" {output.ids}
        mv "$TMP/counts.tsv" {output.counts}
        """
