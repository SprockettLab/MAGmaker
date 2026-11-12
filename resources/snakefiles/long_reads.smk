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


# ---------------------------------------------------------------------------
# Assembly: every assembler in long_read_assemblers, for every long-read sample.
# Each writes output/assemble/<assembler>/<sample>.draft.fasta (contigs renamed
# <sample>__<assembler>__<original>) and <sample>.contig_info.tsv (length,
# circular, depth), kept beside the polished contigs. The assembler's own
# folder is kept under <sample>_<assembler>/ while it runs, so a retry resumes.
# ---------------------------------------------------------------------------

LONG_READ_ASSEMBLERS = config.get("long_read_assemblers", []) or []
_unknown_lr = [a for a in LONG_READ_ASSEMBLERS if a not in ("flye", "myloasm")]
if _unknown_lr:
    _fail("long_read_assemblers names unknown assembler(s): %s\n"
          "Known: flye, myloasm." % ", ".join(_unknown_lr))
if long_read_samples and not LONG_READ_ASSEMBLERS:
    _fail("the metadata has long reads (Platform ont) but long_read_assemblers\n"
          "is empty; list flye and/or myloasm.")


rule flye:
    """
    metaFlye: flye --meta --nano-hq, defaults otherwise. A retry resumes the
    same run (--resume) when Flye got far enough to write params.json.
    """
    input:
        reads="output/qc/host_filter_ont/nonhost/{sample}.ont.fastq.gz"
    output:
        draft="output/assemble/flye/{sample}.draft.fasta",
        info="output/assemble/flye/{sample}.contig_info.tsv",
        report="output/assemble/flye/{sample}_flye/assembly_info.txt"
    params:
        workdir="output/assemble/flye/{sample}_flye/run",
        script=os.path.join(LONG_READ_SCRIPTS, "long_read_contigs.py"),
        extra=long_read_param("flye_extra", "")
    conda:
        "../env/flye.yaml"
    threads:
        config['threads'].get('flye', 64)
    retries:
        config['retries'].get('flye', 1)
    resources:
        mem_mb=mem_escalate('flye', base_default=500000),
        runtime=runtime_escalate('flye', base_default=10080)
    log:
        "output/logs/assemble/flye/{sample}.log"
    benchmark:
        "output/benchmarks/assemble/flye/{sample}_benchmark.txt"
    shell:
        """
        if [ -s {params.workdir}/params.json ] && [ ! -s {params.workdir}/assembly.fasta ]; then
            flye --nano-hq {input.reads} --meta --threads {threads} --out-dir {params.workdir} \
                 --resume {params.extra} > {log} 2>&1
        elif [ ! -s {params.workdir}/assembly.fasta ]; then
            rm -rf {params.workdir}
            flye --nano-hq {input.reads} --meta --threads {threads} --out-dir {params.workdir} \
                 {params.extra} > {log} 2>&1
        fi
        python {params.script} --sample {wildcards.sample} --assembler flye \
            --assembly {params.workdir}/assembly.fasta --flye-info {params.workdir}/assembly_info.txt \
            --out-fasta {output.draft}.partial --out-info {output.info}.partial
        cp {params.workdir}/assembly_info.txt {params.workdir}/assembly_graph.gfa \
           {params.workdir}/flye.log $(dirname {output.report})/
        mv {output.draft}.partial {output.draft}
        mv {output.info}.partial {output.info}
        rm -rf {params.workdir}
        """


rule myloasm:
    """
    myloasm with default settings; its primary assembly (assembly_primary.fa).
    Circular, length and depth come from its contig names; only circular-yes
    counts as circular downstream.
    """
    input:
        reads="output/qc/host_filter_ont/nonhost/{sample}.ont.fastq.gz"
    output:
        draft="output/assemble/myloasm/{sample}.draft.fasta",
        info="output/assemble/myloasm/{sample}.contig_info.tsv"
    params:
        workdir="output/assemble/myloasm/{sample}_myloasm",
        script=os.path.join(LONG_READ_SCRIPTS, "long_read_contigs.py"),
        extra=long_read_param("myloasm_extra", "")
    conda:
        "../env/myloasm.yaml"
    threads:
        config['threads'].get('myloasm', 64)
    retries:
        config['retries'].get('myloasm', 1)
    resources:
        mem_mb=mem_escalate('myloasm', base_default=500000),
        runtime=runtime_escalate('myloasm', base_default=10080)
    log:
        "output/logs/assemble/myloasm/{sample}.log"
    benchmark:
        "output/benchmarks/assemble/myloasm/{sample}_benchmark.txt"
    shell:
        """
        if [ ! -s {params.workdir}/assembly_primary.fa ]; then
            rm -rf {params.workdir}
            myloasm {input.reads} -o {params.workdir} -t {threads} {params.extra} > {log} 2>&1
        fi
        python {params.script} --sample {wildcards.sample} --assembler myloasm \
            --assembly {params.workdir}/assembly_primary.fa \
            --out-fasta {output.draft}.partial --out-info {output.info}.partial
        mv {output.draft}.partial {output.draft}
        mv {output.info}.partial {output.info}
        """
