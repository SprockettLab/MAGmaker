# =============================================================================
# long_reads.smk -- Oxford Nanopore long reads (Platform ont rows)
#
# A sample's long reads are optional extra input. Everything here is planned
# only for long_read_samples, so a sheet without ont rows plans exactly the
# same jobs as before. Steps, in order:
#
#   host_minimap2_index   the run's host reference, indexed for map-ont
#                         (BoarderCollie's shared index when it is on)
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


# The host's minimap2 map-ont index. With BoarderCollie on (the default) it
# is BoarderCollie's shared index, built once per host key in its index_dir
# and linked to host_minimap2_out (host_reference.smk), so no run folder
# builds its own. With BoarderCollie off (host_filter.genome used as it is),
# MAGmaker indexes that FASTA in the run folder.
if USE_BOARDERCOLLIE:
    HOST_MINIMAP2_INDEX = host_minimap2_out
else:
    HOST_MINIMAP2_INDEX = "output/qc/host_filter_ont/index/host.map-ont.mmi"

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
            HOST_MINIMAP2_INDEX
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
        index=HOST_MINIMAP2_INDEX
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


# ---------------------------------------------------------------------------
# Short-read polishing: Polypolish --careful, then pypolca --careful, with the
# sample's QC'd, host-filtered short reads (nonhost_reads). Nothing else.
# Careful mode ignores reads that align equally well to several places and
# changes a base only with strong agreement, so related strains in one sample
# are less likely to pull a contig toward the wrong one.
# ---------------------------------------------------------------------------

_single_end_lr = [s for s in long_read_samples if is_single_end(s)]
if _single_end_lr:
    _fail("long-read sample(s) with single-end short reads: %s\n"
          "Polypolish and pypolca need paired short reads." % ", ".join(_single_end_lr))


rule polish_long_read_assembly:
    """
    bwa mem -a for each mate separately, polypolish filter, polypolish polish
    --careful, pypolca run --careful. The draft is kept. Polishing must keep
    every contig name (checked by polish_changes.py, which also writes the
    bases each stage changed per contig).
    """
    input:
        draft="output/assemble/{lr_assembler}/{sample}.draft.fasta",
        reads=lambda wildcards: nonhost_reads(wildcards.sample)
    output:
        polished="output/assemble/{lr_assembler}/{sample}.polished.fasta",
        changes="output/assemble/{lr_assembler}/{sample}.polish_changes.tsv"
    params:
        workdir="output/assemble/{lr_assembler}/{sample}_polish",
        script=os.path.join(LONG_READ_SCRIPTS, "polish_changes.py")
    wildcard_constraints:
        lr_assembler="flye|myloasm"
    conda:
        "../env/polish.yaml"
    threads:
        config['threads'].get('polish_long_read_assembly', 32)
    retries:
        config['retries'].get('polish_long_read_assembly', 1)
    resources:
        mem_mb=mem_escalate('polish_long_read_assembly', base_default=64000),
        runtime=runtime_escalate('polish_long_read_assembly', base_default=2880)
    log:
        "output/logs/assemble/polish/{lr_assembler}/{sample}.log"
    benchmark:
        "output/benchmarks/assemble/polish/{lr_assembler}/{sample}_benchmark.txt"
    shell:
        """
        W={params.workdir}
        rm -rf "$W" && mkdir -p "$W"
        {{
        echo "polypolish $(polypolish --version) pypolca $(pypolca --version 2>&1 | tail -1) bwa $(bwa 2>&1 | grep -i '^version' | cut -d' ' -f2)"
        bwa index -p "$W/draft" {input.draft}
        bwa mem -t {threads} -a "$W/draft" {input.reads[0]} > "$W/align_R1.sam"
        bwa mem -t {threads} -a "$W/draft" {input.reads[1]} > "$W/align_R2.sam"
        polypolish filter --in1 "$W/align_R1.sam" --in2 "$W/align_R2.sam" \
            --out1 "$W/filtered_R1.sam" --out2 "$W/filtered_R2.sam"
        polypolish polish --careful {input.draft} "$W/filtered_R1.sam" "$W/filtered_R2.sam" \
            > "$W/polypolish.fasta"
        rm -f "$W"/align_R?.sam "$W"/filtered_R?.sam "$W"/draft.*
        pypolca run -a "$W/polypolish.fasta" -1 {input.reads[0]} -2 {input.reads[1]} \
            -t {threads} -o "$W/pypolca" --careful --force
        # Keep only the contig name in each header; the polishers append text.
        awk '/^>/ {{print $1; next}} {{print}}' "$W/pypolca/pypolca_corrected.fasta" > "$W/final.fasta"
        awk '/^>/ {{print $1; next}} {{print}}' "$W/polypolish.fasta" > "$W/polypolish.clean.fasta"
        python {params.script} --draft {input.draft} --polypolish "$W/polypolish.clean.fasta" \
            --final "$W/final.fasta" --out {output.changes}.partial
        }} > {log} 2>&1
        mv "$W/final.fasta" {output.polished}
        mv {output.changes}.partial {output.changes}
        rm -rf "$W"
        """


# ---------------------------------------------------------------------------
# DNA Control Strand removal, right after polishing and before anything else
# reads the assembly: lambda is a phage, so the viral track would report it.
# The result is output/assemble/<assembler>/<sample>.contigs.fasta, the same
# pattern as MEGAHIT's, which binning and the viral track read.
# ---------------------------------------------------------------------------

rule remove_control_strand:
    """
    Contigs at least params.long_reads.control_strand_min_aligned_fraction
    (0.5) aligned to the ligation kit's DNA Control Strand (lambda,
    NC_001416.1; params.long_reads.control_strand_fasta) are removed and kept
    aside with a table. Checks kept + removed = contigs in.
    """
    input:
        polished="output/assemble/{lr_assembler}/{sample}.polished.fasta",
        control=long_read_param("control_strand_fasta",
            "/isilon/datalake/sprockett_lab/original/WF00SprockettLab/dbs/lambda/NC_001416.1.fasta")
    output:
        contigs="output/assemble/{lr_assembler}/{sample}.contigs.fasta",
        removed="output/assemble/{lr_assembler}/{sample}.control_strand.fasta",
        table="output/assemble/{lr_assembler}/{sample}.control_strand.tsv"
    params:
        min_fraction=long_read_param("control_strand_min_aligned_fraction", 0.5),
        script=os.path.join(LONG_READ_SCRIPTS, "control_strand.py")
    wildcard_constraints:
        lr_assembler="flye|myloasm"
    conda:
        "../env/long_read_qc.yaml"
    threads: 4
    resources:
        mem_mb=16000,
        runtime=240
    log:
        "output/logs/assemble/control_strand/{lr_assembler}/{sample}.log"
    shell:
        """
        minimap2 -c -t {threads} {input.control} {input.polished} 2> {log} > {output.table}.paf
        python {params.script} --contigs {input.polished} --paf {output.table}.paf \
            --min-fraction {params.min_fraction} --kept {output.contigs}.partial \
            --removed {output.removed} --table {output.table} 2>> {log}
        rm -f {output.table}.paf
        mv {output.contigs}.partial {output.contigs}
        """
