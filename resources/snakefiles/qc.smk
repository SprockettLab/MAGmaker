from os.path import splitext, basename

host_base = join(config['host_filter']['db_dir'],
                 splitext(basename(config['host_filter']['genome']))[0])

trimmer = config['trimmer']


# Single-end trimmer output lives in its own subdirectory. The per-read
# FASTQ names could not collide (SE versus R1/R2), but the fastp JSON and
# HTML reports carry no read identifier, so a single-end rule writing beside
# the paired one would make two rules claim the same path and Snakemake
# would refuse the DAG as ambiguous. Keeping paired paths exactly where they
# were also means existing runs are not invalidated by this change.
SE_SUBDIR = "se/"


def trimmed_read(sample, seqrun, read):
    """in : sample, sequencing run, read id
       out: path to that run's trimmed FASTQ"""
    return "output/qc/{t}/{d}{s}.{u}.{r}.fastq.gz".format(
        t=trimmer,
        d=SE_SUBDIR if read == SINGLE_END_READ else "",
        s=sample, u=seqrun, r=read
    )


def fastp_report(sample, seqrun):
    """in : sample and sequencing run
       out: path to that run's fastp JSON

    Single-end reports live in the se/ subdirectory alongside the reads,
    so anything reading them has to ask rather than assume. Hardcoding the
    paired-end path made prototype_selection fail on every single-end
    sample with a FileNotFoundError."""
    return "output/qc/fastp/{d}{s}.{u}.fastp.json".format(
        d=SE_SUBDIR if is_single_end(sample) else "", s=sample, u=seqrun)


def trimmer_output(wildcards):
    """Trimmed reads path for the configured trimmer."""
    return trimmed_read(wildcards.sample, wildcards.seqrun, wildcards.read)


def merged_reads(sample, read):
    """Trimmed reads path for host_filter.

    A sample sequenced once needs no concatenation, so its trimmed FASTQ is
    used directly and merge_seqruns never runs for it -- nothing is copied or
    linked. Only multi-run samples pay the cat."""
    runs = seqruns_for(sample)
    if len(runs) == 1:
        return trimmed_read(sample, runs[0], read)
    return "output/qc/merge_seqruns/{s}.combined.{r}.fastq.gz".format(
        s=sample, r=read
    )


if trimmer == 'cutadapt' and any(is_single_end(s) for s in samples):
    _fail(
        "trimmer: cutadapt does not support the single-end samples in this\n"
        "metadata file. The cutadapt parameters are paired-end specific\n"
        "(-U trims the reverse read), so a single-end run would fail or\n"
        "trim the wrong thing.\n\n"
        "Set  trimmer: fastp  to process single-end data."
    )


def trimmer_qc_logs(metadata_table):
    """Trimmer-specific QC files collected by MultiQC."""
    if trimmer == 'fastp':
        return [fastp_report(row.Index[0], row.Index[1])
                for row in metadata_table.itertuples()]
    else:
        return expand(
            "output/logs/qc/cutadapt/{u.Index[0]}.{u.Index[1]}.txt",
            u=metadata_table.itertuples()
        )


rule fastqc_pre_trim:
    input:
        lambda wildcards: get_read(wildcards.sample,
                                   wildcards.seqrun,
                                   wildcards.read)
    output:
        html="output/qc/fastqc_pre_trim/{sample}.{seqrun}.{read}.html",
        zip="output/qc/fastqc_pre_trim/{sample}.{seqrun}.{read}_fastqc.zip"
    benchmark:
        "output/benchmarks/qc/fastqc_pre_trim/{sample}.{seqrun}.{read}_benchmark.txt"
    threads:
        config['threads']['fastqc']
    conda:
        "../env/fastqc.yaml"
    log:
        "output/logs/qc/fastqc_pre_trim/{sample}.{seqrun}.{read}.log"
    shell:
        """
        OUTDIR=$(dirname {output.html})
        mkdir -p $OUTDIR
        fastqc {input} --outdir $OUTDIR --threads {threads} 2> {log}
        # fastqc names its output after the input minus a recognised
        # extension, so the stem has to be derived the same way rather
        # than by stripping one hardcoded suffix. basename only strips a
        # suffix that matches: given a plain .fastq it left .fastq.gz
        # unstripped, the mv below looked for a file that never existed,
        # and the rule failed after fastqc had already succeeded.
        STEM=$(basename {input})
        STEM=${{STEM%.gz}}
        STEM=${{STEM%.fastq}}
        STEM=${{STEM%.fq}}
        mv $OUTDIR/${{STEM}}_fastqc.html {output.html}
        SRC_ZIP="$OUTDIR/${{STEM}}_fastqc.zip"
        [ "$SRC_ZIP" = "{output.zip}" ] || mv "$SRC_ZIP" "{output.zip}"
        """


rule fastp_pe:
    input:
        R1=lambda wildcards: get_read(wildcards.sample, wildcards.seqrun, 'R1'),
        R2=lambda wildcards: get_read(wildcards.sample, wildcards.seqrun, 'R2')
    output:
        R1=temp("output/qc/fastp/{sample}.{seqrun}.R1.fastq.gz"),
        R2=temp("output/qc/fastp/{sample}.{seqrun}.R2.fastq.gz"),
        json="output/qc/fastp/{sample}.{seqrun}.fastp.json",
        html="output/qc/fastp/{sample}.{seqrun}.fastp.html"
    params:
        extra=config['params']['fastp']['extra']
    threads:
        config['threads']['fastp']
    conda:
        "../env/fastp.yaml"
    log:
        "output/logs/qc/fastp/{sample}.{seqrun}.log"
    benchmark:
        "output/benchmarks/qc/fastp/{sample}.{seqrun}_benchmark.txt"
    shell:
        """
        fastp \
            --in1 {input.R1} --in2 {input.R2} \
            --out1 {output.R1} --out2 {output.R2} \
            --json {output.json} --html {output.html} \
            --thread {threads} \
            {params.extra} \
            2> {log}
        """


rule fastp_se:
    input:
        SE=lambda wildcards: get_read(wildcards.sample, wildcards.seqrun,
                                      SINGLE_END_READ)
    output:
        SE=temp("output/qc/fastp/se/{sample}.{seqrun}.SE.fastq.gz"),
        json="output/qc/fastp/se/{sample}.{seqrun}.fastp.json",
        html="output/qc/fastp/se/{sample}.{seqrun}.fastp.html"
    params:
        # params.fastp.extra is written for paired input, where
        # --detect_adapter_for_pe is needed because fastp infers adapters
        # from read overlap rather than detecting the sequence. For
        # single-end input fastp detects adapter sequence by default and
        # the flag is meaningless, so it is stripped while the adapter
        # panel and homopolymer trimming are kept. A single-end specific
        # extra can be set in config as params.fastp.extra_se.
        extra=config['params']['fastp'].get(
            'extra_se',
            config['params']['fastp']['extra'].replace(
                '--detect_adapter_for_pe', '').strip()
        )
    threads:
        config['threads']['fastp']
    conda:
        "../env/fastp.yaml"
    log:
        "output/logs/qc/fastp/se/{sample}.{seqrun}.log"
    benchmark:
        "output/benchmarks/qc/fastp/se/{sample}.{seqrun}_benchmark.txt"
    shell:
        """
        fastp \
            --in1 {input.SE} \
            --out1 {output.SE} \
            --json {output.json} --html {output.html} \
            --thread {threads} \
            {params.extra} \
            2> {log}
        """


rule cutadapt_pe:
    input:
        R1=lambda wildcards: get_read(wildcards.sample, wildcards.seqrun, 'R1'),
        R2=lambda wildcards: get_read(wildcards.sample, wildcards.seqrun, 'R2')
    output:
        fastq1=temp("output/qc/cutadapt/{sample}.{seqrun}.R1.fastq.gz"),
        fastq2=temp("output/qc/cutadapt/{sample}.{seqrun}.R2.fastq.gz"),
        qc="output/logs/qc/cutadapt/{sample}.{seqrun}.txt"
    params:
        adapters=config["params"]["cutadapt"]['adapter'],
        extra=config["params"]["cutadapt"]['other']
    benchmark:
        "output/benchmarks/qc/cutadapt/{sample}.{seqrun}_benchmark.txt"
    log:
        "output/logs/qc/cutadapt/{sample}.{seqrun}.log"
    threads:
        config['threads']['cutadapt_pe']
    conda:
        "../env/cutadapt.yaml"
    shell:
        """
        cutadapt \
            {params.adapters} \
            {params.extra} \
            -j {threads} \
            -o {output.fastq1} \
            -p {output.fastq2} \
            {input.R1} {input.R2} \
            > {output.qc} 2> {log}
        """


rule fastqc_post_trim:
    input:
        trimmer_output
    output:
        html="output/qc/fastqc_post_trim/{sample}.{seqrun}.{read}.html",
        zip="output/qc/fastqc_post_trim/{sample}.{seqrun}.{read}_fastqc.zip"
    benchmark:
        "output/benchmarks/qc/fastqc_post_trim/{sample}.{seqrun}.{read}_benchmark.txt"
    threads:
        config['threads']['fastqc']
    conda:
        "../env/fastqc.yaml"
    log:
        "output/logs/qc/fastqc_post_trim/{sample}.{seqrun}.{read}.log"
    shell:
        """
        OUTDIR=$(dirname {output.html})
        mkdir -p $OUTDIR
        fastqc {input} --outdir $OUTDIR --threads {threads} 2> {log}
        # fastqc names its output after the input minus a recognised
        # extension, so the stem has to be derived the same way rather
        # than by stripping one hardcoded suffix. basename only strips a
        # suffix that matches: given a plain .fastq it left .fastq.gz
        # unstripped, the mv below looked for a file that never existed,
        # and the rule failed after fastqc had already succeeded.
        STEM=$(basename {input})
        STEM=${{STEM%.gz}}
        STEM=${{STEM%.fastq}}
        STEM=${{STEM%.fq}}
        mv $OUTDIR/${{STEM}}_fastqc.html {output.html}
        SRC_ZIP="$OUTDIR/${{STEM}}_fastqc.zip"
        [ "$SRC_ZIP" = "{output.zip}" ] || mv "$SRC_ZIP" "{output.zip}"
        """


rule merge_seqruns:
    input:
        lambda wildcards: [
            trimmed_read(wildcards.sample, seqrun, wildcards.read)
            for seqrun in seqruns_for(wildcards.sample)
        ]
    output:
        temp("output/qc/merge_seqruns/{sample}.combined.{read}.fastq.gz")
    benchmark:
        "output/benchmarks/qc/merge_seqruns/{sample}.combined.{read}_benchmark.txt"
    log:
        "output/logs/qc/merge_seqruns/{sample}.combined.{read}.log"
    threads: 1
    shell:
        "cat {input} > {output[0]} 2> {log}"


rule host_bowtie2_build:
    input:
        reference=config['host_filter']['genome']
    output:
        # bowtie2-build writes small-index (.bt2) files for references
        # under ~4 Gbp and large-index (.bt2l) files above it. Declaring
        # fixed .bt2 outputs breaks on large hosts (e.g. multi-genome
        # vervet references), where only .bt2l files appear. Track a
        # sentinel instead; `bowtie2 -x` auto-detects the index layout at
        # align time.
        touch(host_base + ".bowtie2_build.done")
    log:
        "output/logs/qc/host_bowtie2_build/host_bowtie2_build.log"
    benchmark:
        "output/benchmarks/qc/host_bowtie2_build/host_bowtie2_build_benchmark.txt"
    conda:
        "../env/qc.yaml"
    params:
        extra="",
        indexbase=host_base
    threads:
        config['threads']['host_filter']
    shell:
        # A complete index next to the reference is reused rather than
        # rebuilt. Shared host indexes commonly predate the sentinel, and
        # rebuilding a human genome costs hours -- silently, since the only
        # symptom is host_bowtie2_build appearing in the job list.
        #
        # All six parts must be present before the build is skipped: a
        # partial index (an interrupted build, an incomplete copy) would
        # otherwise be accepted and fail later inside bowtie2, which is a
        # much harder failure to trace back to its cause. Small (.bt2) and
        # large (.bt2l) layouts are checked separately, since a reference
        # produces one or the other, never a mixture.
        """
        IDX="{params.indexbase}"
        FOUND=0
        for suffix in bt2 bt2l; do
            COMPLETE=1
            for part in 1 2 3 4 rev.1 rev.2; do
                [ -f "$IDX.$part.$suffix" ] || COMPLETE=0
            done
            if [ "$COMPLETE" -eq 1 ]; then FOUND=1; fi
        done

        if [ "$FOUND" -eq 1 ]; then
            echo "Complete bowtie2 index already present at $IDX" > {log}
            echo "Skipping bowtie2-build; touching sentinel instead." >> {log}
        else
            bowtie2-build --threads {threads} {params.extra} \
            {input.reference} "$IDX" 2> {log} 1>&2
        fi
        """


rule host_filter:
    input:
        fastq1=lambda wildcards: merged_reads(wildcards.sample, 'R1'),
        fastq2=lambda wildcards: merged_reads(wildcards.sample, 'R2'),
        db=rules.host_bowtie2_build.output
    output:
        nonhost_R1="output/qc/host_filter/nonhost/{sample}.R1.fastq.gz",
        nonhost_R2="output/qc/host_filter/nonhost/{sample}.R2.fastq.gz",
        host="output/qc/host_filter/host/{sample}.bam",
    params:
        ref=host_base
    conda:
        "../env/qc.yaml"
    threads:
        config['threads']['host_filter']
    log:
        "output/logs/qc/host_filter/{sample}.log"
    benchmark:
        "output/benchmarks/qc/host_filter/{sample}_benchmark.txt"
    shell:
        """
        bowtie2 -p {threads} -x {params.ref} \
          -1 {input.fastq1} -2 {input.fastq2} \
          --un-conc-gz {wildcards.sample}_nonhost \
          --no-unal \
          2> {log} | samtools view -bS - > {output.host}

        mv {wildcards.sample}_nonhost.1 output/qc/host_filter/nonhost/{wildcards.sample}.R1.fastq.gz
        mv {wildcards.sample}_nonhost.2 output/qc/host_filter/nonhost/{wildcards.sample}.R2.fastq.gz
        """


rule host_filter_se:
    """
    Host filtering for single-end samples.

    Separate from host_filter rather than a branch inside it because the
    two produce different numbers of output files, which a single rule
    cannot declare. The host BAM goes to its own directory so the two rules
    never match the same path.
    """
    input:
        fastq=lambda wildcards: merged_reads(wildcards.sample,
                                             SINGLE_END_READ),
        db=rules.host_bowtie2_build.output
    output:
        nonhost_SE="output/qc/host_filter/nonhost/{sample}.SE.fastq.gz",
        host="output/qc/host_filter/host_se/{sample}.bam",
    params:
        ref=host_base
    conda:
        "../env/qc.yaml"
    threads:
        config['threads']['host_filter']
    log:
        "output/logs/qc/host_filter/{sample}.se.log"
    benchmark:
        "output/benchmarks/qc/host_filter/{sample}.se_benchmark.txt"
    shell:
        """
        bowtie2 -p {threads} -x {params.ref} \
          -U {input.fastq} \
          --un-gz {output.nonhost_SE} \
          --no-unal \
          2> {log} | samtools view -bS - > {output.host}
        """


rule fastqc_post_host:
    input:
        "output/qc/host_filter/nonhost/{sample}.{read}.fastq.gz"
    output:
        html="output/qc/fastqc_post_host/{sample}.{read}.html",
        zip="output/qc/fastqc_post_host/{sample}.{read}_fastqc.zip"
    benchmark:
        "output/benchmarks/qc/fastqc_post_host/{sample}.{read}_benchmark.txt"
    threads:
        config['threads']['fastqc']
    conda:
        "../env/fastqc.yaml"
    log:
        "output/logs/qc/fastqc_post_host/{sample}.{read}.log"
    shell:
        """
        OUTDIR=$(dirname {output.html})
        mkdir -p $OUTDIR
        fastqc {input} --outdir $OUTDIR --threads {threads} 2> {log}
        # fastqc names its output after the input minus a recognised
        # extension, so the stem has to be derived the same way rather
        # than by stripping one hardcoded suffix. basename only strips a
        # suffix that matches: given a plain .fastq it left .fastq.gz
        # unstripped, the mv below looked for a file that never existed,
        # and the rule failed after fastqc had already succeeded.
        STEM=$(basename {input})
        STEM=${{STEM%.gz}}
        STEM=${{STEM%.fastq}}
        STEM=${{STEM%.fq}}
        mv $OUTDIR/${{STEM}}_fastqc.html {output.html}
        SRC_ZIP="$OUTDIR/${{STEM}}_fastqc.zip"
        [ "$SRC_ZIP" = "{output.zip}" ] || mv "$SRC_ZIP" "{output.zip}"
        """


def fastqc_run_targets(stage):
    """in : 'pre_trim' or 'post_trim'
       out: per-run FastQC targets, one per read the sample actually has

    Iterating reads_for(sample) rather than the global `reads` is what keeps
    an R2 report from being requested for a single-end sample."""
    return [
        "output/qc/fastqc_{stage}/{s}.{u}.{r}.html".format(
            stage=stage, s=row.Index[0], u=row.Index[1], r=read)
        for row in metadata_table.itertuples()
        for read in reads_for(row.Index[0])
    ]


def fastqc_host_targets():
    """Per-sample FastQC targets after host filtering."""
    return [
        "output/qc/fastqc_post_host/{s}.{r}.html".format(s=sample, r=read)
        for sample in samples
        for read in reads_for(sample)
    ]


def host_filter_logs():
    """Host filter logs; the two rules write to different names."""
    return [
        "output/logs/qc/host_filter/{s}.se.log".format(s=sample)
        if is_single_end(sample)
        else "output/logs/qc/host_filter/{s}.log".format(s=sample)
        for sample in samples
    ]


rule multiqc:
    input:
        lambda wildcards: fastqc_run_targets('pre_trim'),
        lambda wildcards: trimmer_qc_logs(metadata_table),
        lambda wildcards: fastqc_run_targets('post_trim'),
        lambda wildcards: fastqc_host_targets(),
        lambda wildcards: host_filter_logs()
    output:
        "output/qc/multiqc/multiqc.html"
    params:
        "--dirs " + config['params']['multiqc']
    log:
        "output/logs/qc/multiqc/multiqc.log"
    benchmark:
        "output/benchmarks/qc/multiqc/multiqc_benchmark.txt"
    resources:
        mem_mb=config['mem_mb']['multiqc']
    wrapper:
        "v3.1.0/bio/multiqc"
