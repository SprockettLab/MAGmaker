from os.path import splitext, basename

host_base = join(config['host_filter']['db_dir'],
                 splitext(basename(config['host_filter']['genome']))[0])

trimmer = config['trimmer']


def trimmer_output(wildcards):
    """Trimmed reads path for the configured trimmer."""
    return "output/qc/{t}/{s}.{u}.{r}.fastq.gz".format(
        t=trimmer, s=wildcards.sample, u=wildcards.unit, r=wildcards.read
    )


def merged_reads(sample, read):
    """Trimmed reads path for host_filter: skip merge_units for single-unit samples."""
    units = list(units_table.loc[sample].index)
    if len(units) == 1:
        return "output/qc/{t}/{s}.{u}.{r}.fastq.gz".format(
            t=trimmer, s=sample, u=units[0], r=read
        )
    return "output/qc/merge_units/{s}.combined.{r}.fastq.gz".format(
        s=sample, r=read
    )


def trimmer_qc_logs(units_table):
    """Trimmer-specific QC files collected by MultiQC."""
    if trimmer == 'fastp':
        return expand(
            "output/qc/fastp/{u.Index[0]}.{u.Index[1]}.fastp.json",
            u=units_table.itertuples()
        )
    else:
        return expand(
            "output/logs/qc/cutadapt/{u.Index[0]}.{u.Index[1]}.txt",
            u=units_table.itertuples()
        )


rule fastqc_pre_trim:
    input:
        lambda wildcards: get_read(wildcards.sample,
                                   wildcards.unit,
                                   wildcards.read)
    output:
        html="output/qc/fastqc_pre_trim/{sample}.{unit}.{read}.html",
        zip="output/qc/fastqc_pre_trim/{sample}.{unit}.{read}_fastqc.zip"
    params: ""
    benchmark:
        "output/benchmarks/qc/fastqc_pre_trim/{sample}.{unit}.{read}_benchmark.txt"
    threads:
        config['threads']['fastqc']
    wrapper:
        "0.72.0/bio/fastqc"


rule fastp_pe:
    input:
        R1=lambda wildcards: get_read(wildcards.sample, wildcards.unit, 'R1'),
        R2=lambda wildcards: get_read(wildcards.sample, wildcards.unit, 'R2')
    output:
        R1=temp("output/qc/fastp/{sample}.{unit}.R1.fastq.gz"),
        R2=temp("output/qc/fastp/{sample}.{unit}.R2.fastq.gz"),
        json="output/qc/fastp/{sample}.{unit}.fastp.json",
        html="output/qc/fastp/{sample}.{unit}.fastp.html"
    params:
        extra=config['params']['fastp']['extra']
    threads:
        config['threads']['fastp']
    conda:
        "../env/fastp.yaml"
    log:
        "output/logs/qc/fastp/{sample}.{unit}.log"
    benchmark:
        "output/benchmarks/qc/fastp/{sample}.{unit}_benchmark.txt"
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


rule cutadapt_pe:
    input:
        lambda wildcards: get_read(wildcards.sample, wildcards.unit, 'R1'),
        lambda wildcards: get_read(wildcards.sample, wildcards.unit, 'R2')
    output:
        fastq1=temp("output/qc/cutadapt/{sample}.{unit}.R1.fastq.gz"),
        fastq2=temp("output/qc/cutadapt/{sample}.{unit}.R2.fastq.gz"),
        qc="output/logs/qc/cutadapt/{sample}.{unit}.txt"
    params:
        "{} {}".format(config["params"]["cutadapt"]['adapter'],
                       config["params"]["cutadapt"]['other'])
    benchmark:
        "output/benchmarks/qc/cutadapt/{sample}.{unit}_benchmark.txt"
    log:
        "output/logs/qc/cutadapt/{sample}.{unit}.log"
    threads:
        config['threads']['cutadapt_pe']
    wrapper:
        "0.17.4/bio/cutadapt/pe"


rule fastqc_post_trim:
    input:
        trimmer_output
    output:
        html="output/qc/fastqc_post_trim/{sample}.{unit}.{read}.html",
        zip="output/qc/fastqc_post_trim/{sample}.{unit}.{read}_fastqc.zip"
    params: ""
    benchmark:
        "output/benchmarks/qc/fastqc_post_trim/{sample}.{unit}.{read}_benchmark.txt"
    threads:
        config['threads']['fastqc']
    wrapper:
        "0.72.0/bio/fastqc"


rule merge_units:
    input:
        lambda wildcards: expand(
            "output/qc/{t}/{s}.{u}.{r}.fastq.gz",
            t=trimmer,
            s=wildcards.sample,
            u=list(units_table.loc[wildcards.sample].index),
            r=wildcards.read
        )
    output:
        temp("output/qc/merge_units/{sample}.combined.{read}.fastq.gz")
    benchmark:
        "output/benchmarks/qc/merge_units/{sample}.combined.{read}_benchmark.txt"
    log:
        "output/logs/qc/merge_units/{sample}.combined.{read}.log"
    threads: 1
    shell:
        "cat {input} > {output[0]} 2> {log}"


rule host_bowtie2_build:
    input:
        reference=config['host_filter']['genome']
    output:
        multiext(host_base,
                 ".1.bt2",
                 ".2.bt2",
                 ".3.bt2",
                 ".4.bt2",
                 ".rev.1.bt2",
                 ".rev.2.bt2")
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
        """
        bowtie2-build --threads {threads} {params.extra} \
        {input.reference} {params.indexbase} 2> {log} 1>&2
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


rule fastqc_post_host:
    input:
        "output/qc/host_filter/nonhost/{sample}.{read}.fastq.gz"
    output:
        html="output/qc/fastqc_post_host/{sample}.{read}.html",
        zip="output/qc/fastqc_post_host/{sample}.{read}_fastqc.zip"
    params: ""
    benchmark:
        "output/benchmarks/qc/fastqc_post_host/{sample}.{read}_benchmark.txt"
    threads:
        config['threads']['fastqc']
    wrapper:
        "0.72.0/bio/fastqc"


rule multiqc:
    input:
        expand("output/qc/fastqc_pre_trim/{units.Index[0]}.{units.Index[1]}.{read}.html",
               units=units_table.itertuples(), read=reads),
        lambda wildcards: trimmer_qc_logs(units_table),
        expand("output/qc/fastqc_post_trim/{units.Index[0]}.{units.Index[1]}.{read}.html",
               units=units_table.itertuples(), read=reads),
        expand("output/qc/fastqc_post_host/{units.Index[0]}.{read}.html",
               units=units_table.itertuples(), read=reads),
        lambda wildcards: expand(rules.host_filter.log, sample=samples)
    output:
        "output/qc/multiqc/multiqc.html"
    params:
        "--dirs " + config['params']['multiqc']
    log:
        "output/logs/qc/multiqc/multiqc.log"
    benchmark:
        "output/benchmarks/qc/multiqc/multiqc_benchmark.txt"
    wrapper:
        "v3.1.0/bio/multiqc"
