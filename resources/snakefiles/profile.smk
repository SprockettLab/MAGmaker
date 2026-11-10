rule taxonomy_kraken:
    """
    Runs Kraken with Bracken to construct taxonomic profiles.
    """
    input:
        reads=lambda wildcards: nonhost_reads(wildcards.sample)
    output:
        report = "output/profile/kraken2/{sample}.report.txt"
    params:
        db = config['params']['kraken2']['db'],
        levels = config['params']['kraken2']['levels'],
        bracken_db = config['params']['kraken2']['bracken-db'],
        # --paired describes the input, so it must not be passed when a
        # single FASTQ is given; kraken2 errors out rather than ignoring it.
        paired = lambda wildcards, input: (
            "--paired" if len(input.reads) == 2 else "")
    conda:
        "../env/profile.yaml"
    threads:
        config['threads']['kraken2']
    retries:
        config['retries'].get('taxonomy_kraken', 2)
    resources:
        # kraken2 reads a ~200 GB index per job, so the cost is dominated by
        # that read rather than by classification. Left uncapped, every
        # sample in an arm competes for the same shared storage at once and
        # jobs are cancelled on wall clock while doing almost no work. This
        # is a counter, not a reservation: set `resources: [kraken_slots=N]`
        # in the profile to allow N concurrent kraken jobs per workflow.
        kraken_slots=1,
        runtime=runtime_escalate('taxonomy_kraken', base_default=480)
    log:
        "output/logs/profile/kraken2/taxonomy_kraken/{sample}.log"
    benchmark:
        "output/benchmarks/profile/kraken2/taxonomy_kraken/{sample}_benchmark.txt"
    shell:
        """
          # get stem file path
          stem={output.report}
          stem=${{stem%.report.txt}}

          # run Kraken to align reads against reference genomes
          kraken2 {input.reads} \
            --db {params.db} \
            {params.paired} \
            --gzip-compressed \
            --only-classified-output \
            --threads {threads} \
            --report {output.report} \
            --output - \
            2> {log}

          # run Bracken to re-estimate abundance at given rank
          if [[ ! -z {params.levels} ]]
          then
            IFS=',' read -r -a levels <<< "{params.levels}"
            for level in "${{levels[@]}}"
            do
              bracken \
                -d {params.bracken_db} \
                -i {output.report} \
                -t 10 \
                -l $(echo $level | head -c 1 | tr a-z A-Z) \
                -o $stem.redist.$level.txt \
                2>> {log} 1>&2
            done
          fi
          """

rule krona:
    input:
        rules.taxonomy_kraken.output.report
    output:
        "output/profile/krona/{sample}.report.html"
    conda:
        "../env/profile.yaml"
    threads:
        1
    log:
        "output/logs/profile/krona/{sample}.log"
    benchmark:
        "output/benchmarks/profile/krona/{sample}_benchmark.txt"
    shell:
        """
        perl resources/scripts/kraken2-translate.pl {input} > {input}.temp
        ktImportText -o {output} {input}.temp
        rm {input}.temp
        """

rule kraken:
    input:
        expand("output/profile/kraken2/{sample}.report.txt",
               sample=samples),
        expand("output/profile/krona/{sample}.report.html",
               sample=samples)

rule download_metaphlan_db:
    output:
        directory(config['params']['metaphlan']['db_path'])
    conda:
        "../env/profile.yaml"
    log:
        "output/logs/profile/download_metaphlan_db/download_metaphlan_db.log"
    benchmark:
        "output/benchmarks/profile/download_metaphlan_db/download_metaphlan_db_benchmark.txt"
    params:
        db_name=config['params']['metaphlan']['db_name']
    shell:
        # The index is named explicitly rather than downloaded as 'latest'.
        # With 'latest' this rule fetched whatever was current on the day it
        # ran, while the metaphlan rule below profiles against db_name from
        # the config, so the downloaded index and the one used could differ
        # with nothing recording that they had.
        """
        metaphlan --install \
            --db_dir {output} \
            -x {params.db_name} \
            2> {log} 1>&2
        """

rule metaphlan:
    """
    Performs taxonomic profiling using MetaPhlAn4.

    Viral Sequence Cluster (VSC) profiling runs in the SAME invocation: after
    the taxonomic profile is written, metaphlan maps the same reads a second
    time against the viral marker sequences of the index ({index}_VSG.fna and
    {index}_VINFO.csv, both shipped with the database) and writes a per-cluster
    breadth-of-coverage table. Columns: M-Group/Cluster, genomeName, len,
    breadth_of_coverage, mapping_reads_count, RPKM, depth_of_coverage_mean,
    depth_of_coverage_median, plus the descriptive fields of VINFO.

    --vsc_breadth is deliberately NOT passed, so MetaPhlAn's own default
    applies: 0.75 for short reads, 0.5 for long reads. Passing it would fix
    one value for both.

    This detects KNOWN viral clusters from reads, which is a different
    question from assembly-based viral discovery (geNomad in the virus track).
    It finds viruses that never assembled, but only ones already in the VSC
    database.
    """
    input:
        reads=lambda wildcards: nonhost_reads(wildcards.sample),
        db_path=rules.download_metaphlan_db.output
    output:
        bt2="output/profile/metaphlan/bowtie2s/{sample}.bowtie2.bz2",
        sam="output/profile/metaphlan/sams/{sample}.sam.bz2",
        profile="output/profile/metaphlan/profiles/{sample}.txt",
        # Viral Sequence Cluster table, written beside the taxonomic
        # profile. Declared only when VSC profiling is on, so a run with it
        # off does not promise a file it will never write.
        **({"vsc": "output/profile/metaphlan/profiles/{sample}_vsc.csv"}
           if vsc_enabled() else {})
    retries:
        config['retries'].get('metaphlan', 2)
    resources:
        # Peaked at 29.7 GB against a 32 GB allocation on a 341-sample
        # gut cohort -- 93%, i.e. the next OOM waiting to happen. The
        # footprint is driven by the bowtie2 index and read volume,
        # both of which vary by sample.
        mem_mb=mem_escalate('metaphlan', base_default=32000),
        runtime=runtime_escalate('metaphlan', base_default=480)
    conda:
        "../env/profile.yaml"
    threads:
        config['threads']['metaphlan']
    params:
        db_name=config['params']['metaphlan']['db_name'],
        other=config['params']['metaphlan']['other'],
        # Empty when VSC profiling is off, so the same shell serves both.
        vsc_args=lambda wildcards: (
            "--profile_vsc --vsc_out "
            f"output/profile/metaphlan/profiles/{wildcards.sample}_vsc.csv"
            if vsc_enabled() else ""
        ),
        # MetaPhlAn takes its inputs as one comma-separated argument.
        read_arg=lambda wildcards, input: ",".join(input.reads)
    benchmark:
        "output/benchmarks/profile/metaphlan/{sample}_benchmark.txt"
    log:
        "output/logs/profile/metaphlan/{sample}.log"
    shell:
        """
        mkdir -p output/profile/metaphlan/bowtie2s output/profile/metaphlan/sams

        metaphlan {params.read_arg} \
        --input_type fastq \
        --nproc {threads} \
        --db_dir {input.db_path} \
        -x {params.db_name} \
        --mapout {output.bt2} \
        -s {output.sam} \
        -o {output.profile} \
        {params.vsc_args} \
        {params.other} \
        2> {log} 1>&2
        """

rule merge_metaphlan_tables:
    """

    Merges MetaPhlAn3 profiles into a single table.

    """
    input:
        expand(rules.metaphlan.output.profile,
               sample=samples)
    output:
        merged_abundance_table="output/profile/metaphlan/merged_abundance_table.txt"
    conda:
        "../env/profile.yaml"
    log:
        "output/logs/profile/metaphlan/merge_metaphlan_tables/merged_abundance_table.log"
    benchmark:
        "output/benchmarks/profile/metaphlan/merge_metaphlan_tables/merged_abundance_table_benchmark.txt"
    shell:
        """
        merge_metaphlan_tables.py {input} \
        -o {output.merged_abundance_table} \
        2> {log} 1>&2
        """
