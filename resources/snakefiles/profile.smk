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
    resources:
        # kraken2 reads a ~200 GB index per job, so the cost is dominated by
        # that read rather than by classification. Left uncapped, every
        # sample in an arm competes for the same shared storage at once and
        # jobs are cancelled on wall clock while doing almost no work. This
        # is a counter, not a reservation: set `resources: [kraken_slots=N]`
        # in the profile to allow N concurrent kraken jobs per workflow.
        kraken_slots=1
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
    shell:
        """
        metaphlan --install \
            --bowtie2db {output} \
            --index latest \
            2> {log} 1>&2
        """

rule metaphlan:
    """
    Performs taxonomic profiling using MetaPhlAn4.
    """
    input:
        reads=lambda wildcards: nonhost_reads(wildcards.sample),
        db_path=rules.download_metaphlan_db.output
    output:
        bt2="output/profile/metaphlan/bowtie2s/{sample}.bowtie2.bz2",
        sam="output/profile/metaphlan/sams/{sample}.sam.bz2",
        profile="output/profile/metaphlan/profiles/{sample}.txt"
    retries:
        config['retries'].get('metaphlan', 2)
    resources:
        # Peaked at 29.7 GB against a 32 GB allocation on a 341-sample
        # gut cohort -- 93%, i.e. the next OOM waiting to happen. The
        # footprint is driven by the bowtie2 index and read volume,
        # both of which vary by sample.
        mem_mb=mem_escalate('metaphlan', base_default=32000)
    conda:
        "../env/profile.yaml"
    threads:
        config['threads']['metaphlan']
    params:
        db_name=config['params']['metaphlan']['db_name'],
        other=config['params']['metaphlan']['other'],
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
        --bowtie2db {input.db_path} \
        --index {params.db_name} \
        --bowtie2out {output.bt2} \
        -s {output.sam} \
        -o {output.profile} \
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
