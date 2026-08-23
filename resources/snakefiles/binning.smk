
def get_bam_list(sample, mapper, contig_pairings):
    fps = expand("output/mapping/{mapper}/sorted_bams/{contig_pairings}_Mapped_To_{sample}.bam",
    mapper = mapper,
    sample = sample,
    contig_pairings = contig_pairings[sample])
    return(fps)

def get_index_list(sample, mapper, contig_pairings):
    fps = expand("output/mapping/{mapper}/sorted_bams/{contig_pairings}_Mapped_To_{sample}.bam.bai",
    mapper = mapper,
    sample = sample,
    contig_pairings = contig_pairings[sample])
    return(fps)

rule make_metabat2_coverage_table:
    """
    Uses jgi_summarize_bam_contig_depths to generate a depth.txt file.
    """
    input:
        bams = lambda wildcards: get_bam_list(wildcards.contig_sample, wildcards.mapper, contig_pairings)
    output:
        coverage_table="output/binning/metabat2/{mapper}/coverage_tables/{contig_sample}_coverage_table.txt"
    conda:
        "../env/binning.yaml"
    benchmark:
        "output/benchmarks/binning/metabat2/{mapper}/make_metabat2_coverage_table/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/metabat2/{mapper}/make_metabat2_coverage_table/{contig_sample}.log"
    retries:
        config['retries'].get('make_metabat2_coverage_table', 2)
    resources:
        mem_mb=config['mem_mb'].get('make_metabat2_coverage_table', 64000),
        runtime=runtime_escalate('make_metabat2_coverage_table', base_default=720)
    shell:
        """
            jgi_summarize_bam_contig_depths --outputDepth {output.coverage_table} {input.bams} 2> {log}
        """

rule run_metabat2:
    """
    Runs Metabat2:
    MetaBAT2 clusters metagenomic contigs into different "bins", each of which should correspond to a putative genome.

    MetaBAT2 uses nucleotide composition information and source strain abundance (measured by depth-of-coverage by aligning the reads to the contigs) to perform binning.
    """
    input:
        contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
                assembler = config['assemblers'],
                contig_sample = wildcards.contig_sample),
        coverage_table = lambda wildcards: expand("output/binning/metabat2/{mapper}/coverage_tables/{contig_sample}_coverage_table.txt",
                mapper=config['mappers'],
                contig_sample=wildcards.contig_sample)
    output:
        bins = directory("output/binning/metabat2/{mapper}/run_metabat2/{contig_sample}/")
    params:
        basename = "output/binning/metabat2/{mapper}/run_metabat2/{contig_sample}/{contig_sample}_bin",
        min_contig_length = config['params']['metabat2']['min_contig_length'],
        # MetaBAT2 defaults --seed to 0, which its source then replaces with
        # time(0), so an unseeded run is not reproducible. Passed explicitly
        # so two runs of this pipeline on the same data agree.
        seed = config['params']['metabat2'].get(
            'seed', config.get('seed', 8675309)),
        extra = config['params']['metabat2']['extra']  # optional parameters
    threads:
        config['threads']['run_metabat2']
    conda:
        "../env/binning.yaml"
    retries:
        config['retries'].get('run_metabat2', 2)
    benchmark:
        "output/benchmarks/binning/metabat2/{mapper}/run_metabat2/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/metabat2/{mapper}/run_metabat2/{contig_sample}.log"
    shell:
        """
            mkdir -p {output.bins}

            # metabat2 exits non-zero with "[Error!] There were no large
            # target contigs. Cannot proceed." when a sample's assembly has
            # zero contigs >= --minContig. Same class of property-of-the-
            # sample failure already tolerated for run_maxbin2, run_semibin2,
            # and run_concoct -- confirmed 2026-08-22 on the same
            # Ferretti_2018/SAMN06350074 assembly that also tripped
            # concoct's equivalent check (0 contigs >=1500bp, only 1
            # >=1000bp). Retries can't fix a property of the input.
            #
            # Only that one cause is tolerated. Every other metabat2
            # failure stays fatal, same rationale as the other binners.
            if ! metabat2 {params.extra} --numThreads {threads} \
                --inFile {input.contigs} \
                --outFile {params.basename} \
                --abdFile {input.coverage_table} \
                --minContig {params.min_contig_length} \
                --seed {params.seed} \
                2> {log} 1>&2; then
                if grep -q "There were no large target contigs" {log}; then
                    echo "metabat2: assembly too sparse to bin (no contig >= {params.min_contig_length}bp); sample yields no bins" >> {log}
                    exit 0
                fi
                echo "metabat2 failed for a reason other than assembly sparsity" >> {log}
                exit 1
            fi
        """


rule make_maxbin2_coverage_table:
    """
       Commands to generate a coverage table using `samtools coverage` for input into maxbin2
    """
    input:
        bams="output/mapping/{mapper}/sorted_bams/{read_sample}_Mapped_To_{contig_sample}.bam"
    output:
        coverage_table="output/binning/maxbin2/{mapper}/coverage_tables/{read_sample}_Mapped_To_{contig_sample}_coverage.txt"
    conda:
        "../env/binning.yaml"
    benchmark:
        "output/benchmarks/binning/maxbin2/{mapper}/make_maxbin2_coverage_table/{read_sample}_Mapped_To_{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/maxbin2/{mapper}/make_maxbin2_coverage_table/{read_sample}_Mapped_To_{contig_sample}.log"
    retries:
        config['retries'].get('make_maxbin2_coverage_table', 3)
    resources:
        runtime=runtime_escalate('make_maxbin2_coverage_table', base_default=240)
    shell:
        """
          samtools coverage {input.bams} | \
          tail -n +2 | \
          sort -k1 | \
          cut -f1,6 > {output.coverage_table} 2> {log}
       """

rule make_maxbin2_abund_list:
    """
       Combines the file paths from 'make_maxbin2_coverage_table' for MaxBin2
    """
    input:
        lambda wildcards: expand("output/binning/maxbin2/{mapper}/coverage_tables/{read_sample}_Mapped_To_{contig_sample}_coverage.txt",
                mapper = wildcards.mapper,
                contig_sample = wildcards.contig_sample,
                read_sample = contig_pairings[wildcards.contig_sample])
    output:
        abund_list = "output/binning/maxbin2/{mapper}/abundance_lists/{contig_sample}_abund_list.txt"
    benchmark:
        "output/benchmarks/binning/maxbin2/{mapper}/make_maxbin2_abund_list/{contig_sample}_abund_list_benchmark.txt"
    log:
        "output/logs/binning/maxbin2/{mapper}/make_maxbin2_abund_list/{contig_sample}_abund_list.log"
    run:
        with open(output.abund_list, 'w') as f:
            for fp in input:
                f.write('%s\n' % fp)


rule run_maxbin2:
    """
    Runs MaxBin2:
    MaxBin2 clusters metagenomic contigs (assembled contiguous genome fragments) into different "bins", each of which corresponds to a putative population genome. It uses nucleotide composition information, source strain abundance (measured by depth-of-coverage by aligning the reads to the contigs), and phylogenetic marker genes to perform binning through an Expectation-Maximization (EM) algorithm.
    """
    input:
        contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
                assembler = config['assemblers'],
                contig_sample = wildcards.contig_sample),
        abund_list = lambda wildcards: expand("output/binning/maxbin2/{mapper}/abundance_lists/{contig_sample}_abund_list.txt",
                mapper=config['mappers'],
                contig_sample=wildcards.contig_sample)
    output:
        bins = directory("output/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}/")
    params:
        basename = "output/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}/{contig_sample}_bin",
        prob = config['params']['maxbin2']['prob_threshold'],  # optional parameters
        min_contig_length = config['params']['maxbin2']['min_contig_length'],
        extra = config['params']['maxbin2']['extra']  # optional parameters
    threads:
        config['threads']['run_maxbin2']
    conda:
        "../env/binning.yaml"
    benchmark:
        "output/benchmarks/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}.log"
    retries:
        config['retries'].get('run_maxbin2', 2)
    resources:
        # Never had an entry at all, in this rule or any profile -- silently
        # ran on the bare 120 min default since this project began. Confirmed
        # TIMEOUT 2026-08-20 (Pan_troglodytes_troglodytes/SAMN28679416, SLURM
        # job 1563464, CANCELLED ... DUE TO TIME LIMIT at exactly 120 min):
        # a 22-sample all-vs-all binning group means MaxBin2's own marker-
        # gene search (FragGeneScan + HMMER over ~100k contigs) plus EM
        # clustering across 15 abundance columns per contig, comparable in
        # scope to make_concoct_coverage_table's ~90-BAM case. No direct
        # timing measurement of a successful large-cohort run exists yet, so
        # this is a moderate, evidence-motivated starting point rather than
        # a guessed-large one -- retries + escalation cover the rest if it's
        # still not enough, rather than guessing higher up front.
        runtime=runtime_escalate('run_maxbin2', base_default=240)
    shell:
        """
            mkdir -p {output.bins}

            # MaxBin2 exits non-zero when an assembly carries too few
            # single-copy marker genes to seed its EM, reporting "the medium of
            # marker gene number <= 1". That is a property of the sample rather
            # than a failure: host-dominated libraries reach it routinely, and
            # metabat2 and concoct bin the same assembly without complaint.
            # DAS_Tool needs only one bin set, so the sample is left with no
            # maxbin2 bins instead of stopping the workflow. Every other
            # non-zero exit is still a real error.
            if ! run_MaxBin.pl -thread {threads} -prob_threshold {params.prob} \
            -min_contig_length {params.min_contig_length} {params.extra} \
            -contig {input.contigs} \
            -abund_list {input.abund_list} \
            -out {params.basename} \
            2> {log} 1>&2; then
                if grep -q "cannot be binned" {log}; then
                    echo "MaxBin2 found too few marker genes to bin this sample; continuing with no maxbin2 bins." >> {log}
                else
                    exit 1
                fi
            fi
        """

rule cut_up_fasta:
    """
    Cut up fasta file in non-overlapping or overlapping parts of equal length.
    Optionally creates a BED-file where the cutup contigs are specified in terms
    of the original contigs. This can be used as input to concoct_coverage_table.py.
    """
    input:
        contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
                assembler = config['assemblers'],
                contig_sample = wildcards.contig_sample)
    output:
        bed="output/binning/concoct/{mapper}/contigs_10K/{contig_sample}.bed",
        contigs_10K="output/binning/concoct/{mapper}/contigs_10K/{contig_sample}.fa"
    conda:
        "../env/concoct_linux.yaml"
    params:
        chunk_size=config['params']['concoct']['chunk_size'],
        overlap_size=config['params']['concoct']['overlap_size']
    benchmark:
        "output/benchmarks/binning/concoct/{mapper}/cut_up_fasta/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/concoct/{mapper}/cut_up_fasta/{contig_sample}.log"
    shell:
        """
          cut_up_fasta.py {input.contigs} \
          -c {params.chunk_size} \
          -o {params.overlap_size} \
          --merge_last \
          -b {output.bed} > {output.contigs_10K} 2> {log}
        """

rule make_concoct_coverage_table:
    """
    Generates table with per sample coverage depth.
    Assumes the directory "/output/binning/{mapper}/mapped_reads/" contains sorted and indexed bam files where each contig file has has reads mapped against it from the selected prototypes.

    """
    input:
        bed = "output/binning/concoct/{mapper}/contigs_10K/{contig_sample}.bed",
        bam = lambda wildcards: get_bam_list(wildcards.contig_sample, config['mappers'], contig_pairings),
        index = lambda wildcards: get_index_list(wildcards.contig_sample, config['mappers'], contig_pairings)
    output:
        coverage_table = "output/binning/concoct/{mapper}/coverage_tables/{contig_sample}_coverage_table.txt"
    conda:
        "../env/concoct_linux.yaml"
    benchmark:
        "output/benchmarks/binning/concoct/{mapper}/make_concoct_coverage_table/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/concoct/{mapper}/make_concoct_coverage_table/{contig_sample}.log"
    retries:
        config['retries'].get('make_concoct_coverage_table', 2)
    resources:
        mem_mb=config['mem_mb'].get('make_concoct_coverage_table', 64000),
        runtime=runtime_escalate('make_concoct_coverage_table', base_default=720)
    shell:
        """
          concoct_coverage_table.py {input.bed} \
          {input.bam} > {output.coverage_table} 2> {log}
        """

rule run_concoct:
    """
    CONCOCT - Clustering cONtigs with COverage and ComposiTion
    CONCOCT does unsupervised binning of metagenomic contigs by using nucleotide composition - kmer frequencies - and coverage data for multiple samples.
    """
    input:
        contigs_10K=rules.cut_up_fasta.output.contigs_10K,
        coverage_table=rules.make_concoct_coverage_table.output.coverage_table
    output:
        clustering = "output/binning/concoct/{mapper}/run_concoct/{contig_sample}/{contig_sample}_bins_clustering.csv"
    params:
        bins = "output/binning/concoct/{mapper}/run_concoct/{contig_sample}/{contig_sample,[A-Za-z0-9_]+}_bins",
        min_contig_length=config['params']['concoct']['min_contig_length']
    conda:
        "../env/concoct_linux.yaml"
    threads:
        config['threads']['run_concoct']
    benchmark:
        "output/benchmarks/binning/concoct/{mapper}/run_concoct/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/concoct/{mapper}/run_concoct/{contig_sample}.log"
    retries:
        config['retries'].get('run_concoct', 2)
    resources:
        runtime=runtime_escalate('run_concoct', base_default=720)
    shell:
        """
            mkdir -p output/binning/concoct/{wildcards.mapper}/run_concoct/{wildcards.contig_sample}

            # concoct exits non-zero with "Not enough contigs pass the
            # threshold filter" when a sample's assembly is too sparse for
            # even one contig to clear -l {params.min_contig_length} --
            # the same class of property-of-the-sample failure already
            # tolerated for run_maxbin2 ("cannot be binned") and
            # run_semibin2 ("no must-link pairs"/"basepairs"), confirmed
            # 2026-08-22 on Ferretti_2018/SAMN06350074 after two identical
            # retries (retries can't fix a property of the input). concoct's
            # own clustering CSV format is headerless contig,cluster pairs,
            # so an empty file is a well-formed zero-contigs-clustered
            # result, not a guessed schema -- merge_cutup_clustering.py is
            # concoct's own script and should treat it as the normal empty
            # case.
            #
            # Only that one cause is tolerated. Every other concoct failure
            # stays fatal, same rationale as the maxbin2/semibin2 branches.
            if ! concoct --threads {threads} -l {params.min_contig_length} \
                --composition_file {input.contigs_10K} \
                --coverage_file {input.coverage_table} \
                -b {params.bins} \
                2> {log} 1>&2; then
                if grep -q "Not enough contigs pass the threshold filter" {log}; then
                    echo "concoct: assembly too sparse to bin (no contig >= {params.min_contig_length}bp); sample yields no bins" >> {log}
                    : > {output.clustering}
                    exit 0
                fi
                echo "concoct failed for a reason other than assembly sparsity" >> {log}
                exit 1
            fi

            mv output/binning/concoct/{wildcards.mapper}/run_concoct/{wildcards.contig_sample}/{wildcards.contig_sample}_bins_clustering_gt{params.min_contig_length}.csv output/binning/concoct/{wildcards.mapper}/run_concoct/{wildcards.contig_sample}/{wildcards.contig_sample}_bins_clustering.csv
        """

rule merge_cutup_clustering:
    """
    Merges subcontig clustering into original contig clustering.
    """
    input:
        bins = lambda wildcards: expand("output/binning/concoct/{mapper}/run_concoct/{contig_sample}/{contig_sample}_bins_clustering.csv",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        merged = "output/binning/concoct/{mapper}/merge_cutup_clustering/{contig_sample}_clustering_merged.csv"
    conda:
        "../env/concoct_linux.yaml"
    benchmark:
        "output/benchmarks/binning/concoct/{mapper}/merge_cutup_clustering/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/concoct/{mapper}/merge_cutup_clustering/{contig_sample}.log"
    shell:
        """
            merge_cutup_clustering.py {input.bins} > {output.merged} 2> {log}
        """

rule extract_fasta_bins:
    """
    Extracts bins as individual FASTA.
    """
    input:
        original_contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
                    assembler = config['assemblers'],
                    contig_sample = wildcards.contig_sample),
        clustering_merged = rules.merge_cutup_clustering.output.merged
    output:
        fasta_bins = directory("output/binning/concoct/{mapper}/extract_fasta_bins/{contig_sample}_bins/")
    conda:
        "../env/concoct_linux.yaml"
    benchmark:
        "output/benchmarks/binning/concoct/{mapper}/extract_fasta_bins/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/concoct/{mapper}/extract_fasta_bins/{contig_sample}.log"
    shell:
        """
            mkdir -p {output.fasta_bins}
            extract_fasta_bins.py \
            {input.original_contigs} \
            {input.clustering_merged} \
            --output_path {output.fasta_bins} \
            2> {log}
        """


rule run_semibin2:
    """
    Bins contigs with SemiBin2, a fourth binner using self-supervised deep
    learning over composition and multi-sample coverage.

    Why a fourth binner at all: consolidation can only gain from a
    candidate no other binner produced, so the value is in failing
    DIFFERENTLY rather than in failing less. MetaBAT2, MaxBin2 and CONCOCT
    are all composition-plus-coverage clustering with different distance
    measures; a learned embedding is a different mechanism.

    It takes the same BAMs the other binners' coverage tables are built
    from, so multi-sample differential coverage -- which is the thing this
    pipeline's prototype mapping exists to provide -- reaches it unchanged
    and no new mapping is needed.

    SEED. SemiBin2 documents --random-seed as reproducing results across
    runs, unlike VAMB which states determinism is not guaranteed even when
    seeded. That claim is why this tool was chosen and it is not taken on
    trust: run the binner three times and compare before using it for
    anything. --engine cpu is set for the same reason as much as for the
    absence of a GPU, since GPU kernels are a common source of run-to-run
    variation.

    Its output layout differs between versions and options, so the bins are
    normalised into the rule's own output directory rather than leaving
    downstream rules to guess which of output_bins,
    output_recluster_bins or output_prerecluster_bins was written.
    """
    input:
        contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
                assembler = config['assemblers'],
                contig_sample = wildcards.contig_sample),
        bams = lambda wildcards: get_bam_list(wildcards.contig_sample, wildcards.mapper, contig_pairings),
        bais = lambda wildcards: get_index_list(wildcards.contig_sample, wildcards.mapper, contig_pairings)
    output:
        bins = directory("output/binning/semibin2/{mapper}/run_semibin2/{contig_sample}/")
    params:
        work = "output/binning/semibin2/{mapper}/work/{contig_sample}",
        seed = config['params'].get('semibin2', {}).get(
            'random_seed', config.get('seed', 8675309)),
        engine = config['params'].get('semibin2', {}).get('engine', 'cpu'),
        environment = config['params'].get('semibin2', {}).get('environment', ''),
        min_len = config['params'].get('semibin2', {}).get('min_contig_length', 1000),
        extra = config['params'].get('semibin2', {}).get('extra', '')
    threads:
        config['threads'].get('run_semibin2', 16)
    resources:
        mem_mb = mem_escalate('run_semibin2', base_default=32000),
        runtime = runtime_escalate('run_semibin2', base_default=360)
    retries:
        config['retries'].get('run_semibin2', 2)
    conda:
        "../env/semibin.yaml"
    benchmark:
        "output/benchmarks/binning/semibin2/{mapper}/run_semibin2/{contig_sample}_benchmark.txt"
    log:
        "output/logs/binning/semibin2/{mapper}/run_semibin2/{contig_sample}.log"
    shell:
        """
            rm -rf {params.work}
            mkdir -p {params.work} {output.bins}

            # With --environment SemiBin2 uses a pretrained model and skips
            # training entirely, which is faster and removes the stochastic
            # step. Without it the model is trained from this sample, which
            # makes no assumption about which published habitat the data
            # resembles. Left empty by default: the built-in models are
            # human, dog, cat, mouse, pig, chicken, ocean, soil and similar,
            # and asserting that a wild primate gut is one of those is a
            # claim about the biology, not a tuning choice.
            if [ -n "{params.environment}" ]; then
                MODEL="--environment {params.environment}"
            else
                MODEL="--self-supervised"
            fi

            # SemiBin2 needs at least one contig of >=4000 bp to form
            # must-link pairs, and exits non-zero when an assembly has none.
            # That is a property of the sample, not a pipeline failure: the
            # sample yields no bins, exactly like the empty-output case
            # below. Without this branch the non-zero exit propagates under
            # `bash -euo pipefail` and takes down the whole arm, and the
            # graceful "no bin directory" path further down is never
            # reached. Retries cannot help -- the input is the problem, so
            # all of them fail identically.
            #
            # Two distinct messages have now been seen for the same
            # underlying "assembly too sparse to bin" situation, confirmed
            # 2026-08-22 on Ferretti_2018/SAMN06350118: "no must-link pairs
            # can be generated" (no contig >=4000bp) and "but only N
            # contain(s) at least 1000 basepairs" (an assembly of just 7
            # contigs, only 1 over even the 1000bp floor). The second
            # message wasn't matched by the original single-string check,
            # so this exact tolerated case was being retried as fatal --
            # pointlessly, since retries can't fix a property of the input.
            #
            # Only these two causes are tolerated. Every other SemiBin2
            # failure stays fatal, because a silent `|| true` here would
            # turn real crashes into samples that quietly contribute
            # nothing.
            if ! SemiBin2 single_easy_bin \
                -i {input.contigs} \
                -b {input.bams} \
                -o {params.work} \
                --threads {threads} \
                --min-len {params.min_len} \
                --random-seed {params.seed} \
                --engine {params.engine} \
                ${{MODEL}} {params.extra} \
                2> {log} 1>&2; then
                if grep -qE "no must-link pairs can be generated|contain\(s\) at least [0-9]+ basepairs" {log}; then
                    echo "SemiBin2: assembly too fragmented/sparse to bin; sample yields no bins" >> {log}
                    rm -rf {params.work}
                    exit 0
                fi
                echo "SemiBin2 failed for a reason other than assembly fragmentation" >> {log}
                exit 1
            fi

            # Normalise: whichever directory this version wrote into, the
            # bins end up in the rule's declared output as plain .fa.
            SRC=""
            for d in output_bins output_recluster_bins output_prerecluster_bins bins; do
                if [ -d "{params.work}/$d" ]; then SRC="{params.work}/$d"; break; fi
            done
            if [ -z "$SRC" ]; then
                echo "SemiBin2 produced no bin directory; sample yields no bins" >> {log}
                exit 0
            fi
            n=0
            for f in "$SRC"/*.fa "$SRC"/*.fa.gz "$SRC"/*.fna "$SRC"/*.fna.gz; do
                [ -e "$f" ] || continue
                b=$(basename "$f"); b=${{b%.gz}}; b=${{b%.fna}}; b=${{b%.fa}}
                case "$f" in
                    *.gz) gunzip -c "$f" > {output.bins}/{wildcards.contig_sample}_$b.fa ;;
                    *)    cp "$f" {output.bins}/{wildcards.contig_sample}_$b.fa ;;
                esac
                n=$((n+1))
            done
            echo "normalised $n bins from $SRC" >> {log}
            rm -rf {params.work}
        """
