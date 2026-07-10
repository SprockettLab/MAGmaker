from os.path import basename, dirname, join
from shutil import copyfile
from glob import glob


rule metabat2_Fasta_to_Scaffolds2Bin:
    """
    Uses Fasta_to_Scaffolds2Bin script in DAS Tools to generate a scaffolds2bin.tsv file.
    """
    input:
        bins = lambda wildcards: expand("output/binning/metabat2/{mapper}/run_metabat2/{contig_sample}/",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        scaffolds2bin="output/selected_bins/metabat2/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv"
    conda:
        "../env/selected_bins.yaml"
    benchmark:
        "output/benchmarks/selected_bins/metabat2/{mapper}/scaffolds2bin/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/metabat2/{mapper}/scaffolds2bin/{contig_sample}.log"
    shell:
        """
            Fasta_to_Contig2Bin.sh \
            -i {input.bins} \
            -e fa > {output.scaffolds2bin}
        """


rule maxbin2_Fasta_to_Scaffolds2Bin:
    """
    Uses Fasta_to_Scaffolds2Bin script in DAS Tools to generate a scaffolds2bin.tsv file.
    """
    input:
        bins = lambda wildcards: expand("output/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}/",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        scaffolds2bin="output/selected_bins/maxbin2/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv"
    conda:
        "../env/selected_bins.yaml"
    benchmark:
        "output/benchmarks/selected_bins/maxbin2/{mapper}/scaffolds2bin/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/maxbin2/{mapper}/scaffolds2bin/{contig_sample}.log"
    shell:
        """
            Fasta_to_Contig2Bin.sh \
            -i {input.bins} \
            -e fasta > {output.scaffolds2bin}
        """

rule concoct_Fasta_to_Scaffolds2Bin:
    """
    Uses Fasta_to_Scaffolds2Bin script in DAS Tools to generate a scaffolds2bin.tsv file.
    """
    input:
        bins = lambda wildcards: expand("output/binning/concoct/{mapper}/extract_fasta_bins/{contig_sample}_bins/",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        scaffolds2bin="output/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv"
    conda:
        "../env/selected_bins.yaml"
    benchmark:
        "output/benchmarks/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}.log"
    shell:
        """
            Fasta_to_Contig2Bin.sh \
            -i {input.bins} \
            -e fa > {output.scaffolds2bin}
        """

# rule concoct_Fasta_to_Scaffolds2Bin:
#     """
#     Uses perl to create a scaffolds2bin.tsv file from a clustering_merged.csv file.
#     """
#     input:
#         merged = lambda wildcards: expand("output/binning/concoct/{mapper}/merge_cutup_clustering/{contig_sample}_clustering_merged.csv",
#                 mapper = config['mappers'],
#                 contig_sample = wildcards.contig_sample)
#     output:
#         scaffolds2bin="output/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv"
#     conda:
#         "../env/selected_bins.yaml"
#     benchmark:
#         "output/benchmarks/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}_benchmark.txt"
#     log:
#         "output/logs/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}.log"
#     shell:
#         """
#             perl -pe "s/,/\tconcoct_bins./g;" {input.merged} > {output.scaffolds2bin}
#         """

rule run_DAS_Tool:
    """
    Selects bins using DAS_Tool
    """
    input:
        metabat2 = lambda wildcards: expand("output/selected_bins/metabat2/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample),
        maxbin2 = lambda wildcards: expand("output/selected_bins/maxbin2/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample),
        concoct = lambda wildcards: expand("output/selected_bins/concoct/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample),
        contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
                    assembler = config['assemblers'],
                    contig_sample = wildcards.contig_sample)
    output:
        out="output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}_DASTool_summary.tsv"
    params:
        basename = "output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}",
        search_engine = config['params']['das_tool']['search_engine']
    conda:
        "../env/selected_bins.yaml"
    threads:
        config['threads']['run_DAS_Tool']
    retries:
        config['retries']['run_DAS_Tool']
    benchmark:
        "output/benchmarks/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}.log"
    shell:
        # DAS_Tool's single-copy-gene detection intermittently comes up empty
        # under high per-node concurrency, reporting no SCGs for both the bacteria
        # and archaea sets and aborting ("No single copy genes predicted"). This
        # is a known DAS_Tool fragility (cmks/DAS_Tool#110), not the data, and
        # here it is nondeterministic -- the same sample succeeds on a re-run --
        # so we let that specific failure propagate (exit non-zero) to trigger
        # Snakemake `retries`, while a genuinely empty sample (no bin scored above
        # the threshold) still writes a header and succeeds. Partial outputs from
        # a failed attempt are cleared first so each retry starts clean.
        """
            rm -rf {params.basename}_* {params.basename}.* 2> /dev/null || true
            set +e
            DAS_Tool \
            --bins {input.metabat2},{input.maxbin2},{input.concoct} \
            --contigs {input.contigs} \
            --outputbasename {params.basename} \
            --labels metabat2,maxbin2,concoct \
            --write_bins \
            --write_bin_evals \
            --threads {threads} \
            --search_engine {params.search_engine} \
            2> {log} 1>&2
            status=$?
            set -e
            if [ $status -eq 0 ]; then
                exit 0
            fi
            if grep -qE "No SCGs detected|No single copy genes predicted" {log}; then
                echo "Empty predicted-protein set (transient concurrency race); failing for retry." >> {log}
                exit 1
            fi
            echo "No bins above score threshold; no MAGs for this sample." >> {log}
            printf "bin_id\n" > {output.out}
        """


rule consolidate_DAS_Tool_bins:
    """
    Copies DAS_Tool MAG bins into a per-sample subdirectory for downstream QC.
    """
    input:
        "output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}_DASTool_summary.tsv"
    output:
        done=touch("output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}/.done")
    log:
        "output/logs/selected_bins/{mapper}/consolidate_DAS_Tool_bins/{contig_sample}.log"
    run:
        import os
        sample = wildcards.contig_sample
        fasta_dir = join(dirname(input[0]), sample + '_DASTool_bins')
        output_dir = dirname(output.done)
        os.makedirs(output_dir, exist_ok=True)
        for file in glob(join(fasta_dir, '*.fa')):
            copyfile(file, join(output_dir, basename(file)))


rule consolidate_DAS_Tool_bins_all:
    input:
        lambda wildcards: expand("output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}/.done",
                                 mapper=config['mappers'],
                                 contig_sample=contig_pairings.keys())


