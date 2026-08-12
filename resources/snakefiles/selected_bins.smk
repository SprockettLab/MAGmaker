import os
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
        search_engine = config['params']['das_tool']['search_engine'],
        score_threshold = config['params']['das_tool'].get(
            'score_threshold', 0.5),
        # Only bin sets that actually contain bins reach DAS_Tool. A binner may
        # legitimately decline a sample -- MaxBin2 stops when an assembly has
        # too few single-copy marker genes to seed its EM -- and passing it an
        # empty bin set aborts the run for every other binner too. Filtering
        # here means one binner declining costs that binner's contribution
        # rather than the whole sample.
        bin_sets = lambda wildcards, input: ",".join(
            path for path, label in zip(
                [input.metabat2[0], input.maxbin2[0], input.concoct[0]],
                ["metabat2", "maxbin2", "concoct"])
            if os.path.getsize(path) > 0),
        bin_labels = lambda wildcards, input: ",".join(
            label for path, label in zip(
                [input.metabat2[0], input.maxbin2[0], input.concoct[0]],
                ["metabat2", "maxbin2", "concoct"])
            if os.path.getsize(path) > 0)
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
        # DAS_Tool's parallel gene-prediction/SCG stage fails nondeterministically
        # under high per-node concurrency (a known DAS_Tool fragility,
        # cmks/DAS_Tool#110), with several signatures: "No SCGs detected", pullseq
        # "failed to open names file", "Gene prediction failed", and NFS "stale
        # file handle" reads of its own intermediates. All are transient and clear
        # on a re-run. So the rule treats any non-zero exit as retryable UNLESS the
        # log shows the one genuine-empty signature ("No bins with bin-score"),
        # i.e. a real sample with no bin above the score threshold. Enumerating the
        # empty case and retrying everything else avoids masking a crash as an
        # empty sample (which previously dropped richly-binnable samples). Partial
        # outputs from a failed attempt are cleared first so each retry is clean.
        """
            rm -rf {params.basename}_* {params.basename}.* 2> /dev/null || true

            # Every binner declined this sample. Nothing to reconcile, so
            # record it as empty rather than handing DAS_Tool no input.
            if [ -z "{params.bin_sets}" ]; then
                echo "No binner produced bins for this sample; no MAGs." >> {log}
                printf "bin_id\n" > {output.out}
                exit 0
            fi

            set +e
            DAS_Tool \
            --bins {params.bin_sets} \
            --contigs {input.contigs} \
            --outputbasename {params.basename} \
            --labels {params.bin_labels} \
            --write_bins \
            --write_bin_evals \
            --threads {threads} \
            --score_threshold {params.score_threshold} \
            --search_engine {params.search_engine} \
            2> {log} 1>&2
            status=$?
            set -e
            if [ $status -eq 0 ]; then
                # Drop the large gene-prediction/DIAMOND intermediates (whole-
                # sample proteins, chunked DasToolParallel tmp files, SCG hit
                # tables) that DAS_Tool leaves behind. They bloat run_DAS_Tool/
                # into thousands of files and slow NFS metadata ops for
                # downstream rules. Keep the summary, bins, and evals.
                rm -rf {params.basename}_proteins.faa* {params.basename}.seqlength \
                    2> /dev/null || true
                exit 0
            fi
            # Non-zero exit: only a genuine "no bin scored above the threshold" is
            # a real empty sample -- write a header and succeed. Everything else
            # (SCG race, pullseq/gene-prediction/NFS crashes in DAS_Tool's parallel
            # plumbing) is transient, so fail and let Snakemake retry rather than
            # masking a richly-binnable sample as empty.
            if grep -qE "No bins with bin-score" {log}; then
                # "No bins with bin-score" is ambiguous. It is the genuine
                # empty-sample signature, but it is ALSO how the transient SCG
                # race surfaces: when DIAMOND's single-copy-gene search comes
                # up empty, DAS_Tool does not abort -- it proceeds to scoring,
                # every bin scores 0 for want of markers, and the run reports
                # no bin above the threshold. Treating that as an empty sample
                # discards a fully binnable one.
                #
                # The two are told apart by whether the BACTERIAL SCG set came
                # back empty. "No SCGs detected for SCG set: archaea" is
                # unremarkable -- most samples have no archaea -- but an empty
                # bacterial set alongside real input bins means the search
                # failed, not that the sample is empty.
                #
                # Observed on a 341-sample gut cohort: 87 samples (25%) were
                # written off as empty this way, discarding 6430 candidate
                # bins. Their assemblies were not weak; median N50 was 16%
                # HIGHER than the samples that succeeded.
                if grep -qE "No SCGs detected for SCG set: bacteria" {log} \
                   && [ -s "$(echo {input.metabat2} | cut -d' ' -f1)" ]; then
                    echo "Bacterial SCG set empty despite non-empty input bins: DAS_Tool's SCG search failed rather than the sample being empty. Failing for retry." >> {log}
                    exit 1
                fi
                echo "No bins above score threshold; no MAGs for this sample." >> {log}
                printf "bin_id\n" > {output.out}
                exit 0
            fi
            echo "DAS_Tool failed before scoring (transient/infrastructure); failing for retry." >> {log}
            exit 1
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
    retries:
        config['retries']['consolidate_DAS_Tool_bins']
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


