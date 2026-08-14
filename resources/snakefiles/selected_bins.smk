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


# =============================================================================
# Binette -- the alternative to DAS_Tool, selected by params.consolidation.tool
#
# Same job, different estimator. DAS_Tool scores candidate bins with 51
# bacterial single-copy genes and a redundancy heuristic; Binette scores them
# with CheckM2, which is the same estimator that judges the result in mag_qc.
# Under DAS_Tool the selection is made with the weakest instrument available
# and then audited with the strongest ones, and the two can disagree.
#
# Binette also generates candidates differently. For every pair of input bins
# sharing at least one contig it forms the intersection, both differences and
# the union, scores all of them, and greedily keeps a non-redundant set. Its
# score is `completeness - contamination_weight * contamination`, default
# weight 2, so the preference for clean genomes over complete ones is a
# stated parameter rather than a property of whichever marker set a tool
# happens to ship. See docs/bin_consolidation_options.md in the analysis repo.
#
# It consumes the same scaffolds2bin tables DAS_Tool does, so nothing above
# this point in the pipeline changes or reruns.
# =============================================================================

rule stage_binette_inputs:
    """
    Collects a sample's per-binner scaffolds2bin tables under binner-named
    filenames for Binette.

    Binette has no --labels option. It infers a name per bin set by stripping
    the common prefix and suffix from the input paths, and MAGmaker's tables
    differ only in a middle path component while sharing the same basename
    (<sample>_scaffolds2bin.tsv). Staging them as <binner>.tsv in one
    directory makes the inferred names exactly metabat2/maxbin2/concoct,
    which is what lets `origin` in Binette's report populate Winning_Binner.

    Empty tables are skipped rather than staged. A binner may legitimately
    decline a sample -- MaxBin2 stops when an assembly has too few marker
    genes to seed its EM -- and an empty bin set is not the same as a bin set
    with no bins in it.
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
                contig_sample = wildcards.contig_sample)
    output:
        staged = directory("output/selected_bins/{mapper}/binette_input/{contig_sample}")
    log:
        "output/logs/selected_bins/{mapper}/stage_binette_inputs/{contig_sample}.log"
    run:
        import os
        from shutil import copyfile
        os.makedirs(output.staged, exist_ok=True)
        staged = []
        with open(log[0], 'w') as fh:
            for label, paths in (('metabat2', input.metabat2),
                                 ('maxbin2', input.maxbin2),
                                 ('concoct', input.concoct)):
                path = paths[0]
                if os.path.getsize(path) == 0:
                    fh.write("%s produced no bins; not staged\n" % label)
                    continue
                copyfile(path, os.path.join(output.staged, label + '.tsv'))
                staged.append(label)
            fh.write("staged: %s\n" % (", ".join(staged) or "none"))


rule run_binette:
    """
    Reconciles the per-binner bin sets into one set of MAGs using CheckM2.
    """
    input:
        staged = "output/selected_bins/{mapper}/binette_input/{contig_sample}",
        contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
                    assembler = config['assemblers'],
                    contig_sample = wildcards.contig_sample)
    output:
        report = "output/selected_bins/{mapper}/run_binette/{contig_sample}/final_bins_quality_reports.tsv",
        contig2bin = "output/selected_bins/{mapper}/run_binette/{contig_sample}/final_contig_to_bin.tsv"
    params:
        out_dir = "output/selected_bins/{mapper}/run_binette/{contig_sample}",
        # Binette's own default is "bin", which would make bin names collide
        # across samples in the pooled mag_summary. Prefixing with the sample
        # keeps every MAG name unique without a rename step.
        prefix = "{contig_sample}",
        db_path = config['params']['checkm2']['db_path'],
        contamination_weight = config['params']['binette']['contamination_weight'],
        min_completeness = config['params']['binette']['min_completeness'],
        max_contamination = config['params']['binette']['max_contamination'],
        min_length = config['params']['binette']['min_length'],
        max_length = config['params']['binette']['max_length'],
        extra = config['params']['binette'].get('extra', '')
    threads:
        config['threads']['run_binette']
    resources:
        mem_mb = mem_escalate('run_binette', base_default=64000)
    retries:
        config['retries'].get('run_binette', 2)
    conda:
        "../env/binette.yaml"
    benchmark:
        "output/benchmarks/selected_bins/{mapper}/run_binette/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/{mapper}/run_binette/{contig_sample}.log"
    shell:
        # The empty case is written out rather than left to fail, mirroring
        # run_DAS_Tool: a sample every binner declined is a real outcome and
        # must be distinguishable from a crash. Headers match Binette's own
        # so make_mag_summary reads both the same way.
        """
            tables=$(ls {input.staged}/*.tsv 2>/dev/null || true)
            if [ -z "$tables" ]; then
                echo "No binner produced bins for this sample; no MAGs." >> {log}
                mkdir -p {params.out_dir}/final_bins
                printf "name\torigin\tis_original\toriginal_name\tcompleteness\tcontamination\tcheckm2_model\tscore\tsize\tN50\tcoding_density\tcontig_count\n" > {output.report}
                printf "contig\tbin\n" > {output.contig2bin}
                exit 0
            fi

            db_flag=""
            if [ -n "{params.db_path}" ]; then
                db_flag="--checkm2_db {params.db_path}"
            fi

            # Cleared rather than resumed. --resume reuses temporary_files/,
            # which after a killed job may be half-written, and a retry that
            # silently reuses a truncated DIAMOND result is worse than one
            # that redoes the work.
            rm -rf {params.out_dir}

            binette \
                --contig2bin_tables $tables \
                --contigs {input.contigs} \
                --outdir {params.out_dir} \
                --prefix {params.prefix} \
                --threads {threads} \
                --contamination_weight {params.contamination_weight} \
                --min_completeness {params.min_completeness} \
                --max_contamination {params.max_contamination} \
                --min_length {params.min_length} \
                --max_length {params.max_length} \
                $db_flag {params.extra} \
                2> {log} 1>&2

            # Binette writes no contig2bin table when it selects nothing.
            # An absent file would fail the rule; an empty one records that
            # the run succeeded and found nothing.
            if [ ! -s {output.contig2bin} ]; then
                printf "contig\tbin\n" > {output.contig2bin}
            fi

            # Self-clean. temporary_files/ holds the assembly proteins and
            # DIAMOND output, tens of GB across an arm, and directory bloat
            # on isilon has slowed this pipeline before.
            rm -rf {params.out_dir}/temporary_files
        """


rule consolidate_binette_bins:
    """
    Copies Binette MAG bins into a per-sample subdirectory for downstream QC.
    """
    input:
        report = "output/selected_bins/{mapper}/run_binette/{contig_sample}/final_bins_quality_reports.tsv"
    output:
        done = touch("output/selected_bins/{mapper}/Binette_Fastas/{contig_sample}/.done")
    log:
        "output/logs/selected_bins/{mapper}/consolidate_binette_bins/{contig_sample}.log"
    retries:
        config['retries']['consolidate_DAS_Tool_bins']
    run:
        import os
        sample = wildcards.contig_sample
        fasta_dir = join(dirname(input.report), 'final_bins')
        output_dir = dirname(output.done)
        os.makedirs(output_dir, exist_ok=True)
        copied = 0
        for file in glob(join(fasta_dir, '*.fa')):
            copyfile(file, join(output_dir, basename(file)))
            copied += 1
        with open(log[0], 'w') as fh:
            fh.write("%d bins copied from %s\n" % (copied, fasta_dir))


rule consolidate_binette_bins_all:
    input:
        lambda wildcards: expand("output/selected_bins/{mapper}/Binette_Fastas/{contig_sample}/.done",
                                 mapper=config['mappers'],
                                 contig_sample=contig_pairings.keys())




rule semibin2_Fasta_to_Scaffolds2Bin:
    """
    Contig-to-bin table for SemiBin2's bins, in the format both DAS_Tool and
    Binette consume.
    """
    input:
        bins = lambda wildcards: expand("output/binning/semibin2/{mapper}/run_semibin2/{contig_sample}/",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        scaffolds2bin="output/selected_bins/semibin2/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv"
    conda:
        "../env/selected_bins.yaml"
    benchmark:
        "output/benchmarks/selected_bins/semibin2/{mapper}/scaffolds2bin/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/semibin2/{mapper}/scaffolds2bin/{contig_sample}.log"
    shell:
        """
            # A sample where SemiBin2 produced nothing gets an empty table
            # rather than a failure. Both consolidation tools already treat
            # an empty bin set as a binner declining the sample.
            if ! ls {input.bins}/*.fa 2>/dev/null 1>/dev/null; then
                echo "no SemiBin2 bins for this sample" > {log}
                : > {output.scaffolds2bin}
            else
                Fasta_to_Contig2Bin.sh \
                -i {input.bins} \
                -e fa > {output.scaffolds2bin} 2> {log}
            fi
        """
