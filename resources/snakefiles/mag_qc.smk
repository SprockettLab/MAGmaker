import os
import glob as _glob
from os.path import join, dirname, basename


rule run_checkm2:
    input:
        done=f"output/selected_bins/{{mapper}}/{SELECTED_FASTAS}/{{contig_sample}}/.done"
    output:
        report="output/mag_qc/checkm2/{mapper}/{contig_sample}/quality_report.tsv"
    params:
        bins_dir=f"output/selected_bins/{{mapper}}/{SELECTED_FASTAS}/{{contig_sample}}",
        out_dir="output/mag_qc/checkm2/{mapper}/{contig_sample}",
        db_path=config['params']['checkm2']['db_path']
    threads:
        config['threads']['checkm2']
    retries:
        # checkm2 runs DIAMOND, which fails transiently under concurrency on
        # this cluster the same way gunc and DAS_Tool do. It was the only
        # DIAMOND-using rule without retries, and one such failure killed a
        # whole comparison run partway through.
        config['retries'].get('run_checkm2', 3)
    resources:
        mem_mb=config['mem_mb']['checkm2'],
        runtime=runtime_escalate('run_checkm2', base_default=360)
    conda:
        "../env/checkm2.yaml"
    log:
        "output/logs/mag_qc/checkm2/{mapper}/{contig_sample}.log"
    benchmark:
        "output/benchmarks/mag_qc/checkm2/{mapper}/{contig_sample}.txt"
    shell:
        """
        if ! ls {params.bins_dir}/*.fa 2>/dev/null 1>/dev/null; then
            echo "No bins in {params.bins_dir}; skipping checkm2." >> {log}
            printf "Name\tCompleteness\tContamination\n" > {output.report}
        else
            # Clear the output directory OURSELVES rather than letting
            # checkm2 --force do it. --force calls shutil.rmtree, which walks
            # the tree and unlinks entry by entry; on isilon a file can
            # vanish between the listing and the unlink and rmtree dies with
            # FileNotFoundError on its own intermediate, typically
            # protein_files. Observed killing a run outright.
            #
            # Retrying does not help, because every attempt meets the same
            # half-deleted directory. run_gunc has cleared its output for
            # this reason since the DIAMOND concurrency work; checkm2 was
            # missed.
            rm -rf {params.out_dir}
            mkdir -p {params.out_dir}

            db_flag=""
            if [ -n "{params.db_path}" ]; then
                db_flag="--database_path {params.db_path}"
            fi
            checkm2 predict \
                --input {params.bins_dir} \
                --output-directory {params.out_dir} \
                --threads {threads} \
                --extension fa \
                --force \
                $db_flag \
                2> {log} 1>&2
        fi
        """


rule run_gunc:
    input:
        done=f"output/selected_bins/{{mapper}}/{SELECTED_FASTAS}/{{contig_sample}}/.done"
    output:
        done=touch("output/mag_qc/gunc/{mapper}/{contig_sample}/.done")
    params:
        bins_dir=f"output/selected_bins/{{mapper}}/{SELECTED_FASTAS}/{{contig_sample}}",
        out_dir="output/mag_qc/gunc/{mapper}/{contig_sample}",
        db_flag=lambda wildcards: (
            f"--db_file {config['params']['gunc']['db_path']}"
            if config['params']['gunc']['db_path']
            else ""
        )
    threads:
        config['threads']['gunc']
    retries:
        config['retries']['run_gunc']
    resources:
        runtime=runtime_escalate('run_gunc', base_default=240)
    conda:
        "../env/gunc.yaml"
    log:
        "output/logs/mag_qc/gunc/{mapper}/{contig_sample}.log"
    benchmark:
        "output/benchmarks/mag_qc/gunc/{mapper}/{contig_sample}.txt"
    shell:
        """
        if ! ls {params.bins_dir}/*.fa 2>/dev/null 1>/dev/null; then
            echo "No bins in {params.bins_dir}; skipping gunc." >> {log}
        else
            # Clear any half-written output so a retry (gunc calls DIAMOND, which
            # races under concurrency like run_DAS_Tool) starts from a clean dir.
            rm -rf {params.out_dir}
            mkdir -p {params.out_dir}
            gunc run \
                --input_dir {params.bins_dir} \
                --out_dir {params.out_dir} \
                --threads {threads} \
                --file_suffix .fa \
                {params.db_flag} \
                2> {log} 1>&2
        fi
        """


# =============================================================================
# GTDB-Tk, as its three real stages
#
#   identify  Prodigal gene calls, then an HMM search for the bac120/ar53
#             markers
#   align     the masked concatenated alignment, and with it MSA_Percent
#   classify  pplacer placement in the reference tree
#
# Split because the stages differ enormously in cost and in what they are
# for. GTDB-Tk's documentation puts the bacterial memory requirement at
# about 140 GB and attributes it to pplacer, which only classify runs;
# identify and align fit in a fraction of that. Splitting them means an
# analysis that needs MSA_Percent but not taxonomy never schedules a job
# that only the largest nodes can take.
#
# The other reason is resumability. --skip-gtdbtk-classify then a later run
# without it adds ONLY classify: identify and align are already done and
# Snakemake sees them as up to date. Running this as one rule could not do
# that -- a single sentinel is either present or absent, so the second run
# would either redo all three stages or, worse, do nothing at all.
#
# Existing runs made before the split have classify_wf's output but none of
# the per-stage sentinels. Both upstream rules therefore detect complete
# output and reuse it instead of recomputing, and stamp the sentinel with
# that output's own mtime rather than the current time. Without the stamp a
# freshly created sentinel would be newer than the .done beside it and would
# trigger a pointless re-classification of every finished arm.
# =============================================================================

rule gtdbtk_identify:
    """
    Calls genes and finds the GTDB marker genes in each MAG.
    """
    input:
        done=f"output/selected_bins/{{mapper}}/{SELECTED_FASTAS}/{{contig_sample}}/.done"
    output:
        done=touch("output/mag_qc/gtdbtk/{mapper}/{contig_sample}/.identify.done")
    params:
        bins_dir=f"output/selected_bins/{{mapper}}/{SELECTED_FASTAS}/{{contig_sample}}",
        out_dir="output/mag_qc/gtdbtk/{mapper}/{contig_sample}",
        db_path=config['params']['gtdbtk']['db_path']
    threads:
        config['threads']['gtdbtk']
    retries:
        config['retries'].get('run_gtdbtk', 2)
    resources:
        mem_mb=mem_escalate('gtdbtk_align', base_default=32000, cap_multiple=2),
        runtime=runtime_escalate('gtdbtk_identify', base_default=480)
    conda:
        "../env/gtdbtk.yaml"
    log:
        "output/logs/mag_qc/gtdbtk/{mapper}/{contig_sample}.identify.log"
    benchmark:
        "output/benchmarks/mag_qc/gtdbtk/{mapper}/{contig_sample}.identify.txt"
    shell:
        """
        if ! ls {params.bins_dir}/*.fa 2>/dev/null 1>/dev/null; then
            echo "No bins in {params.bins_dir}; skipping gtdbtk identify." >> {log}
            exit 0
        fi

        # Output from a classify_wf run made before this rule existed is
        # complete and correct; recomputing it would cost hours for nothing.
        existing=$(ls {params.out_dir}/identify/*markers_summary.tsv 2>/dev/null | head -1 || true)
        if [ -n "$existing" ]; then
            echo "Reusing existing identify output in {params.out_dir}/identify" >> {log}
            touch -r "$existing" {output.done}
            exit 0
        fi

        export GTDBTK_DATA_PATH="{params.db_path}"
        gtdbtk identify \
            --genome_dir {params.bins_dir} \
            --out_dir {params.out_dir} \
            --cpus {threads} \
            --extension fa \
            2> {log} 1>&2
        """


rule gtdbtk_align:
    """
    Builds the concatenated marker alignment, which is where MSA_Percent
    comes from.
    """
    input:
        identify="output/mag_qc/gtdbtk/{mapper}/{contig_sample}/.identify.done"
    output:
        done=touch("output/mag_qc/gtdbtk/{mapper}/{contig_sample}/.align.done")
    params:
        bins_dir=f"output/selected_bins/{{mapper}}/{SELECTED_FASTAS}/{{contig_sample}}",
        out_dir="output/mag_qc/gtdbtk/{mapper}/{contig_sample}",
        db_path=config['params']['gtdbtk']['db_path'],
        # Two values because the two modes want different things. When
        # classify follows, GTDB-Tk's own default of 10 applies and the run
        # behaves as it always has. When it does not, the point is to have a
        # number for every genome, and 10 would drop the weakest ones into
        # filtered.tsv instead of reporting a low value for them -- and a
        # genome missing from the output is not the same as one measured at
        # 4%.
        min_perc_aa=(
            config['params']['gtdbtk'].get('min_perc_aa_no_classify', 0)
            if gtdbtk_mode() == 'align_only'
            else config['params']['gtdbtk'].get('min_perc_aa', 10)
        )
    threads:
        config['threads']['gtdbtk']
    retries:
        config['retries'].get('run_gtdbtk', 2)
    resources:
        mem_mb=mem_escalate('gtdbtk_align', base_default=32000, cap_multiple=2),
        runtime=runtime_escalate('gtdbtk_align', base_default=240)
    conda:
        "../env/gtdbtk.yaml"
    log:
        "output/logs/mag_qc/gtdbtk/{mapper}/{contig_sample}.align.log"
    benchmark:
        "output/benchmarks/mag_qc/gtdbtk/{mapper}/{contig_sample}.align.txt"
    shell:
        """
        if ! ls {params.bins_dir}/*.fa 2>/dev/null 1>/dev/null; then
            echo "No bins in {params.bins_dir}; skipping gtdbtk align." >> {log}
            exit 0
        fi

        existing=$(ls {params.out_dir}/align/*user_msa.fasta* 2>/dev/null | head -1 || true)
        if [ -n "$existing" ]; then
            echo "Reusing existing align output in {params.out_dir}/align" >> {log}
            touch -r "$existing" {output.done}
            exit 0
        fi

        export GTDBTK_DATA_PATH="{params.db_path}"
        gtdbtk align \
            --identify_dir {params.out_dir} \
            --out_dir {params.out_dir} \
            --cpus {threads} \
            --min_perc_aa {params.min_perc_aa} \
            2> {log} 1>&2
        """


rule run_gtdbtk:
    """
    Places each MAG in the GTDB reference tree.

    Keeps the name `run_gtdbtk` because the demon profile sets its wall
    clock and config sets retries.run_gtdbtk by that name, and this is the
    stage those values were measured for: pplacer is what takes hours and
    hundreds of gigabytes.

    Never requested under --skip-gtdbtk-classify. Nothing depends on its
    output in that mode, so it is not that the rule is skipped -- it is that
    the DAG never reaches it.
    """
    input:
        align="output/mag_qc/gtdbtk/{mapper}/{contig_sample}/.align.done"
    output:
        done=touch("output/mag_qc/gtdbtk/{mapper}/{contig_sample}/.done")
    params:
        bins_dir=f"output/selected_bins/{{mapper}}/{SELECTED_FASTAS}/{{contig_sample}}",
        out_dir="output/mag_qc/gtdbtk/{mapper}/{contig_sample}",
        db_path=config['params']['gtdbtk']['db_path']
    threads:
        config['threads']['gtdbtk']
    retries:
        config['retries'].get('run_gtdbtk', 2)
    resources:
        # pplacer's footprint scales with the number of bins, and asking for
        # the worst case up front would pin every sample to the few largest
        # nodes. Start at the value that suits most samples and escalate on
        # retry, capped by gtdbtk_max. The .get keeps configs predating
        # gtdbtk_max working.
        mem_mb=lambda wildcards, attempt: min(
            config['mem_mb'].get('gtdbtk_max', config['mem_mb']['gtdbtk']),
            config['mem_mb']['gtdbtk'] * attempt
        ),
        runtime=runtime_escalate('run_gtdbtk', base_default=1440)
    conda:
        "../env/gtdbtk.yaml"
    log:
        "output/logs/mag_qc/gtdbtk/{mapper}/{contig_sample}.log"
    benchmark:
        "output/benchmarks/mag_qc/gtdbtk/{mapper}/{contig_sample}.txt"
    shell:
        """
        if ! ls {params.bins_dir}/*.fa 2>/dev/null 1>/dev/null; then
            echo "No bins in {params.bins_dir}; skipping gtdbtk classify." >> {log}
            exit 0
        fi

        # Already classified by a classify_wf run predating the split.
        if ls {params.out_dir}/gtdbtk.*.summary.tsv 2>/dev/null 1>/dev/null; then
            echo "Reusing existing classification in {params.out_dir}" >> {log}
            exit 0
        fi

        export GTDBTK_DATA_PATH="{params.db_path}"
        gtdbtk classify \
            --genome_dir {params.bins_dir} \
            --align_dir {params.out_dir} \
            --out_dir {params.out_dir} \
            --cpus {threads} \
            --extension fa \
            2> {log} 1>&2
        """


rule make_mag_summary:
    """
    Collects all MAGs across samples, assigns globally sequential IDs (MAG_0001...N),
    looks up GTDB-tk taxonomy, CheckM2 and GUNC QC metrics, and writes mag_summary.tsv.
    Edit new_name in that table before running rename_mags to customize MAG names.
    """
    input:
        checkm2=lambda wildcards: expand(
            "output/mag_qc/checkm2/{mapper}/{contig_sample}/quality_report.tsv",
            mapper=config['mappers'],
            contig_sample=list(contig_pairings.keys())
        ),
        gunc=lambda wildcards: expand(
            "output/mag_qc/gunc/{mapper}/{contig_sample}/.done",
            mapper=config['mappers'],
            contig_sample=list(contig_pairings.keys())
        ),
        # Which GTDB-Tk stage the summary waits for. Requiring only the
        # stage that is actually wanted is what keeps classify out of the
        # DAG entirely under --skip-gtdbtk-classify, rather than running it
        # and discarding the result.
        #
        # Empty under --skip-gtdbtk, in which case the taxonomy and
        # MSA_Percent columns are written as NA rather than omitted, so a run
        # without GTDB-Tk stays distinguishable from a MAG it could not place.
        gtdbtk=lambda wildcards: ([] if gtdbtk_mode() == 'none' else expand(
            "output/mag_qc/gtdbtk/{mapper}/{contig_sample}/" + (
                ".align.done" if gtdbtk_mode() == 'align_only' else ".done"),
            mapper=config['mappers'],
            contig_sample=list(contig_pairings.keys())
        )),
        # The selection tool's own report. Which one depends on
        # params.consolidation.tool; only the active tool's outputs are
        # required, so switching tools does not demand the other's files.
        selection=lambda wildcards: expand(
            ("output/selected_bins/{mapper}/run_DAS_Tool/"
             "{contig_sample}_DASTool_summary.tsv")
            if consolidation_tool() == 'das_tool' else
            ("output/selected_bins/{mapper}/run_binette/{contig_sample}/"
             "final_bins_quality_reports.tsv"),
            mapper=config['mappers'],
            contig_sample=list(contig_pairings.keys())
        ),
        # Empty when --skip-cmseq was given, in which case the two strain
        # heterogeneity columns are written as NA rather than omitted, so a
        # run without CMSeq is distinguishable from a MAG it could not score.
        cmseq=lambda wildcards: (expand(
            "output/mag_qc/cmseq/{mapper}/{contig_sample}/strain_heterogeneity.tsv",
            mapper=config['mappers'],
            contig_sample=list(contig_pairings.keys())
        ) if cmseq_enabled() else [])
    output:
        summary="output/mag_qc/mag_summary.tsv"
    params:
        bins_base="output/selected_bins",
        gtdbtk_base="output/mag_qc/gtdbtk",
        checkm2_base="output/mag_qc/checkm2",
        gunc_base="output/mag_qc/gunc",
        cmseq_base="output/mag_qc/cmseq",
        mappers=config['mappers'],
        contig_samples=list(contig_pairings.keys()),
        assemblers=config['assemblers'],
        consolidation_tool=consolidation_tool(),
        selected_fastas=SELECTED_FASTAS
    retries:
        config['retries'].get('make_mag_summary', 2)
    resources:
        runtime=runtime_escalate('make_mag_summary', base_default=480)
    conda:
        "../env/mag_qc.yaml"
    log:
        "output/logs/mag_qc/make_mag_summary.log"
    script:
        "../scripts/make_mag_summary.py"


rule rename_mags:
    """
    Copies MAG FASTAs into output/mag_qc/renamed_mags/ using new_name from mag_summary.tsv.
    To rename MAGs: edit new_name in the table, then re-run: snakemake rename_mags
    """
    input:
        summary="output/mag_qc/mag_summary.tsv"
    output:
        done=touch("output/mag_qc/renamed_mags/.done")
    params:
        renamed_dir="output/mag_qc/renamed_mags"
    conda:
        "../env/mag_qc.yaml"
    log:
        "output/logs/mag_qc/rename_mags.log"
    script:
        "../scripts/rename_mags.py"
