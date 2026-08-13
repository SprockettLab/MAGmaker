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
    resources:
        mem_mb=config['mem_mb']['checkm2']
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


rule run_gtdbtk:
    input:
        done=f"output/selected_bins/{{mapper}}/{SELECTED_FASTAS}/{{contig_sample}}/.done"
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
        )
    conda:
        "../env/gtdbtk.yaml"
    log:
        "output/logs/mag_qc/gtdbtk/{mapper}/{contig_sample}.log"
    benchmark:
        "output/benchmarks/mag_qc/gtdbtk/{mapper}/{contig_sample}.txt"
    shell:
        """
        if ! ls {params.bins_dir}/*.fa 2>/dev/null 1>/dev/null; then
            echo "No bins in {params.bins_dir}; skipping gtdbtk." >> {log}
        else
            export GTDBTK_DATA_PATH="{params.db_path}"
            gtdbtk classify_wf \
                --genome_dir {params.bins_dir} \
                --out_dir {params.out_dir} \
                --cpus {threads} \
                --extension fa \
                2> {log} 1>&2
        fi
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
        gtdbtk=lambda wildcards: expand(
            "output/mag_qc/gtdbtk/{mapper}/{contig_sample}/.done",
            mapper=config['mappers'],
            contig_sample=list(contig_pairings.keys())
        ),
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
