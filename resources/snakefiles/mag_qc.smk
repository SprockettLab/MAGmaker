import os
import glob as _glob
from os.path import join, dirname, basename


rule run_checkm2:
    input:
        done="output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}/.done"
    output:
        report="output/mag_qc/checkm2/{mapper}/{contig_sample}/quality_report.tsv"
    params:
        bins_dir="output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}",
        out_dir="output/mag_qc/checkm2/{mapper}/{contig_sample}",
        db_path=config['params']['checkm2']['db_path']
    threads:
        config['threads']['checkm2']
    resources:
        mem_mb=config['mem_mb']['checkm2']
    conda:
        "../env/mag_qc.yaml"
    log:
        "output/logs/mag_qc/checkm2/{mapper}/{contig_sample}.log"
    benchmark:
        "output/benchmarks/mag_qc/checkm2/{mapper}/{contig_sample}.txt"
    shell:
        """
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
        """


rule run_gunc:
    input:
        done="output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}/.done"
    output:
        done=touch("output/mag_qc/gunc/{mapper}/{contig_sample}/.done")
    params:
        bins_dir="output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}",
        out_dir="output/mag_qc/gunc/{mapper}/{contig_sample}",
        db_flag=lambda wildcards: (
            f"--db_file {config['params']['gunc']['db_path']}"
            if config['params']['gunc']['db_path']
            else ""
        )
    threads:
        config['threads']['gunc']
    conda:
        "../env/mag_qc.yaml"
    log:
        "output/logs/mag_qc/gunc/{mapper}/{contig_sample}.log"
    benchmark:
        "output/benchmarks/mag_qc/gunc/{mapper}/{contig_sample}.txt"
    shell:
        """
        gunc run \
            --input_dir {params.bins_dir} \
            --out_dir {params.out_dir} \
            --threads {threads} \
            --file_suffix .fa \
            {params.db_flag} \
            2> {log} 1>&2
        """


rule run_gtdbtk:
    input:
        done="output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}/.done"
    output:
        done=touch("output/mag_qc/gtdbtk/{mapper}/{contig_sample}/.done")
    params:
        bins_dir="output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}",
        out_dir="output/mag_qc/gtdbtk/{mapper}/{contig_sample}",
        db_path=config['params']['gtdbtk']['db_path']
    threads:
        config['threads']['gtdbtk']
    resources:
        mem_mb=config['mem_mb']['gtdbtk']
    conda:
        "../env/gtdbtk.yaml"
    log:
        "output/logs/mag_qc/gtdbtk/{mapper}/{contig_sample}.log"
    benchmark:
        "output/benchmarks/mag_qc/gtdbtk/{mapper}/{contig_sample}.txt"
    shell:
        """
        gtdbtk classify_wf \
            --genome_dir {params.bins_dir} \
            --out_dir {params.out_dir} \
            --cpus {threads} \
            --extension fa \
            --skip_ani_screen \
            --data_dir {params.db_path} \
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
        gtdbtk=lambda wildcards: expand(
            "output/mag_qc/gtdbtk/{mapper}/{contig_sample}/.done",
            mapper=config['mappers'],
            contig_sample=list(contig_pairings.keys())
        ),
        das_tool=lambda wildcards: expand(
            "output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}_DASTool_summary.tsv",
            mapper=config['mappers'],
            contig_sample=list(contig_pairings.keys())
        )
    output:
        summary="output/mag_qc/mag_summary.tsv"
    params:
        bins_base="output/selected_bins",
        gtdbtk_base="output/mag_qc/gtdbtk",
        checkm2_base="output/mag_qc/checkm2",
        gunc_base="output/mag_qc/gunc",
        mappers=config['mappers'],
        contig_samples=list(contig_pairings.keys())
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
