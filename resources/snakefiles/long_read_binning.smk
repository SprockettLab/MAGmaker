# =============================================================================
# long_read_binning.smk -- what a long-read binning row adds (Snakefile-bin)
#
# A long-read assembly (<sample>__flye, <sample>__myloasm) is binned like every
# other assembly, with the prototypes' short-read coverage. This adds:
#
#   map_long_reads_to_contigs     the sample's own long reads mapped to it: one
#                                 more coverage column for MetaBAT2, CONCOCT and
#                                 SemiBin2 (get_bam_list in binning.smk)
#   make_maxbin2_long_read_coverage  the same coverage as a MaxBin2 table
#   circular_candidates           one bin per circular contig of at least
#                                 params.long_reads.circular_bin_min_length,
#                                 offered to DAS_Tool/Binette beside the
#                                 binners' bins (binner_tables in common.smk)
#
# Nothing is removed before binning: every contig, circular or not, goes to
# every binner. A circular contig offered as its own bin is kept only if it
# scores better than the overlapping candidates and passes the consolidation
# tool's thresholds, like any other candidate.
# =============================================================================


rule map_long_reads_to_contigs:
    """
    minimap2 -ax map-ont of the sample's host-filtered long reads to the row's
    contigs, sorted and indexed. Every read is used.
    """
    input:
        contigs=lambda wildcards: row_contigs(wildcards.contig_sample),
        reads=lambda wildcards: "output/qc/host_filter_ont/nonhost/{0}.ont.fastq.gz".format(
            sample_of(wildcards.contig_sample))
    output:
        bam="output/mapping/long_reads/{contig_sample}.bam",
        bai="output/mapping/long_reads/{contig_sample}.bam.bai"
    conda:
        "../env/long_read_qc.yaml"
    threads:
        config['threads'].get('map_long_reads', 16)
    retries:
        config['retries'].get('map_long_reads', 2)
    resources:
        mem_mb=mem_escalate('map_long_reads', base_default=32000),
        runtime=runtime_escalate('map_long_reads', base_default=720)
    log:
        "output/logs/mapping/long_reads/{contig_sample}.log"
    benchmark:
        "output/benchmarks/mapping/long_reads/{contig_sample}_benchmark.txt"
    shell:
        """
        minimap2 -ax map-ont -t {threads} {input.contigs} {input.reads} 2> {log} \
          | samtools sort -@ 4 -m 2G -o {output.bam}.partial - 2>> {log}
        mv {output.bam}.partial {output.bam}
        samtools index {output.bam} 2>> {log}
        """


rule make_maxbin2_long_read_coverage:
    """MaxBin2's per-contig mean depth from the long-read BAM, as for short reads."""
    input:
        bam="output/mapping/long_reads/{contig_sample}.bam"
    output:
        coverage_table="output/binning/maxbin2/long_read_coverage/{contig_sample}_coverage.txt"
    conda:
        "../env/binning.yaml"
    log:
        "output/logs/binning/maxbin2/long_read_coverage/{contig_sample}.log"
    shell:
        """
        samtools coverage {input.bam} | tail -n +2 | sort -k1 | cut -f1,6 > {output.coverage_table} 2> {log}
        """


rule circular_candidates:
    """
    Each contig the assembler called circular ("yes"; myloasm's "possibly" is
    not counted) and at least circular_bin_min_length long, still present
    after DNA Control Strand removal, as its own bin. Written as a
    scaffolds2bin table (contig, bin), empty when there are none; empty tables
    are skipped by DAS_Tool and Binette like a binner that declined a sample.
    """
    input:
        info=lambda wildcards: "output/assemble/{0}/{1}.contig_info.tsv".format(
            assembler_of(wildcards.contig_sample), sample_of(wildcards.contig_sample)),
        contigs=lambda wildcards: row_contigs(wildcards.contig_sample)
    output:
        table="output/selected_bins/circular/{mapper}/scaffolds2bin/{contig_sample}_scaffolds2bin.tsv"
    params:
        min_length=config.get('params', {}).get('long_reads', {}).get('circular_bin_min_length', 500000)
    log:
        "output/logs/selected_bins/circular/{mapper}/{contig_sample}.log"
    shell:
        """
        grep '^>' {input.contigs} | cut -c2- | awk '{{print $1}}' > {output.table}.present
        awk -F'\\t' -v m={params.min_length} -v row={wildcards.contig_sample} '
            FILENAME == ARGV[1] {{present[$1] = 1; next}}
            FNR > 1 && $4 == "yes" && $3 >= m && ($1 in present) {{n++; print $1 "\\t" row "_circular_" n}}
            ' {output.table}.present {input.info} > {output.table}
        rm -f {output.table}.present
        echo "$(wc -l < {output.table}) circular contigs >= {params.min_length} bp offered as bins" > {log}
        """
