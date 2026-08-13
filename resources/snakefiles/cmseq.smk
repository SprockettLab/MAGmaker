# Strain heterogeneity, per MAG, with CMSeq.
#
# CheckM2 detects material from another species and GUNC detects chimerism
# against a reference database. Neither sees a MAG that is a clean
# consensus of several co-resident strains of one species: everything in it
# belongs where it appears to. CMSeq measures that directly, by mapping the
# sample's own reads back to its own MAGs and counting positions carrying
# more than one allele.
#
# It matters at the tips of a phylogeny rather than for clade membership. A
# consensus of co-resident strains has a blended sequence, and blended tips
# blur exactly the short branches a co-diversification test depends on.
#
# Thresholds follow Pasolli 2019 and Sanders 2023 so the numbers are
# comparable with the published genome databases: minimum coverage 10, base
# quality above 30, and a MAG is only reported when more than 100 positions
# clear both.
#
# This cannot reuse the binning alignments. Those BAMs are temp() and
# deleted, and the all-vs-all maps PROTOTYPE reads to each sample's
# contigs, so a non-prototype sample's own reads are never mapped to its
# own contigs. One new mapping per sample is required and is what this
# costs.


def cmseq_enabled():
    """CMSeq runs unless --skip-cmseq was given."""
    return str(config.get('skip_cmseq', '')).strip().lower() not in (
        'true', '1', 'yes')


rule cmseq_reference:
    """One FASTA per sample holding that sample's selected MAGs.

    Mapping once against all of a sample's MAGs, rather than once per MAG,
    keeps this to a single alignment per sample. Reads are assigned to
    whichever MAG they fit best instead of being forced into each MAG in
    turn, which is also the more honest assignment."""
    input:
        done="output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}/.done"
    output:
        ref=temp("output/mag_qc/cmseq/{mapper}/{contig_sample}/mags.fa"),
        map=("output/mag_qc/cmseq/{mapper}/{contig_sample}/contig_to_mag.tsv")
    params:
        bins_dir="output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}"
    log:
        "output/logs/mag_qc/cmseq/{mapper}/{contig_sample}.reference.log"
    threads: 1
    shell:
        """
        : > {output.ref}
        : > {output.map}
        shopt -s nullglob
        for f in {params.bins_dir}/*.fa {params.bins_dir}/*.fna; do
            BIN=$(basename "$f"); BIN=${{BIN%.*}}
            cat "$f" >> {output.ref}
            grep '^>' "$f" | sed 's/^>//' | cut -d' ' -f1 \
                | awk -v b="$BIN" '{{print $1"\\t"b}}' >> {output.map}
        done
        echo "$(grep -c '^>' {output.ref}) contigs from \
$(cut -f2 {output.map} | sort -u | wc -l) MAGs" > {log}
        """


rule cmseq_map:
    """Map a sample's own host-filtered reads against its own MAGs."""
    input:
        ref="output/mag_qc/cmseq/{mapper}/{contig_sample}/mags.fa",
        reads=lambda wildcards: nonhost_reads(wildcards.contig_sample)
    output:
        bam=temp("output/mag_qc/cmseq/{mapper}/{contig_sample}/reads.bam"),
        bai=temp("output/mag_qc/cmseq/{mapper}/{contig_sample}/reads.bam.bai")
    threads:
        config['threads'].get('cmseq', 8)
    conda:
        "../env/cmseq.yaml"
    log:
        "output/logs/mag_qc/cmseq/{mapper}/{contig_sample}.map.log"
    benchmark:
        "output/benchmarks/mag_qc/cmseq/{mapper}/{contig_sample}.map.txt"
    resources:
        mem_mb=mem_escalate('cmseq_map', base_default=16000)
    shell:
        """
        # An empty reference means the sample produced no MAGs. Write a
        # valid empty BAM rather than failing the whole arm for it.
        if [ ! -s {input.ref} ]; then
            echo "no MAGs for this sample; empty BAM" > {log}
            samtools view -H -b /dev/null > {output.bam} 2>> {log} || \
                printf '' > {output.bam}
            touch {output.bai}
            exit 0
        fi
        MEM_PER_THREAD=$(( {resources.mem_mb} * 6 / 10 / {threads} ))
        minimap2 -ax sr -t {threads} {input.ref} {input.reads} 2> {log} \
            | samtools sort -m ${{MEM_PER_THREAD}}M -@ {threads} \
                -o {output.bam} - 2>> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        """


rule cmseq_poly:
    """CMSeq's per-reference polymorphic rate.

    Thresholds match Pasolli 2019 and Sanders 2023 so the values can be
    compared with the published databases. The raw output is kept rather
    than piped straight into the aggregator: it is small, and it is the
    only record of what CMSeq actually reported."""
    input:
        bam="output/mag_qc/cmseq/{mapper}/{contig_sample}/reads.bam",
        bai="output/mag_qc/cmseq/{mapper}/{contig_sample}/reads.bam.bai"
    output:
        raw="output/mag_qc/cmseq/{mapper}/{contig_sample}/poly_raw.tsv"
    params:
        mincov=config['params'].get('cmseq', {}).get('mincov', 10),
        minqual=config['params'].get('cmseq', {}).get('minqual', 30),
        dom=config['params'].get('cmseq', {}).get('dominant_frq_thrsh', 0.8)
    threads: 1
    conda:
        "../env/cmseq.yaml"
    log:
        "output/logs/mag_qc/cmseq/{mapper}/{contig_sample}.poly.log"
    shell:
        """
        if [ ! -s {input.bam} ]; then
            : > {output.raw}
            echo "empty BAM; nothing to evaluate" > {log}
            exit 0
        fi
        poly.py -f \
            --mincov {params.mincov} \
            --minqual {params.minqual} \
            --dominant_frq_thrsh {params.dom} \
            {input.bam} > {output.raw} 2> {log} || {{
                echo "poly.py failed; leaving an empty table" >> {log}
                : > {output.raw}
            }}
        """


rule cmseq_aggregate:
    """Roll CMSeq's per-contig output up to one row per MAG."""
    input:
        raw="output/mag_qc/cmseq/{mapper}/{contig_sample}/poly_raw.tsv",
        map="output/mag_qc/cmseq/{mapper}/{contig_sample}/contig_to_mag.tsv"
    output:
        tsv="output/mag_qc/cmseq/{mapper}/{contig_sample}/strain_heterogeneity.tsv"
    params:
        min_positions=config['params'].get('cmseq', {}).get(
            'min_positions', 100)
    log:
        "output/logs/mag_qc/cmseq/{mapper}/{contig_sample}.aggregate.log"
    threads: 1
    conda:
        "../env/cmseq.yaml"
    script:
        "../scripts/aggregate_cmseq.py"
