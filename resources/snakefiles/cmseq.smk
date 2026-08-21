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


rule cmseq_reference:
    """One FASTA per sample holding that sample's selected MAGs.

    Mapping once against all of a sample's MAGs, rather than once per MAG,
    keeps this to a single alignment per sample. Reads are assigned to
    whichever MAG they fit best instead of being forced into each MAG in
    turn, which is also the more honest assignment."""
    input:
        done=f"output/selected_bins/{{mapper}}/{SELECTED_FASTAS}/{{contig_sample}}/.done"
    output:
        ref=temp("output/mag_qc/cmseq/{mapper}/{contig_sample}/mags.fa"),
        map=("output/mag_qc/cmseq/{mapper}/{contig_sample}/contig_to_mag.tsv")
    params:
        bins_dir=f"output/selected_bins/{{mapper}}/{SELECTED_FASTAS}/{{contig_sample}}"
    log:
        "output/logs/mag_qc/cmseq/{mapper}/{contig_sample}.reference.log"
    threads: 1
    retries:
        config['retries'].get('cmseq_reference', 2)
    resources:
        runtime=runtime_escalate('cmseq_reference', base_default=240)
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
    params:
        # Overridable from the command line as well as the config file, the
        # same way consolidation_tool is, so the effect of the filters can
        # be measured by running one arm twice without editing a tracked
        # file between the two. Setting both to 0 restores the unfiltered
        # behaviour these replaced.
        min_mapq=config.get(
            'cmseq_min_mapq',
            config['params']['cmseq'].get('min_mapq', 20)),
        min_identity=config.get(
            'cmseq_min_read_identity',
            config['params']['cmseq'].get('min_read_identity', 0.95)),
        # The mismatch budget as a fraction of read length, computed HERE in
        # Python rather than in the shell. The shell version called python3,
        # which is not guaranteed to be in mapping.yaml -- that environment
        # exists to provide minimap2 and samtools. When it is absent the
        # substitution yields an empty string and samtools receives
        # `[NM] <=  * qlen`, which is a filter-expression syntax error rather
        # than an unfiltered run, so it fails loudly, but it fails for a
        # reason that has nothing to do with the data.
        max_nm_frac=lambda wildcards: 1.0 - float(config.get(
            'cmseq_min_read_identity',
            config['params']['cmseq'].get('min_read_identity', 0.95)))
    threads:
        config['threads'].get('cmseq', 8)
    conda:
        # mapping.yaml already provides minimap2 and samtools, and is built
        # in every arm. A dedicated environment for this rule would be a
        # fourth copy of both and one more environment to create.
        "../env/mapping.yaml"
    log:
        "output/logs/mag_qc/cmseq/{mapper}/{contig_sample}.map.log"
    benchmark:
        "output/benchmarks/mag_qc/cmseq/{mapper}/{contig_sample}.map.txt"
    retries:
        config['retries'].get('cmseq_map', 2)
    resources:
        mem_mb=mem_escalate('cmseq_map', base_default=16000),
        runtime=runtime_escalate('cmseq_map', base_default=480)
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

        # Two filters, both aimed at reads that are not from the organism
        # whose MAG they land on. Without them, strain heterogeneity partly
        # measures how much of the community failed to bin: every read in
        # the sample is mapped to that sample's MAGs, so reads from species
        # that were never binned pile onto the nearest available reference
        # and read out as polymorphism. Measured across five arms before
        # this was added, strain heterogeneity ranked exactly inversely to
        # how well each arm's MAGs were covered, and every within-arm rank
        # correlation was negative, which is that mechanism's signature.
        #
        # -q {params.min_mapq} drops reads that map about as well somewhere
        # else in the reference. The reference is a concatenation of the
        # sample's own MAGs, so this is what removes ambiguity BETWEEN this
        # sample's genomes.
        #
        # The NM expression keeps only reads aligning at roughly
        # {params.min_identity} identity or better. That is what removes
        # reads from relatives with no MAG of their own, which map uniquely
        # -- so MAPQ cannot catch them -- but divergently. Within-species
        # strain variation sits above 95% by definition, so a filter there
        # removes cross-species reads while keeping the signal being
        # measured.
        #
        # -F 0x904 drops unmapped, secondary and supplementary alignments.
        # Unmapped reads are excluded FIRST because they carry no NM tag and
        # the identity expression below has nothing to evaluate on them;
        # secondary and supplementary would otherwise let one read vote at
        # several positions.
        #
        # The identity filter is applied with awk rather than samtools'
        # -e filter expression. mapping.yaml pins tbb, libzlib and minimap2
        # to old versions, which drags conda down to samtools 1.6 or 1.9 --
        # and -e did not exist until 1.12. An older samtools does not reject
        # the option, it treats the expression as a positional FILENAME and
        # dies with "failed to open ... for reading", which points nowhere
        # near the real cause. Loosening the pins to get a newer samtools
        # would rebuild an environment used by the main mapping step across
        # every arm, and concurrent rebuilds of a shared environment have
        # already corrupted the package cache on this cluster once.
        #
        # awk needs no version of anything. Flag and MAPQ filtering stay in
        # samtools, where they work in every release.
        #
        # A read with no NM tag is dropped, since its identity cannot be
        # judged. After -F 0x904 every surviving record is a primary
        # alignment, and minimap2 emits NM for those, so in practice this
        # removes nothing.
        minimap2 -ax sr -t {threads} {input.ref} {input.reads} 2> {log} \
            | samtools view -h -F 0x904 -q {params.min_mapq} - 2>> {log} \
            | awk -v maxfrac={params.max_nm_frac} '
                /^@/ {{ print; next }}
                {{
                    nm = -1
                    for (i = 12; i <= NF; i++)
                        if ($i ~ /^NM:i:/) {{ split($i, a, ":"); nm = a[3] + 0; break }}
                    if (nm >= 0 && nm <= maxfrac * length($10)) print
                }}' FS='\t' OFS='\t' \
            | samtools sort -m ${{MEM_PER_THREAD}}M -@ {threads} \
                -o {output.bam} - 2>> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}

        # Report what the filters cost, so a run that removed most of its
        # reads is visible rather than silently producing clean-looking
        # genomes. A MAG whose reads were nearly all filtered will report
        # low heterogeneity for want of evidence, not for want of strains.
        echo "reads retained after filtering: $(samtools view -c {output.bam} 2>/dev/null || echo 0)" >> {log}
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
    retries:
        config['retries'].get('cmseq_poly', 2)
    resources:
        runtime=runtime_escalate('cmseq_poly', base_default=1440)
    conda:
        # metaphlan depends on cmseq, so profile.yaml already provides
        # poly.py and that environment is built in every arm. If metaphlan
        # ever drops that dependency, add `cmseq` to profile.yaml
        # explicitly rather than creating a separate environment: building
        # a new one across many arms at once corrupts the shared conda
        # package cache, which is how this was discovered.
        "../env/profile.yaml"
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
    retries:
        config['retries'].get('cmseq_aggregate', 2)
    resources:
        runtime=runtime_escalate('cmseq_aggregate', base_default=240)
    conda:
        # the script is standard library only
        "../env/mag_qc.yaml"
    script:
        "../scripts/aggregate_cmseq.py"
