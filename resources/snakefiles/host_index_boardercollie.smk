# Host index, BoarderCollie on. Included by host_reference.smk, which sets
# host_reference, host_base, BOARDERCOLLIE_SNAKEFILE and BOARDERCOLLIE_CONFIG.
#
# BoarderCollie's steps are loaded with a boardercollie_ prefix. Asking for
# host_reference runs everything needed to make it; nothing else of
# BoarderCollie's is requested (its own bowtie2 and minimap2 indexes are not
# built).

module boardercollie:
    snakefile: BOARDERCOLLIE_SNAKEFILE
    config: BOARDERCOLLIE_CONFIG

use rule * from boardercollie as boardercollie_*


rule host_bowtie2_build:
    input:
        # ancient: the reference is a hard link of BoarderCollie's shared
        # library FASTA, and every other run folder with this host touches
        # that same file (Snakemake touches outputs after a job). Without
        # ancient, any such run made this folder's index look stale, and the
        # rerun of this rule would have rerun host_filter for every sample. The index's md5 check below still catches a changed
        # reference, and a rerun of reference_out here still reruns this.
        reference=ancient(host_reference),
        members=host_reference + ".members.tsv"
    output:
        # Holds the reference's md5 and the index stem, so the record of
        # which reference an index was built from sits beside it.
        host_base + ".bowtie2_build.done"
    log:
        "output/logs/qc/host_bowtie2_build/host_bowtie2_build.log"
    benchmark:
        "output/benchmarks/qc/host_bowtie2_build/host_bowtie2_build_benchmark.txt"
    conda:
        "../env/qc.yaml"
    params:
        indexbase=host_base
    threads:
        config['threads']['host_filter']
    retries:
        config['retries'].get('host_bowtie2_build', 2)
    resources:
        # Memory scales with the reference: the profile's 8 GB default held a
        # build of the ~8.8 Gbp three-species vervet reference at its limit
        # (MaxRSS 8,190,112 KB of 8,192,000) with almost no CPU progress.
        mem_mb=mem_escalate('host_bowtie2_build', base_default=64000),
        runtime=runtime_escalate('host_bowtie2_build', base_default=1440)
    # DECISION: the index is tied to the reference's contents through its md5.
    # BoarderCollie writes the joined reference to the same path whatever its
    # contents (other genomes, other settings), so a complete index at the
    # same stem says nothing about whether it was built from this reference.
    # The md5 is written only after a build finishes, and removed before one
    # starts, so an interrupted build is never taken for a finished one.
    shell:
        """
        REF="{input.reference}"
        IDX="{params.indexbase}"
        STAMP="$IDX.reference.md5"
        if command -v md5sum >/dev/null 2>&1; then
            MD5=$(md5sum "$REF" | cut -d' ' -f1)
        else
            MD5=$(md5 -q "$REF")
        fi

        FOUND=0
        for suffix in bt2 bt2l; do
            COMPLETE=1
            for part in 1 2 3 4 rev.1 rev.2; do
                [ -f "$IDX.$part.$suffix" ] || COMPLETE=0
            done
            if [ "$COMPLETE" -eq 1 ]; then FOUND=1; fi
        done

        if [ "$FOUND" -eq 1 ] && [ -f "$STAMP" ] && [ "$(cat "$STAMP")" = "$MD5" ]; then
            echo "Complete bowtie2 index of $REF (md5 $MD5) already at $IDX; reused." > {log}
        else
            echo "Building the bowtie2 index of $REF (md5 $MD5) at $IDX." > {log}
            rm -f "$STAMP"
            for suffix in bt2 bt2l; do
                for part in 1 2 3 4 rev.1 rev.2; do
                    rm -f "$IDX.$part.$suffix"
                done
            done
            bowtie2-build --threads {threads} "$REF" "$IDX" >> {log} 2>&1
            echo "$MD5" > "$STAMP"
        fi

        cp {input.members} "$IDX.members.tsv"
        printf 'reference\\t%s\\nmd5\\t%s\\nindex\\t%s\\n' "$REF" "$MD5" "$IDX" > {output}
        """
