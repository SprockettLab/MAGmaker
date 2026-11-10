# Host index, BoarderCollie off: host_filter.genome is an existing FASTA used
# as it is, and its bowtie2 index lives in host_filter.db_dir under the
# FASTA's file name stem. Included by host_reference.smk.
#
# This rule is kept exactly as it was before BoarderCollie was added.
# Snakemake reruns a job whose shell text has changed, and a rerun here would
# make every host_filter job of a finished project run again.

rule host_bowtie2_build:
    input:
        reference=config['host_filter']['genome']
    output:
        # bowtie2-build writes small-index (.bt2) files for references
        # under ~4 Gbp and large-index (.bt2l) files above it. Declaring
        # fixed .bt2 outputs breaks on large hosts (e.g. multi-genome
        # vervet references), where only .bt2l files appear. Track a
        # sentinel instead; `bowtie2 -x` auto-detects the index layout at
        # align time.
        touch(host_base + ".bowtie2_build.done")
    log:
        "output/logs/qc/host_bowtie2_build/host_bowtie2_build.log"
    benchmark:
        "output/benchmarks/qc/host_bowtie2_build/host_bowtie2_build_benchmark.txt"
    conda:
        "../env/qc.yaml"
    params:
        extra="",
        indexbase=host_base
    threads:
        config['threads']['host_filter']
    retries:
        config['retries'].get('host_bowtie2_build', 2)
    resources:
        runtime=runtime_escalate('host_bowtie2_build', base_default=1440)
    shell:
        # A complete index next to the reference is reused rather than
        # rebuilt. Shared host indexes commonly predate the sentinel, and
        # rebuilding a human genome costs hours -- silently, since the only
        # symptom is host_bowtie2_build appearing in the job list.
        #
        # All six parts must be present before the build is skipped: a
        # partial index (an interrupted build, an incomplete copy) would
        # otherwise be accepted and fail later inside bowtie2, which is a
        # much harder failure to trace back to its cause. Small (.bt2) and
        # large (.bt2l) layouts are checked separately, since a reference
        # produces one or the other, never a mixture.
        """
        IDX="{params.indexbase}"
        FOUND=0
        for suffix in bt2 bt2l; do
            COMPLETE=1
            for part in 1 2 3 4 rev.1 rev.2; do
                [ -f "$IDX.$part.$suffix" ] || COMPLETE=0
            done
            if [ "$COMPLETE" -eq 1 ]; then FOUND=1; fi
        done

        if [ "$FOUND" -eq 1 ]; then
            echo "Complete bowtie2 index already present at $IDX" > {log}
            echo "Skipping bowtie2-build; touching sentinel instead." >> {log}
        else
            bowtie2-build --threads {threads} {params.extra} \
            {input.reference} "$IDX" 2> {log} 1>&2
        fi
        """
