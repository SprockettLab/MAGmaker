# Viral and plasmid discovery track.
#
# geNomad classifies viral and plasmid contigs directly from the assemblies
# (no read mapping needed), then CheckV scores viral completeness/contamination.
# This runs parallel to the MAG track and does not touch it.
#
# Scope note: this is a PER-SAMPLE track. Dereplicating viruses into vOTUs
# across samples (95% ANI / 85% AF, MIUViG) is a downstream step, kept out of
# MAGmaker for the same reason dRep is -- it is a cross-sample analysis choice,
# not part of per-sample genome recovery.


rule run_genomad:
    """
    geNomad end-to-end: identify viral + plasmid contigs and classify them.
    """
    input:
        contigs = "output/assemble/{assembler}/{contig_sample}.contigs.fasta"
    output:
        done = touch("output/virus/genomad/{assembler}/{contig_sample}/.done")
    params:
        outdir = "output/virus/genomad/{assembler}/{contig_sample}",
        db = config['params']['genomad']['db_path'],
        extra = config['params']['genomad']['extra']
    threads:
        config['threads']['genomad']
    retries:
        config['retries'].get('run_genomad', 2)
    resources:
        mem_mb=config['mem_mb'].get('run_genomad', 32000),
        runtime=runtime_escalate('run_genomad', base_default=480)
    conda:
        "../env/genomad.yaml"
    benchmark:
        "output/benchmarks/virus/genomad/{assembler}/{contig_sample}.txt"
    log:
        "output/logs/virus/genomad/{assembler}/{contig_sample}.log"
    shell:
        # NB: geNomad's own --cleanup is deliberately NOT used. It calls
        # shutil.rmtree() on its mmseqs2 scratch dir, which races NFS
        # silly-rename on isilon: files deleted while a handle is still open
        # become .nfs* entries, so rmdir fails with OSError ENOTEMPTY
        # ([Errno 39] Directory not empty: .../besthit_db) and the whole job
        # exits non-zero AFTER all real outputs are already written. Instead we
        # let geNomad finish, then remove the heavy module intermediates
        # ourselves once the process has fully exited (no open handles), keeping
        # the *_summary dir (+ virus/plasmid .fna) that downstream rules read.
        # The rm is best-effort (|| true) so a stray .nfs* file can't fail the
        # job; NFS reaps those once their holder is gone.
        #
        # The output directory is cleared FIRST, on every attempt. geNomad
        # resumes from whatever intermediates it finds, so a retry after a
        # killed attempt mixes new results with a half-written table. On the
        # WF18 Jasinska arm that produced
        #   summary.py:52 KeyError: np.str_('k119_...|provirus_...')
        # on 7 of 90 samples -- a provirus present in the virus list but
        # missing from the length lookup built by an earlier attempt. Deleting
        # the directory and rerunning fixed all 7 with no other change.
        # A retry therefore redoes the work instead of resuming it. That is
        # the intended trade: a slower retry beats a wrong table.
        """
            rm -rf {params.outdir}
            mkdir -p {params.outdir}
            genomad end-to-end --threads {threads} {params.extra} \
                {input.contigs} {params.outdir} {params.db} \
                2> {log} 1>&2
            rm -rf {params.outdir}/*_annotate \
                   {params.outdir}/*_aggregated_classification \
                   {params.outdir}/*_find_proviruses \
                   {params.outdir}/*_marker_classification \
                   {params.outdir}/*_nn_classification \
                   2> /dev/null || true
        """


rule run_checkv:
    """
    CheckV quality (completeness/contamination) for the geNomad viral contigs.
    Skips cleanly when a sample has no viral contigs.
    """
    input:
        "output/virus/genomad/{assembler}/{contig_sample}/.done"
    output:
        done = touch("output/virus/checkv/{assembler}/{contig_sample}/.done")
    params:
        virus_fna = ("output/virus/genomad/{assembler}/{contig_sample}/"
                     "{contig_sample}.contigs_summary/"
                     "{contig_sample}.contigs_virus.fna"),
        outdir = "output/virus/checkv/{assembler}/{contig_sample}",
        db = config['params']['checkv']['db_path']
    threads:
        config['threads']['checkv']
    retries:
        config['retries'].get('run_checkv', 2)
    resources:
        mem_mb=config['mem_mb'].get('run_checkv', 16000),
        runtime=runtime_escalate('run_checkv', base_default=240)
    conda:
        "../env/checkv.yaml"
    benchmark:
        "output/benchmarks/virus/checkv/{assembler}/{contig_sample}.txt"
    log:
        "output/logs/virus/checkv/{assembler}/{contig_sample}.log"
    shell:
        # The output directory is cleared FIRST, on every attempt, as a
        # precaution matching run_genomad above. CheckV also reuses its tmp/
        # directory, so a retry can in principle resume from a half-written
        # state, and unlike geNomad it would do so without failing.
        #
        # What was OBSERVED on the WF18 Jasinska 2013 megahit arm: 9 of 90
        # samples produced a quality_summary.tsv with 274 to 569 rows against
        # 3,504 to 7,335 geNomad virus sequences, and exited 0. That left
        # 44,170 of 380,306 viral contigs (11.6%) with an empty
        # checkv_quality field in the merged summary, with nothing in the
        # table to show it was incomplete.
        #
        # What was NOT established: why. A leftover tmp/ is NOT the signal --
        # every sample in every arm has one, because CheckV never cleans it
        # up. The three other arms of the same build show no truncation at
        # all. So this clear-first is a safety measure whose necessity for
        # CheckV is unproven, not a fix for a diagnosed cause.
        """
            rm -rf {params.outdir}
            mkdir -p {params.outdir}
            if [ -s {params.virus_fna} ]; then
                checkv end_to_end {params.virus_fna} {params.outdir} \
                    -t {threads} -d {params.db} 2> {log} 1>&2
            else
                echo "No viral contigs for this sample; skipping CheckV." > {log}
            fi
        """


rule make_virus_summary:
    """
    Merge per-sample geNomad (virus + plasmid) and CheckV outputs into one table.
    """
    input:
        genomad = lambda wildcards: expand(
            "output/virus/genomad/{assembler}/{contig_sample}/.done",
            assembler=config['assemblers'],
            contig_sample=contig_pairings.keys()),
        checkv = lambda wildcards: expand(
            "output/virus/checkv/{assembler}/{contig_sample}/.done",
            assembler=config['assemblers'],
            contig_sample=contig_pairings.keys())
    output:
        summary = "output/virus/virus_summary.tsv"
    params:
        contig_samples = list(contig_pairings.keys()),
        assemblers = config['assemblers'],
        genomad_base = "output/virus/genomad",
        checkv_base = "output/virus/checkv"
    conda:
        "../env/mag_qc.yaml"
    log:
        "output/logs/virus/make_virus_summary.log"
    script:
        "../scripts/make_virus_summary.py"


rule virus_all:
    input:
        "output/virus/virus_summary.tsv"
