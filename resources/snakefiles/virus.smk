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
    conda:
        "../env/genomad.yaml"
    benchmark:
        "output/benchmarks/virus/genomad/{assembler}/{contig_sample}.txt"
    log:
        "output/logs/virus/genomad/{assembler}/{contig_sample}.log"
    shell:
        """
            genomad end-to-end --cleanup --threads {threads} {params.extra} \
                {input.contigs} {params.outdir} {params.db} \
                2> {log} 1>&2
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
    conda:
        "../env/checkv.yaml"
    benchmark:
        "output/benchmarks/virus/checkv/{assembler}/{contig_sample}.txt"
    log:
        "output/logs/virus/checkv/{assembler}/{contig_sample}.log"
    shell:
        """
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
