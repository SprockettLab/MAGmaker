# vRhyme viral binning.
#
# vRhyme is a BINNER, not a predictor. It joins geNomad's viral contigs into
# viral bins (vMAGs), so a phage that assembled in several pieces is reported
# as one genome instead of several partial ones. It computes nothing new: it
# reads geNomad's viral FASTA and the MetaBAT2 coverage table the MAG track
# already wrote, so it needs geNomad on and metabat2 among the binners.
# virus.smk's VRHYME_RUNS checks both and says so once when it cannot run.
#
# Included by virus.smk, which holds the switches and the path helpers.

rule vrhyme_prep:
    """Choose which viral contigs to bin and subset the coverage table.

    Selection, and why each part:
      Element_Type == virus  -- plasmids are a different question.
      not a provirus, by geNomad's name OR CheckV's call -- a dormant
        prophage's coverage tracks its host, so vRhyme would bin it with
        whatever host contigs are nearby. geNomad excises proviruses and
        writes |provirus_<start>_<end> into the name; CheckV then reports
        provirus=No for exactly those, because no host DNA is left, so the
        CheckV column alone misses most of them.
      length >= params.vrhyme.min_length (2000, vRhyme's floor).
      every CheckV quality tier, Low-quality and Not-determined included --
        the pieces of a fragmented phage each score badly alone, and joining
        them is the point. Quality is judged after binning, on the vMAGs.
    """
    input:
        genomad="output/virus/genomad/{assembler}/{contig_sample}/.done",
        checkv="output/virus/checkv/{assembler}/{contig_sample}/.done",
        coverage=lambda w: (
            "output/binning/metabat2/{0}/coverage_tables/{1}_coverage_table.txt".format(
                vrhyme_mapper(), w.contig_sample)),
    output:
        fasta="output/virus/vrhyme/{assembler}/{contig_sample}/input.viral.fna",
        coverage="output/virus/vrhyme/{assembler}/{contig_sample}/input.coverage.jgi.tsv",
        params_file="output/virus/vrhyme/{assembler}/{contig_sample}/prep_params.txt",
    params:
        virus_fna=lambda w: predictor_fasta('genomad', w.assembler, w.contig_sample),
        checkv_dir=lambda w: checkv_dir('genomad', w.assembler, w.contig_sample),
        min_length=config['params'].get('vrhyme', {}).get('min_length', 2000),
    conda:
        "../env/mag_qc.yaml"
    log:
        "output/logs/virus/vrhyme/{assembler}/{contig_sample}.prep.log"
    script:
        "../scripts/vrhyme_prep.py"


rule run_vrhyme:
    """Bin one assembly's viral contigs into vMAGs."""
    input:
        fasta=rules.vrhyme_prep.output.fasta,
        coverage=rules.vrhyme_prep.output.coverage,
    output:
        done=touch("output/virus/vrhyme/{assembler}/{contig_sample}/.done"),
    params:
        outdir="output/virus/vrhyme/{assembler}/{contig_sample}/vrhyme_results",
        converted="output/virus/vrhyme/{assembler}/{contig_sample}/input.coverage.vrhyme.tsv",
        min_length=config['params'].get('vrhyme', {}).get('min_length', 2000),
    threads: config['threads'].get('vrhyme', 16)
    retries: config['retries'].get('run_vrhyme', 2)
    resources:
        mem_mb=mem_escalate('run_vrhyme', base_default=32000),
        runtime=runtime_escalate('run_vrhyme', base_default=720),
    conda:
        tool_env('vrhyme') if VRHYME_RUNS else "../env/mag_qc.yaml"
    benchmark:
        "output/benchmarks/virus/vrhyme/{assembler}/{contig_sample}.txt"
    log:
        "output/logs/virus/vrhyme/{assembler}/{contig_sample}.log"
    shell:
        # The coverage table is converted by vRhyme's own bundled
        # aux/coverage_table_convert.py, found inside the environment, so the
        # conversion is always the one that vRhyme version expects.
        #
        # Flags: only -i, -c, -o, -t and -l. --outliers is incompatible with
        # -c coverage-table input, upstream says.
        #
        # vRhyme is NOT deterministic: it keeps the best of 20 iterations and
        # has no seed option, so an identical input can give a slightly
        # different answer. The result is what this run produced.
        #
        # An empty input means geNomad called nothing bin-worthy for this
        # assembly, which is an outcome, not a failure.
        """
        rm -rf {params.outdir}
        n_in=$(grep -c '^>' {input.fasta} || true)
        if [ "$n_in" -eq 0 ]; then
            echo "No viral contigs to bin for this assembly." > {log}
            mkdir -p {params.outdir}
            exit 0
        fi
        convert=$(find "$CONDA_PREFIX" -path '*/vRhyme/aux/coverage_table_convert.py' \
            -print -quit)
        if [ -z "$convert" ]; then
            echo "ERROR: vRhyme's coverage_table_convert.py not found in $CONDA_PREFIX" > {log}
            exit 1
        fi
        python "$convert" -i {input.coverage} -o {params.converted} > {log} 2>&1
        vRhyme -i {input.fasta} -c {params.converted} -o {params.outdir} \
            -t {threads} -l {params.min_length} >> {log} 2>&1
        """


rule vrhyme_summary:
    """One row per vMAG across all assemblies, plus a membership table."""
    input:
        done=lambda w: ["output/virus/vrhyme/{0}/{1}/.done".format(a, r)
                        for a, r in vrhyme_pairs()],
    output:
        summary="output/virus/vrhyme_summary.tsv",
        membership="output/virus/vrhyme_membership.tsv",
    params:
        pairs=vrhyme_pairs(),
        row_sample={r: sample_of(r) for r in contig_pairings},
        base="output/virus/vrhyme",
    conda:
        "../env/mag_qc.yaml"
    log:
        "output/logs/virus/vrhyme_summary.log"
    script:
        "../scripts/make_vrhyme_summary.py"
