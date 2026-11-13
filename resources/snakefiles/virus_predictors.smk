# Viral predictors other than geNomad, and CheckV on each one's own calls.
#
# Included by virus.smk, which holds the switches (config['virus']), the
# per-predictor output paths (predictor_fasta) and tool_env().
#
# Each tool runs in an environment that already exists, named in
# params.<tool>.env, rather than from a conda yaml. Their dependencies cannot
# be re-solved safely: VIBRANT needs scikit-learn 0.21.3 on a numpy that still
# has np.float, and VirSorter2 needs python 3.8 with scikit-learn 0.22.1, the
# version its classifier was trained under. See docs/viral.md.
#
# Every rule clears its output folder before running, as run_genomad does: all
# three resume from whatever they find, so a retry after a killed attempt would
# otherwise mix new results with half-written ones. A tool that finishes
# without writing its viral FASTA gets an empty one plus a NO_VIRAL_FASTA note,
# so "called nothing" and "wrote nothing" stay apart in the summary.

wildcard_constraints:
    predictor="vibrant|cenotetaker3|virsorter2",


def _predictor_contigs(wildcards):
    return "output/assemble/{0}/{1}.contigs.fasta".format(
        wildcards.assembler, sample_of(wildcards.contig_sample))


def cenotetaker3_mem_mb(wildcards, input, attempt):
    """Cenote-Taker3 memory from the size of the contigs FASTA: the larger of
    config mem_mb.run_cenotetaker3 and _per_mb MB per MiB, times the attempt
    number, capped at _max. A flat 32 GB fitted every MEGAHIT assembly (<= 709
    Mbp) but OOM-killed all six 1.2-2.2 GiB long-read assemblies in
    WF18_FullTest_2026; see config.yaml for the measurements."""
    mem = config["mem_mb"]
    contigs_mib = os.path.getsize(input.contigs) / 1048576
    base = max(mem.get("run_cenotetaker3", 32000),
               int(contigs_mib * mem.get("run_cenotetaker3_per_mb", 60)))
    return min(mem.get("run_cenotetaker3_max", 700000), base * attempt)


rule run_vibrant:
    """VIBRANT phage prediction for one assembly."""
    input:
        contigs=_predictor_contigs,
    output:
        done=touch("output/virus/vibrant/{assembler}/{contig_sample}/.done"),
    params:
        outdir="output/virus/vibrant/{assembler}/{contig_sample}",
        fasta=lambda w: predictor_fasta('vibrant', w.assembler, w.contig_sample),
        db=config['params'].get('vibrant', {}).get('db_dir', ''),
        files=config['params'].get('vibrant', {}).get('files_dir', ''),
    threads: config['threads'].get('vibrant', 16)
    retries: config['retries'].get('run_vibrant', 2)
    resources:
        mem_mb=mem_escalate('run_vibrant', base_default=32000),
        runtime=runtime_escalate('run_vibrant', base_default=480),
    conda:
        tool_env('vibrant') if VIRUS_ON['vibrant'] else "../env/mag_qc.yaml"
    benchmark:
        "output/benchmarks/virus/vibrant/{assembler}/{contig_sample}.txt"
    log:
        "output/logs/virus/vibrant/{assembler}/{contig_sample}.log"
    shell:
        # -d is the databases folder (.HMM files) and -m the files folder
        # (.tsv and model files). Both are passed because the databases live
        # outside VIBRANT's package directory.
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        VIBRANT_run.py -i {input.contigs} -t {threads} \
            -d {params.db} -m {params.files} -folder {params.outdir} \
            > {log} 2>&1
        if [ ! -e {params.fasta} ]; then
            mkdir -p "$(dirname {params.fasta})"
            : > {params.fasta}
            echo "VIBRANT wrote no viral FASTA; recorded as no output" \
                > {params.outdir}/NO_VIRAL_FASTA
        fi
        {{
          echo "tool\tVIBRANT"
          echo "version\t$(VIBRANT_run.py --version 2>&1 | tail -n 1)"
          echo "env\t${{CONDA_PREFIX:-unknown}}"
          echo "databases\t{params.db}"
          echo "files\t{params.files}"
          echo "threads\t{threads}"
          echo "input\t{input.contigs}"
          echo "input_contigs\t$(grep -c '^>' {input.contigs} || true)"
          echo "sequences_called\t$(grep -c '^>' {params.fasta} || true)"
        }} > {params.outdir}/params.txt
        """


rule run_cenotetaker3:
    """Cenote-Taker3 virus discovery for one assembly."""
    input:
        contigs=_predictor_contigs,
    output:
        done=touch("output/virus/cenotetaker3/{assembler}/{contig_sample}/.done"),
    params:
        outdir="output/virus/cenotetaker3/{assembler}/{contig_sample}",
        fasta=lambda w: predictor_fasta('cenotetaker3', w.assembler, w.contig_sample),
        db=config['params'].get('cenotetaker3', {}).get('db_dir', ''),
        prune=config['params'].get('cenotetaker3', {}).get('prune', 'T'),
    threads: config['threads'].get('cenotetaker3', 16)
    retries: config['retries'].get('run_cenotetaker3', 2)
    resources:
        mem_mb=cenotetaker3_mem_mb,
        runtime=runtime_escalate('run_cenotetaker3', base_default=480),
    conda:
        tool_env('cenotetaker3') if VIRUS_ON['cenotetaker3'] else "../env/mag_qc.yaml"
    benchmark:
        "output/benchmarks/virus/cenotetaker3/{assembler}/{contig_sample}.txt"
    log:
        "output/logs/virus/cenotetaker3/{assembler}/{contig_sample}.log"
    shell:
        # Flags verified against `cenotetaker3 --help` (3.4.4): -wd takes ONE
        # dash, and -r must be under 18 characters of letters, digits and
        # underscores, so it is the fixed "ct3" and not the sample name. The
        # input is made absolute because the tool works inside -wd.
        # -p T prunes host flanks from non-circular contigs, which upstream
        # advises for material that is not virus-enriched.
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        export CENOTE_DBS={params.db}
        cenotetaker3 -c "$(realpath {input.contigs})" -r ct3 -p {params.prune} \
            -t {threads} -wd {params.outdir} --cenote-dbs {params.db} \
            > {log} 2>&1
        if [ ! -e {params.fasta} ]; then
            mkdir -p "$(dirname {params.fasta})"
            : > {params.fasta}
            echo "Cenote-Taker3 wrote no viral FASTA; recorded as no output" \
                > {params.outdir}/NO_VIRAL_FASTA
        fi
        {{
          echo "tool\tCenote-Taker3"
          echo "version\t$(cenotetaker3 --version 2>&1 | tail -n 1)"
          echo "env\t${{CONDA_PREFIX:-unknown}}"
          echo "databases\t{params.db}"
          echo "prune\t{params.prune}"
          echo "run_title\tct3"
          echo "threads\t{threads}"
          echo "input\t{input.contigs}"
          echo "input_contigs\t$(grep -c '^>' {input.contigs} || true)"
          echo "sequences_called\t$(grep -c '^>' {params.fasta} || true)"
        }} > {params.outdir}/params.txt
        """


rule run_virsorter2:
    """VirSorter2 virus prediction for one assembly."""
    input:
        contigs=_predictor_contigs,
    output:
        done=touch("output/virus/virsorter2/{assembler}/{contig_sample}/.done"),
    params:
        outdir="output/virus/virsorter2/{assembler}/{contig_sample}",
        fasta=lambda w: predictor_fasta('virsorter2', w.assembler, w.contig_sample),
        db=config['params'].get('virsorter2', {}).get('db_dir', ''),
    threads: config['threads'].get('virsorter2', 16)
    retries: config['retries'].get('run_virsorter2', 2)
    resources:
        mem_mb=mem_escalate('run_virsorter2', base_default=32000),
        # Roughly ten times the other predictors: 27 min on a 10.3 Mb
        # assembly and unfinished at 81 min on three of 297 to 620 Mb
        # (WF30, 2026-11-11, 16 jobs sharing nodes, so contended).
        runtime=runtime_escalate('run_virsorter2', base_default=2880),
    conda:
        tool_env('virsorter2') if VIRUS_ON['virsorter2'] else "../env/mag_qc.yaml"
    benchmark:
        "output/benchmarks/virus/virsorter2/{assembler}/{contig_sample}.txt"
    log:
        "output/logs/virus/virsorter2/{assembler}/{contig_sample}.log"
    shell:
        # --use-conda-off is right ONLY because params.virsorter2.env is the
        # single environment built to VirSorter2's own per-rule pins. Its
        # per-rule envs never activate in a non-interactive shell, so rules
        # otherwise fall through to whatever python is on PATH and die on
        # `import screed`; and --use-conda-off in any other env would classify
        # under a scikit-learn the model was not trained with, turning a crash
        # into a silent wrong answer. The check below fails the job instead.
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        python -c "import screed, sklearn, sys; v = sklearn.__version__; \
sys.exit('' if v.startswith('0.22') else 'scikit-learn %s, expected 0.22.x: \
params.virsorter2.env is not the single env VirSorter2 needs' % v)" \
            > {log} 2>&1
        virsorter run --db-dir {params.db} -i {input.contigs} \
            -w {params.outdir} -j {threads} all --use-conda-off \
            >> {log} 2>&1
        if [ ! -e {params.fasta} ]; then
            : > {params.fasta}
            echo "VirSorter2 wrote no viral FASTA; recorded as no output" \
                > {params.outdir}/NO_VIRAL_FASTA
        fi
        {{
          echo "tool\tVirSorter2"
          echo "version\t$(virsorter --version 2>&1 | tail -n 1)"
          echo "env\t${{CONDA_PREFIX:-unknown}}"
          echo "database\t{params.db}"
          echo "mode\tall --use-conda-off"
          echo "threads\t{threads}"
          echo "input\t{input.contigs}"
          echo "input_contigs\t$(grep -c '^>' {input.contigs} || true)"
          echo "sequences_called\t$(grep -c '^>' {params.fasta} || true)"
        }} > {params.outdir}/params.txt
        """


rule run_checkv_predictor:
    """CheckV on one predictor's own sequences.

    CheckV runs after each predictor, on that predictor's calls, because it
    corrects predicted boundaries and the tools differ in exactly that.
    geNomad's CheckV is run_checkv in virus.smk and keeps its old path, so
    runs that already have it redo nothing.
    """
    input:
        "output/virus/{predictor}/{assembler}/{contig_sample}/.done",
    output:
        done=touch("output/virus/checkv_{predictor}/{assembler}/{contig_sample}/.done"),
    params:
        fasta=lambda w: predictor_fasta(w.predictor, w.assembler, w.contig_sample),
        outdir="output/virus/checkv_{predictor}/{assembler}/{contig_sample}",
        db=config['params']['checkv']['db_path'],
        # Relative to the working directory, as profile.smk does for its own
        # helper: a run folder is a clone, so resources/ sits beside output/.
        script="resources/scripts/unique_fasta_ids.py",
    threads: config['threads']['checkv']
    retries: config['retries'].get('run_checkv', 2)
    resources:
        mem_mb=config['mem_mb'].get('run_checkv', 16000),
        runtime=runtime_escalate('run_checkv', base_default=240),
    conda:
        "../env/checkv.yaml"
    benchmark:
        "output/benchmarks/virus/checkv_{predictor}/{assembler}/{contig_sample}.txt"
    log:
        "output/logs/virus/checkv_{predictor}/{assembler}/{contig_sample}.log"
    shell:
        # unique_fasta_ids.py renames repeated IDs (VIBRANT gives one contig's
        # several phage regions the same ID, which CheckV refuses outright).
        # The row count guards CheckV's silent short tables: on WF18,
        # 2026-11-10, 9 of 90 jobs exited 0 having graded only part of their
        # input, leaving 44,170 of 380,306 contigs with an empty quality.
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        python {params.script} {params.fasta} {params.outdir}/input.fna > {log} 2>&1
        n_in=$(grep -c '^>' {params.outdir}/input.fna || true)
        if [ "$n_in" -gt 0 ]; then
            checkv end_to_end {params.outdir}/input.fna {params.outdir} \
                -t {threads} -d {params.db} >> {log} 2>&1
            n_rows=$(( $(wc -l < {params.outdir}/quality_summary.tsv) - 1 ))
            if [ "$n_rows" -lt "$n_in" ]; then
                echo "ERROR: CheckV graded $n_rows of $n_in sequences; this is" \
                     "the silent short-table failure" >> {log}
                exit 1
            fi
        else
            echo "This predictor called nothing; CheckV skipped." >> {log}
        fi
        """
