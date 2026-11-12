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


def virus_pairs():
    """(assembler, binning row) pairs the viral track runs on: every configured
    short-read assembler for a short-read row, as before; a long-read row's
    own assembler only."""
    pairs = []
    for row in contig_pairings:
        if is_long_read_row(row):
            pairs.append((assembler_of(row), row))
        else:
            pairs += [(a, row) for a in config['assemblers']]
    return pairs


def contigs_stem(row):
    """geNomad names its outputs after the input FASTA, <sample>.contigs; for a
    long-read row <sample>__<assembler> that is the source sample's name."""
    return sample_of(row) + ".contigs"


# ── Which tools run ──────────────────────────────────────────────────────────
#
# config['virus'] has one switch per tool. Per run, --config
# virus_predictors=[genomad,vibrant] replaces the predictor switches and
# skip_vrhyme=True turns vRhyme off (run_magmaker.sh: --virus-predictors=...,
# --skip-vrhyme).

VIRUS_PREDICTORS = ['genomad', 'vibrant', 'cenotetaker3', 'virsorter2']
OTHER_PREDICTORS = ['vibrant', 'cenotetaker3', 'virsorter2']
VIRUS_DEFAULTS = {'genomad': True, 'vrhyme': True, 'vibrant': False,
                  'cenotetaker3': False, 'virsorter2': False}


def _virus_switches():
    switches = dict(VIRUS_DEFAULTS)
    switches.update(config.get('virus') or {})
    unknown = sorted(set(switches) - set(VIRUS_DEFAULTS))
    if unknown:
        raise ValueError(f"config virus: unknown switch(es) {unknown}; "
                         f"known: {sorted(VIRUS_DEFAULTS)}")
    chosen = config.get('virus_predictors')
    if chosen is not None:
        if isinstance(chosen, str):
            chosen = [c for c in chosen.replace(',', ' ').split() if c]
        bad = sorted(set(chosen) - set(VIRUS_PREDICTORS))
        if bad:
            raise ValueError(f"virus_predictors: unknown predictor(s) {bad}; "
                             f"known: {VIRUS_PREDICTORS}")
        for tool in VIRUS_PREDICTORS:
            switches[tool] = tool in chosen
    if config.get('skip_vrhyme'):
        switches['vrhyme'] = False
    return switches


VIRUS_ON = _virus_switches()
PREDICTORS_ON = [t for t in VIRUS_PREDICTORS if VIRUS_ON[t]]


def vrhyme_mapper():
    return (config['params'].get('vrhyme', {}).get('mapper') or config['mappers'][0])


def _vrhyme_runs():
    """vRhyme bins geNomad's calls with MetaBAT2's coverage table, so it needs
    both. Said once at parse time when it cannot run, rather than failing a
    run whose other parts are fine."""
    if not VIRUS_ON['vrhyme']:
        return False
    missing = []
    if not VIRUS_ON['genomad']:
        missing.append("geNomad is off (vRhyme bins geNomad's viral contigs)")
    if 'metabat2' not in config.get('binners', []):
        missing.append("metabat2 is not among the binners (vRhyme uses its "
                       "coverage table; without it the reads would be mapped again)")
    if missing:
        print("MAGmaker viral track: vRhyme is on but will NOT run: "
              + "; ".join(missing), file=sys.stderr)
        return False
    return True


VRHYME_RUNS = _vrhyme_runs()
if VIRUS_ON:
    print("MAGmaker viral track: predictors " + (", ".join(PREDICTORS_ON) or "none")
          + "; vRhyme " + ("on" if VRHYME_RUNS else "off"), file=sys.stderr)


def vrhyme_pairs():
    """vRhyme needs the coverage table of the contigs it bins, which exists only
    for a row's own assembler (the binning table's Contigs column)."""
    if not VRHYME_RUNS:
        return []
    return [(a, r) for a, r in virus_pairs() if a == assembler_of(r)]


def predictor_fasta(predictor, assembler, row):
    """Each predictor's viral sequences, the file CheckV reads."""
    stem = contigs_stem(row)
    base = f"output/virus/{predictor}/{assembler}/{row}"
    if predictor == 'genomad':
        return f"{base}/{stem}_summary/{stem}_virus.fna"
    if predictor == 'vibrant':
        return f"{base}/VIBRANT_{stem}/VIBRANT_phages_{stem}/{stem}.phages_combined.fna"
    if predictor == 'cenotetaker3':
        # -r ct3: Cenote-Taker3 names its run folder and files after it.
        return f"{base}/ct3/ct3_virus_sequences.fna"
    if predictor == 'virsorter2':
        return f"{base}/final-viral-combined.fa"
    raise ValueError(predictor)


def checkv_dir(predictor, assembler, row):
    """geNomad's CheckV keeps the path it always had, so runs that already
    have it rerun nothing; the other predictors' CheckV sits beside it."""
    if predictor == 'genomad':
        return f"output/virus/checkv/{assembler}/{row}"
    return f"output/virus/checkv_{predictor}/{assembler}/{row}"


def tool_env(tool):
    """An existing environment folder named in params.<tool>.env."""
    env = config['params'].get(tool, {}).get('env')
    if not env:
        raise ValueError(f"params.{tool}.env is not set; it must name the "
                         f"environment folder {tool} runs in (docs/viral.md)")
    return env


rule run_genomad:
    """
    geNomad end-to-end: identify viral + plasmid contigs and classify them.
    """
    input:
        contigs = lambda wildcards: "output/assemble/{0}/{1}.contigs.fasta".format(
            wildcards.assembler, sample_of(wildcards.contig_sample))
    output:
        done = touch("output/virus/genomad/{assembler}/{contig_sample}/.done")
    params:
        outdir = "output/virus/genomad/{assembler}/{contig_sample}",
        db = config['params']['genomad']['db_path'],
        extra = config['params']['genomad']['extra'],
        # geNomad reports FDR only when scores are calibrated, and that also
        # turns on its --max-fdr filter. See params.genomad.score_calibration.
        calibration = ('--enable-score-calibration'
                       if config['params']['genomad'].get('score_calibration', True)
                       else '')
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
            genomad end-to-end --threads {threads} {params.calibration} {params.extra} \
                {input.contigs} {params.outdir} {params.db} \
                2> {log} 1>&2
            rm -rf {params.outdir}/*_annotate \
                   {params.outdir}/*_aggregated_classification \
                   {params.outdir}/*_find_proviruses \
                   {params.outdir}/*_marker_classification \
                   {params.outdir}/*_nn_classification \
                   2> /dev/null || true
            {{
              echo "tool	geNomad"
              echo "version	$(genomad --version 2>&1 | tail -n 1)"
              echo "database	{params.db}"
              echo "score_calibration_flag	{params.calibration}"
              echo "extra	{params.extra}"
              echo "threads	{threads}"
              echo "input	{input.contigs}"
              echo "input_contigs	$(grep -c '^>' {input.contigs} || true)"
            }} > {params.outdir}/params.txt
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
        virus_fna = lambda wildcards: (
            "output/virus/genomad/{0}/{1}/{2}_summary/{2}_virus.fna".format(
                wildcards.assembler, wildcards.contig_sample, contigs_stem(wildcards.contig_sample))),
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


include: "virus_predictors.smk"
include: "vrhyme.smk"


def predictor_targets():
    """Every switched-on predictor's .done, plus its own CheckV's."""
    targets = []
    for predictor in PREDICTORS_ON:
        for assembler, row in virus_pairs():
            targets.append("output/virus/{0}/{1}/{2}/.done".format(
                predictor, assembler, row))
            targets.append(checkv_dir(predictor, assembler, row) + "/.done")
    return targets


rule make_virus_summary:
    """
    Merge every predictor's calls, with its own CheckV's grades, into one table.

    One row per call per predictor, carrying Predictor and Source_Contig, so
    nothing downstream has to parse a tool's own naming. Backward compatible:
    with only geNomad on, the rows are the same as before plus those columns.
    """
    input:
        genomad = lambda wildcards: (
            ["output/virus/genomad/{0}/{1}/.done".format(a, r)
             for a, r in virus_pairs()] if VIRUS_ON['genomad'] else []),
        checkv = lambda wildcards: (
            ["output/virus/checkv/{0}/{1}/.done".format(a, r)
             for a, r in virus_pairs()] if VIRUS_ON['genomad'] else []),
        predictors = lambda wildcards: predictor_targets()
    output:
        summary = "output/virus/virus_summary.tsv",
        mapping = "output/virus/virus_name_mapping.tsv"
    params:
        contig_samples = list(contig_pairings.keys()),
        assemblers = config['assemblers'],
        pairs = virus_pairs(),
        predictors = [p for p in PREDICTORS_ON],
        row_sample = {r: sample_of(r) for r in contig_pairings},
        row_stem = {r: contigs_stem(r) for r in contig_pairings},
        long_read_rows = [r for r in contig_pairings if is_long_read_row(r)],
        genomad_base = "output/virus/genomad",
        checkv_base = "output/virus/checkv",
        virus_base = "output/virus"
    conda:
        "../env/mag_qc.yaml"
    log:
        "output/logs/virus/make_virus_summary.log"
    script:
        "../scripts/make_virus_summary.py"


rule virus_all:
    input:
        "output/virus/virus_summary.tsv",
        (["output/virus/vrhyme_summary.tsv"] if VRHYME_RUNS else [])
