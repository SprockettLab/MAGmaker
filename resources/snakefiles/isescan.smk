# Insertion sequence (IS) elements with ISEScan.
#
# Runs on every contig of each binning row's own assembly, not on the MAGs:
# IS elements are multi-copy repeats, so assemblers often break at them and
# binners often leave the resulting contigs out of bins. Running on the bins
# would undercount them, worst in the IS-rich strains. Each call is then
# labelled with the MAG its contig was binned into (if any), and a summary
# by species-level genome bin can be built downstream from that table without
# running ISEScan again.
#
# Off by default. Turn on with `isescan: true` in the config, or for one run
# with run_magmaker.sh --isescan (docs/isescan.md).
#
# Runtime measured on WF18 (one assembly per job, 32 threads, ISEScan 1.7.3):
#   flye v2980          82.1 Mbp,   1,526 contigs   19.6 min   607 MB
#   MEGAHIT (median)   214.8 Mbp, 241,571 contigs   82.8 min   753 MB
#   myloasm v2193     2,261 Mbp,  60,457 contigs   6.46 h     5.4 GB
# 10-23 s per Mbp. The largest of 96 WF18 assemblies fits in one job, so the
# rule does not split assemblies.


def isescan_on():
    return bool(config.get('isescan', False))


def isescan_pairs():
    """(assembler, binning row): each row's own assembly, the one its MAGs
    came from, so every call can be labelled with its MAG."""
    if not isescan_on():
        return []
    return [(assembler_of(r), r) for r in contig_pairings]


rule run_isescan:
    """ISEScan on every contig of one assembly."""
    input:
        contigs = lambda wildcards: "output/assemble/{0}/{1}.contigs.fasta".format(
            wildcards.assembler, sample_of(wildcards.contig_sample))
    output:
        tsv = "output/isescan/{assembler}/{contig_sample}/is.tsv",
        done = touch("output/isescan/{assembler}/{contig_sample}/.done")
    params:
        outdir = "output/isescan/{assembler}/{contig_sample}",
        extra = config['params'].get('isescan', {}).get('extra', '')
    threads:
        config['threads'].get('isescan', 32)
    retries:
        config['retries'].get('run_isescan', 2)
    resources:
        mem_mb = mem_escalate('run_isescan', base_default=16000),
        runtime = runtime_escalate('run_isescan', base_default=1440)
    conda:
        "../env/isescan.yaml"
    benchmark:
        "output/benchmarks/isescan/{assembler}/{contig_sample}.txt"
    log:
        "output/logs/isescan/{assembler}/{contig_sample}.log"
    shell:
        # The output folder is cleared on every attempt: ISEScan skips any
        # HMM search whose result file already exists, so a retry after a
        # killed attempt would build its calls from a half-written search.
        #
        # ISEScan writes a file list into the working directory, named after
        # the input's folder and file, so it runs from inside its own output
        # folder with the input given as an absolute path.
        #
        # Its results land in <out>/run/<input folder>/<input file>.tsv; the
        # one .tsv outside proteome/ and hmm/ is copied to is.tsv. If its HMM
        # search found nothing it writes no table at all, and logs "No hit was
        # returned by HMM search"; that case gets an empty table plus a
        # NO_IS_FOUND note, so "found nothing" and "wrote nothing" stay apart.
        # Any other missing table fails the job.
        #
        # proteome/ and hmm/ (gene calls and raw HMM hits, the bulk of the
        # output) are removed after a successful run.
        """
            rm -rf {params.outdir}
            mkdir -p {params.outdir}/run
            contigs=$(realpath {input.contigs})
            log=$(realpath -m {log})
            tsv=$(realpath -m {output.tsv})
            cd {params.outdir}/run
            isescan.py --seqfile "$contigs" --output . --nthread {threads} \
                {params.extra} > "$log" 2>&1
            found=$(find . -name '*.tsv' -not -path './proteome/*' -not -path './hmm/*')
            n=$(printf '%s' "$found" | grep -c . || true)
            if [ "$n" -eq 1 ]; then
                cp "$found" "$tsv"
            elif [ "$n" -eq 0 ] && grep -q 'No hit was returned by HMM search' "$log"; then
                : > "$tsv"
                echo "ISEScan finished; its HMM search found no transposase hits" > ../NO_IS_FOUND
            else
                echo "MAGmaker: expected one ISEScan table, found $n" >> "$log"
                exit 1
            fi
            rm -rf proteome hmm
            {{
              echo "tool	ISEScan"
              echo "version	$(isescan.py --version 2>&1 | tail -n 1)"
              echo "extra	{params.extra}"
              echo "threads	{threads}"
              echo "input	{input.contigs}"
              echo "input_contigs	$(grep -c '^>' "$contigs" || true)"
              echo "is_rows	$(( $(wc -l < "$tsv") > 0 ? $(wc -l < "$tsv") - 1 : 0 ))"
            }} > ../params.txt
        """


rule make_isescan_summary:
    """
    Every row's ISEScan table in one file, each call labelled with its row,
    sample, assembler and, where its contig was binned, the MAG and taxonomy
    from mag_summary.tsv.
    """
    input:
        tables = lambda wildcards: ["output/isescan/{0}/{1}/is.tsv".format(a, r)
                                    for a, r in isescan_pairs()],
        mag_summary = "output/mag_qc/mag_summary.tsv"
    output:
        summary = "output/isescan/isescan_summary.tsv"
    params:
        pairs = lambda wildcards: isescan_pairs(),
        row_sample = {r: sample_of(r) for r in contig_pairings}
    conda:
        "../env/mag_qc.yaml"
    log:
        "output/logs/isescan/make_isescan_summary.log"
    script:
        "../scripts/make_isescan_summary.py"


rule isescan_all:
    input:
        (["output/isescan/isescan_summary.tsv"] if isescan_on() else [])
