#!/usr/bin/env python3
"""Plans the binning stage with a long-read row and checks every place a row
is not a sample.

A binning table with a short-read row (John) and a long-read row
(John__flye: John's Flye assembly) is dry-run through Snakefile-bin with
commands printed. Checks that the long-read row gets its own contigs, its
long-read mapping, SemiBin2's long-read mode and the circular candidates;
that CMSeq maps John's reads to it; that geNomad reads John's Flye contigs;
and that the short-read row's commands are untouched by any of this.

Needs snakemake on PATH. Run from the repository root:

    python3 resources/test/test_long_read_binning_plan.py
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
FAILURES = []


def check(name, condition, detail=''):
    print(("PASS  " if condition else "FAIL  ") + name)
    if not condition:
        FAILURES.append(name)
        if detail:
            print("      " + str(detail)[-1500:])


def jobs(out):
    """rule name -> list of (wildcards line, shell text) from a printed dry run."""
    found = {}
    for block in re.split(r"\n(?=\[[^\]]+\]\n(?:local)?rule )", out):
        # The rule name only from the job's own header line: printed shell
        # commands can mention other rules in their comments.
        m = re.match(r"\s*\[[^\]]+\]\n(?:local)?rule (\w+):", block)
        if m:
            found.setdefault(m.group(1), []).append(block)
    return found


def main():
    if not shutil.which('snakemake'):
        sys.exit("Needs snakemake on PATH")
    tmp = tempfile.mkdtemp(prefix='magmaker_lr_bin_', dir=os.environ.get('TMPDIR'))
    os.symlink(os.path.join(REPO_ROOT, 'resources'), os.path.join(tmp, 'resources'))
    open(os.path.join(tmp, 'host.fna'), 'w').write('>c1\nACGTACGTACGGTTAACCGG\n')
    open(os.path.join(tmp, 'lambda.fa'), 'w').write('>NC_001416.1\nACGTTGCA\n')
    meta = os.path.join(tmp, 'metadata.txt')
    with open(meta, 'w') as fh:
        fh.write('Sample\tSequencing_Run\tPlatform\tR1_fp\tR2_fp\n')
        for s in ('John', 'Paul'):
            fh.write(f'{s}\tRun_1\tillumina\tresources/test/test_reads/{s}_R1.fastq.gz\t'
                     f'resources/test/test_reads/{s}_R2.fastq.gz\n')
        fh.write(f'John\tONT1\tont\t{tmp}/John.ont.fastq.gz\tNA\n')
    open(os.path.join(tmp, 'John.ont.fastq.gz'), 'w').close()
    binning = os.path.join(tmp, 'binning.txt')
    with open(binning, 'w') as fh:
        fh.write('Sample\tContigs\tRead_Groups\tContig_Groups\tSource_Sample\tAssembler\n')
        fh.write('John\toutput/assemble/megahit/John.contigs.fasta\tGroup1\tGroup1\tJohn\tmegahit\n')
        fh.write('Paul\toutput/assemble/megahit/Paul.contigs.fasta\tGroup1\tGroup1\tPaul\tmegahit\n')
        fh.write('John__flye\toutput/assemble/flye/John.contigs.fasta\t\tGroup1\tJohn\tflye\n')
    with open(os.path.join(tmp, 'over.yaml'), 'w') as fh:
        yaml.safe_dump({'metadata': meta, 'binning': binning, 'long_read_assemblers': ['flye'],
                        'params': {'long_reads': {'control_strand_fasta': os.path.join(tmp, 'lambda.fa')}},
                        'host_filter': {'boardercollie': False,
                                        'genome': os.path.join(tmp, 'host.fna'), 'db_dir': tmp}}, fh)
    p = subprocess.run(['snakemake', '-s', os.path.join(REPO_ROOT, 'Snakefile-bin'),
                        '--configfile', os.path.join(REPO_ROOT, 'config.yaml'),
                        os.path.join(tmp, 'over.yaml'), '--cores', '1', '-d', tmp, '-n', '-p',
                        'rename_mags', 'virus_all'],
                       capture_output=True, text=True)
    out = p.stdout + p.stderr
    check("binning stage plans with a long-read row", p.returncode == 0, out)
    if p.returncode != 0:
        return
    j = jobs(out)
    # A job's row comes from its own wildcards line. Matching row names anywhere
    # in the printed text is not enough: notes Snakemake prints after a job
    # (e.g. temporary files it would remove) can name other rows, and jobs are
    # printed in a different order on each run.
    def row(block):
        m = re.search(r"^\s*wildcards: .*?contig_sample=([^,\s]+)", block, re.M)
        return m.group(1) if m else None
    lr = lambda rule: [b for b in j.get(rule, []) if row(b) == 'John__flye']
    sr = lambda rule: [b for b in j.get(rule, []) if row(b) == 'John']

    check("long-read mapping planned once, for the long-read row only",
          [row(b) for b in j.get('map_long_reads_to_contigs', [])] == ['John__flye'])
    check("circular candidates planned for the long-read row only",
          [row(b) for b in j.get('circular_candidates', [])] == ['John__flye'])
    mb = lr('make_metabat2_coverage_table')
    check("MetaBAT2 depth for the long-read row includes its long-read BAM",
          mb and 'output/mapping/long_reads/John__flye.bam' in mb[0], mb)
    check("MetaBAT2 depth for the short-read row has no long-read BAM",
          sr('make_metabat2_coverage_table') and 'mapping/long_reads' not in sr('make_metabat2_coverage_table')[0])
    sb = lr('run_semibin2')
    check("SemiBin2 runs in long-read mode for the long-read row",
          len(sb) == 1 and '--sequencing-type long_read' in sb[0]
          and 'output/assemble/flye/John.contigs.fasta' in sb[0], sb)
    check("SemiBin2 keeps its default mode for the short-read row",
          sr('run_semibin2') and '--sequencing-type' not in sr('run_semibin2')[0])
    cm = lr('cmseq_map')
    check("CMSeq maps John's own reads to the long-read row",
          cm and 'host_filter/nonhost/John.R1.fastq.gz' in cm[0], cm)
    gn = lr('run_genomad')
    check("geNomad reads John's Flye contigs, once",
          len(gn) == 1 and 'output/assemble/flye/John.contigs.fasta' in gn[0], gn)
    check("geNomad does not run megahit for the long-read row",
          not any('assembler=megahit' in b for b in gn))
    tables = lr('stage_binette_inputs')
    check("the circular candidate table reaches Binette",
          tables and 'selected_bins/circular/' in tables[0], tables)
    shutil.rmtree(tmp)


if __name__ == '__main__':
    main()
    print("\n%d failure(s)" % len(FAILURES))
    sys.exit(1 if FAILURES else 0)
