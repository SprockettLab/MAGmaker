#!/usr/bin/env python3
"""Tests the metadata loader in common.smk, including single-end handling.

common.smk defines no rules, so it can be exec'd directly with a stub
config and a stub expand(). That keeps these checks runnable without a
snakemake install or a cluster.

Run from the repository root:

    python3 resources/test/test_metadata_layout.py

The behaviour under test is a safety property as much as a feature: a
blank R2_fp must stay an error, and only the exact literal NA may declare
single-end input, so that a truncated sample sheet cannot be mistaken for
single-end data and quietly assemble half of a paired library.
"""

import os
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
COMMON = os.path.join(REPO_ROOT, 'resources', 'snakefiles', 'common.smk')

HEADER = "Sample\tSequencing_Run\tR1_fp\tR2_fp\tGroup"

PAIRED_ROW = "A\tRun_1\treads/A_R1.fq.gz\treads/A_R2.fq.gz\tctl"
SINGLE_ROW = "B\tRun_1\treads/B_R1.fq.gz\tNA\tctl"


def load(rows, header=HEADER):
    """in : metadata rows
       out: (namespace, None) on success, (None, error message) on failure"""
    fd, path = tempfile.mkstemp(suffix='.txt')
    os.close(fd)
    with open(path, 'w') as handle:
        handle.write(header + "\n" + "\n".join(rows) + "\n")

    namespace = {
        'config': {'metadata': path, 'reads': ['R1', 'R2']},
        'expand': lambda template, **kw: [
            template.format(sample=kw['sample'], read=read)
            for read in kw['read']
        ],
    }
    try:
        exec(compile(open(COMMON).read(), COMMON, 'exec'), namespace)
        return namespace, None
    except SystemExit as exc:
        return None, str(exc)
    finally:
        os.unlink(path)


FAILURES = []


def check(name, condition):
    print(("PASS  " if condition else "FAIL  ") + name)
    if not condition:
        FAILURES.append(name)


def main():
    # Paired-end behaviour is unchanged.
    ns, err = load([PAIRED_ROW])
    check("paired sample loads", err is None and not ns['is_single_end']('A'))
    check("paired reads_for is R1,R2",
          err is None and ns['reads_for']('A') == ['R1', 'R2'])

    # NA declares single-end.
    ns, err = load([SINGLE_ROW])
    check("NA in R2_fp declares single-end",
          err is None and ns['is_single_end']('B'))
    check("single-end reads_for is SE",
          err is None and ns['reads_for']('B') == ['SE'])
    check("get_read('SE') falls back to R1_fp",
          err is None
          and ns['get_read']('B', 'Run_1', 'SE') == 'reads/B_R1.fq.gz')
    check("nonhost_reads returns one file for single-end",
          err is None and ns['nonhost_reads']('B')
          == ['output/qc/host_filter/nonhost/B.SE.fastq.gz'])

    # A cohort may mix single-end and paired-end samples.
    ns, err = load([PAIRED_ROW, SINGLE_ROW])
    check("mixed cohort loads", err is None)
    check("metaspades skips single-end samples",
          err is None and ns['samples_for_assembler']('metaspades') == ['A'])
    check("megahit keeps every sample",
          err is None
          and sorted(ns['samples_for_assembler']('megahit')) == ['A', 'B'])

    # A blank R2_fp remains an error: this is the whole safety property.
    for value, label in [("", "empty"), ("None", "None"), ("nan", "nan")]:
        ns, err = load(["C\tRun_1\treads/C_R1.fq.gz\t%s\tctl" % value])
        check("blank R2_fp (%s) is rejected" % label,
              err is not None and "R2_fp" in err)

    # NA is only meaningful in R2_fp.
    ns, err = load(["D\tRun_1\tNA\treads/D_R2.fq.gz\tctl"])
    check("NA in R1_fp is rejected", err is not None)

    ns, err = load(["\tRun_1\treads/E_R1.fq.gz\tNA\tctl"])
    check("blank Sample is rejected", err is not None)

    # The marker is exact and case-sensitive.
    ns, err = load(["F\tRun_1\treads/F_R1.fq.gz\tna\tctl"])
    check("lowercase na does not declare single-end",
          err is not None or not ns['is_single_end']('F'))

    # merge_seqruns concatenates per read id, so layout cannot vary within
    # a sample.
    ns, err = load(["G\tRun_1\treads/G1_R1.fq.gz\treads/G1_R2.fq.gz\tctl",
                    "G\tRun_2\treads/G2_R1.fq.gz\tNA\tctl"])
    check("mixed layout within one sample is rejected",
          err is not None and "mixes single-end and paired-end" in err)

    print("\n%d failure(s)" % len(FAILURES))
    return 1 if FAILURES else 0


if __name__ == '__main__':
    sys.exit(main())
