#!/usr/bin/env python3
"""Checks the hybrid assembly arm (rule metaspades_hybrid).

1. Source assertions on resources/snakefiles/assemble.smk, so a later
   rewrite cannot quietly change what the arm is:
   - it passes --nanopore;
   - its long-read input is the HOST-FILTERED ONT FASTQ, not raw reads;
   - it writes to its own assembler directory, not metaspades';
   - it keeps its own temp directory, so a hybrid and a short-read-only
     assembly of the same sample cannot collide;
   - rule metaspades itself is untouched by --nanopore.

2. Config: the hybrid_assemblers switch exists and is OFF by default, and
   the resource keys assembler_mem() dereferences are present. Without
   mem_mb.metaspades_hybrid the rule would KeyError at DAG construction.

3. Target expansion in Snakefile.smk is over long_read_samples, NOT over
   every sample. This is the one that matters most: expanding over all
   samples would assemble short-read-only samples with no long reads,
   producing an assembly identical to the metaspades arm under a different
   name, and any comparison between the arms would then be meaningless
   while looking fine.

Run from the repository root:

    python3 resources/test/test_hybrid_assembly.py
"""

import os
import re
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
ASSEMBLE = os.path.join(REPO_ROOT, 'resources', 'snakefiles', 'assemble.smk')
SNAKEFILE = os.path.join(REPO_ROOT, 'resources', 'snakefiles', 'Snakefile.smk')
CONFIG = os.path.join(REPO_ROOT, 'resources', 'config', 'config.yaml')

FAILURES = []


def check(name, condition, detail=''):
    print(("PASS  " if condition else "FAIL  ") + name)
    if not condition:
        FAILURES.append(name)
        if detail:
            print("      " + str(detail))


def rule_body(text, name):
    """The source of one rule, from its header to the next top-level rule."""
    m = re.search(rf'^rule {name}:$', text, re.M)
    if not m:
        return ''
    rest = text[m.end():]
    nxt = re.search(r'^rule \w+:$', rest, re.M)
    return rest[:nxt.start()] if nxt else rest


def main():
    assemble = open(ASSEMBLE).read()
    snakefile = open(SNAKEFILE).read()
    config_text = open(CONFIG).read()

    print('rule metaspades_hybrid:')
    body = rule_body(assemble, 'metaspades_hybrid')
    check('the rule exists', bool(body))
    if body:
        check('passes --nanopore', '--nanopore' in body)
        check('long reads come from the host-filtered ONT FASTQ',
              'output/qc/host_filter_ont/nonhost/{sample}.ont.fastq.gz' in body,
              'raw long reads would reintroduce host sequence the ONT host '
              'filter exists to remove')
        check('writes to its own assembler directory',
              'output/assemble/metaspades_hybrid/{sample}.contigs.fasta' in body)
        check('uses both short-read mates',
              '--pe1-1' in body and '--pe1-2' in body)
        check('its temp dir cannot collide with rule metaspades',
              '_hybrid_temp' in body,
              'both rules rm -rf their temp dir, so a shared path would let '
              'one delete the other mid-assembly')
        check('clears its temp dir before assembling', 'rm -rf' in body)
        check('memory comes from its own key',
              "assembler_mem('metaspades_hybrid')" in body)

    print()
    print('rule metaspades (must be unchanged by this feature):')
    plain = rule_body(assemble, 'metaspades')
    check('rule metaspades still exists', bool(plain))
    check('rule metaspades does NOT pass --nanopore',
          '--nanopore' not in plain,
          'the two arms must differ only in the assembly, so the short-read '
          'arm cannot quietly become hybrid')
    check('rule metaspades keeps its own output path',
          'output/assemble/metaspades/{sample}.contigs.fasta' in plain)

    print()
    print('config:')
    check('hybrid_assemblers key exists', 'hybrid_assemblers:' in config_text)
    check('hybrid_assemblers is OFF by default',
          re.search(r'^hybrid_assemblers:\s*\[\]\s*$', config_text, re.M)
          is not None,
          'an arm that costs 1.5x metaSPAdes memory must be opt-in')
    for key in ('metaspades_hybrid', 'metaspades_hybrid_min',
                'metaspades_hybrid_per_mb'):
        check(f'mem_mb.{key} is set',
              re.search(rf'^\s+{key}:\s*\d+\s*$', config_text, re.M) is not None,
              'assembler_mem() reads mem_mb[key] directly and KeyErrors '
              'without the bare key')

    print()
    print('target expansion:')
    lrt = snakefile[snakefile.find('long_read_targets'):]
    lrt = lrt[:lrt.find('rule all')]
    check('hybrid_assemblers is expanded', 'hybrid_assemblers' in lrt)
    hybrid_lines = [ln for ln in lrt.splitlines() if 'hybrid_assemblers' in ln]
    check('expanded over long_read_samples, not all samples',
          all('sample=long_read_samples' in ln
              or 'sample=long_read_samples' in lrt.splitlines()[
                  lrt.splitlines().index(ln) + 1]
              for ln in hybrid_lines),
          'expanding over every sample would assemble short-read-only '
          'samples with no long reads and produce a duplicate of the '
          'metaspades arm under a different name')
    check('samples= is never used for the hybrid expansion',
          not re.search(r'hybrid_assemblers.*\n?.*sample=samples', lrt))

    print()
    if FAILURES:
        print(f'{len(FAILURES)} FAILED: ' + ', '.join(FAILURES))
        return 1
    print('test_hybrid_assembly: all checks passed')
    return 0


if __name__ == '__main__':
    sys.exit(main())
