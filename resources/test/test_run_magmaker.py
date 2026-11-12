#!/usr/bin/env python3
"""Checks that run_magmaker.sh runs the viral / plasmid track by default.

Stage 3 asks for rename_mags and virus_all (geNomad, CheckV,
virus_summary.tsv); --skip-virus drops virus_all for that run only.

1. The wrapper's Snakemake command lines, recorded by a stand-in `snakemake`
   on PATH: stage 3 names both targets by default, only rename_mags with
   --skip-virus, and the other --skip-* overrides still reach stage 3.
2. A real Snakemake dry run of Snakefile-bin with both targets on the test
   data plans geNomad and CheckV for every assembly, and a dry run of
   rename_mags alone plans neither (so the first check can fail).

Run from the repository root (needs snakemake on PATH for part 2):

    python3 resources/test/test_run_magmaker.py
"""

import os
import re
import shutil
import stat
import subprocess
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
WRAPPER = os.path.join(REPO_ROOT, 'run_magmaker.sh')
TEST_SAMPLES = ['John', 'Paul', 'George', 'Ringo']

FAKE_SNAKEMAKE = """#!/bin/bash
printf '%s\\n' "$*" >> "$CALLS"
"""


def wrapper_calls(*flags):
    """Run the wrapper in a scratch folder with a stand-in snakemake; return
    one string per Snakemake call."""
    with tempfile.TemporaryDirectory(prefix='magmaker_wrapper_') as tmp:
        for name in ('Snakefile', 'Snakefile-bin'):
            open(os.path.join(tmp, name), 'w').close()
        os.makedirs(os.path.join(tmp, 'output', 'config'))
        open(os.path.join(tmp, 'output', 'config', 'auto_binning.txt'), 'w').close()
        bindir = os.path.join(tmp, 'bin')
        os.makedirs(bindir)
        fake = os.path.join(bindir, 'snakemake')
        with open(fake, 'w') as handle:
            handle.write(FAKE_SNAKEMAKE)
        os.chmod(fake, os.stat(fake).st_mode | stat.S_IEXEC)
        shutil.copy(WRAPPER, os.path.join(tmp, 'run_magmaker.sh'))
        calls = os.path.join(tmp, 'calls.txt')
        env = {**os.environ, 'PATH': bindir + os.pathsep + os.environ['PATH'],
               'CALLS': calls}
        result = subprocess.run(
            ['bash', 'run_magmaker.sh', '--profile', 'prof', *flags],
            cwd=tmp, env=env, capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        with open(calls) as handle:
            return handle.read().splitlines()


def stage3_targets(call):
    """Positional words between Snakefile-bin and the first option."""
    words = call.split()
    assert words[:2] == ['--snakefile', 'Snakefile-bin'], call
    targets = []
    for word in words[2:]:
        if word.startswith('-'):
            break
        targets.append(word)
    return targets


def test_wrapper_targets():
    calls = wrapper_calls()
    assert len(calls) == 3, calls
    assert stage3_targets(calls[2]) == ['rename_mags', 'virus_all'], calls[2]
    assert 'virus_all' not in calls[0] and 'virus_all' not in calls[1]

    calls = wrapper_calls('--skip-virus', '--skip-cmseq')
    assert stage3_targets(calls[2]) == ['rename_mags'], calls[2]
    # --skip-virus is consumed by the wrapper, never passed to Snakemake,
    # and the other overrides still reach stage 3 in its one --config.
    assert all('--skip-virus' not in c for c in calls), calls
    assert calls[2].count('--config') == 1, calls[2]
    assert 'skip_cmseq=True' in calls[2], calls[2]
    assert calls[2].rstrip().endswith('binning=output/config/auto_binning.txt'), calls[2]


def planned(*targets):
    """Job counts from a Snakemake dry run of Snakefile-bin on the test data."""
    with tempfile.TemporaryDirectory(prefix='magmaker_virus_plan_') as tmp:
        # metadata.txt names its reads relative to the repository root.
        os.symlink(os.path.join(REPO_ROOT, 'resources'),
                   os.path.join(tmp, 'resources'))
        # A tiny host genome of its own, not BoarderCollie, whose default
        # paths exist only on the cluster.
        with open(os.path.join(tmp, 'host.fna'), 'w') as handle:
            handle.write(">c1\nACGTACGTACGGTTAACCGGTTAAGGCCAATTGGCCAAGGTTCCAAGG\n")
        over = os.path.join(tmp, 'over.yaml')
        with open(over, 'w') as handle:
            handle.write(f"host_filter:\n  boardercollie: false\n"
                         f"  db_dir: {tmp}/host_db\n  genome: {tmp}/host.fna\n")
        result = subprocess.run(
            ['snakemake', '-s', os.path.join(REPO_ROOT, 'Snakefile-bin'),
             '--configfile', os.path.join(REPO_ROOT, 'config.yaml'), over,
             '--cores', '1', '-d', tmp, '-n', *targets,
             '--config', 'binning=resources/config/binning.txt'],
            capture_output=True, text=True)
        out = result.stdout + result.stderr
        assert result.returncode == 0, out[-3000:]
        stats = out[out.index('Job stats'):]
        return {m.group(1): int(m.group(2))
                for m in re.finditer(r'^(\w+)\s+(\d+)$', stats, re.M)}


def test_plan():
    both = planned('rename_mags', 'virus_all')
    n = len(TEST_SAMPLES)
    assert both.get('run_genomad') == n, both
    assert both.get('run_checkv') == n, both
    assert both.get('make_virus_summary') == 1, both
    assert both.get('rename_mags') == 1, both
    mags_only = planned('rename_mags')
    assert 'run_genomad' not in mags_only and 'run_checkv' not in mags_only, mags_only
    # The viral track adds its own jobs and nothing else.
    extra = {k: both[k] - mags_only.get(k, 0) for k in both if both[k] != mags_only.get(k, 0)}
    assert set(extra) == {'run_genomad', 'run_checkv', 'make_virus_summary', 'virus_all',
                          'total'}, extra
    return both, mags_only


if __name__ == '__main__':
    test_wrapper_targets()
    both, mags_only = test_plan()
    print(f"plan: rename_mags + virus_all = {both['total']} jobs "
          f"(run_genomad {both['run_genomad']}, run_checkv {both['run_checkv']}); "
          f"rename_mags alone = {mags_only['total']} jobs")
    print('test_run_magmaker: 2 passed')
