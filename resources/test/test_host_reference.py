#!/usr/bin/env python3
"""Tests the two ways of giving the host reference (host_reference.smk).

  boardercollie: true   BoarderCollie, loaded as a Snakemake module, joins
                        the host's genomes into <db_dir>/<host>.fna, and the
                        bowtie2 index is tied to that file's md5.
  boardercollie: false  host_filter.genome is used as it is.

Runs the real Snakefile on two tiny FASTA files with BoarderCollie's
download_only setting (no network, no FCS-GX). Needs snakemake, pandas,
scikit-bio and bowtie2-build on PATH, plus the python packages BoarderCollie
imports (pyyaml, edlib), and a BoarderCollie clone:

    BOARDERCOLLIE_DIR=/path/to/BoarderCollie \\
        python3 resources/test/test_host_reference.py

The point most worth testing: BoarderCollie writes the joined reference to
the same path whatever it contains, so an index sitting next to it must be
rebuilt when the contents change, and must not be when they do not.
"""

import hashlib
import os
import shutil
import subprocess
import sys
import tempfile

import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
BC_DIR = os.environ.get('BOARDERCOLLIE_DIR')

FAILURES = []


def check(name, condition, detail=''):
    print(("PASS  " if condition else "FAIL  ") + name)
    if not condition:
        FAILURES.append(name)
        if detail:
            print(detail[-3000:])


def md5(path):
    with open(path, 'rb') as handle:
        return hashlib.md5(handle.read()).hexdigest()


def read(path):
    with open(path) as handle:
        return handle.read()


def write_config(tmp, host_filter, boardercollie=None):
    """Override config on top of the repository's config.yaml."""
    config = {
        'metadata': os.path.join(REPO_ROOT, 'resources', 'config',
                                 'metadata.txt'),
        'host_filter': host_filter,
    }
    if boardercollie is not None:
        config['boardercollie'] = boardercollie
    path = os.path.join(tmp, 'over.yaml')
    with open(path, 'w') as handle:
        yaml.safe_dump(config, handle)
    return path


def bc_section(tmp, genomes, **host):
    return {
        'library_dir': os.path.join(tmp, 'lib'),
        'index_dir': os.path.join(tmp, 'idx'),
        'work_dir': os.path.join(tmp, 'work'),
        'hosts': [{'name': 'H', 'genomes': genomes, 'indexes': [], **host}],
        'settings': {'download_only': True},
    }


def snakemake(tmp, *args):
    env = {**os.environ, 'BOARDERCOLLIE_RUN_ID': 'test'}
    return subprocess.run(
        ['snakemake', '-s', os.path.join(REPO_ROOT, 'Snakefile'),
         '--configfile', os.path.join(REPO_ROOT, 'config.yaml'),
         os.path.join(tmp, 'over.yaml'), '--cores', '1', '-d', tmp, *args],
        env=env, capture_output=True, text=True)


def new_tmp():
    tmp = tempfile.mkdtemp(prefix='magmaker_host_')
    # metadata.txt names its reads relative to the repository root.
    os.symlink(os.path.join(REPO_ROOT, 'resources'),
               os.path.join(tmp, 'resources'))
    with open(os.path.join(tmp, 'host1.fna'), 'w') as handle:
        handle.write(">c1\nACGTACGTACGGTTAACCGGTTAAGGCCAATTGGCCAAGGTTCCAAGG\n"
                     ">c2\nGGGCCCAAATTTGGGCCCAAATTTAACCGGTTAACCGGTTAAGGCCTT\n")
    with open(os.path.join(tmp, 'host2.fna'), 'w') as handle:
        handle.write(">d1\nTTGACCATGGCATTACGGATCCAGTTACGATCCGGATACCAGGTTCA\n"
                     ">d2\nCAGTTGACCAGTAGGCATTACCGATGGACCATTAGCAGTACCGGATT\n")
    return tmp


def test_boardercollie_mode():
    tmp = new_tmp()
    db_dir = os.path.join(tmp, 'hostref')
    ref = os.path.join(db_dir, 'H.fna')
    done = os.path.join(db_dir, 'H.bowtie2_build.done')
    stamp = os.path.join(db_dir, 'H.reference.md5')
    host1 = os.path.join(tmp, 'host1.fna')
    host2 = os.path.join(tmp, 'host2.fna')
    hf = {'boardercollie': True, 'boardercollie_dir': BC_DIR,
          'db_dir': db_dir, 'genome': ''}

    host1_md5, host1_mtime = md5(host1), os.stat(host1).st_mtime_ns
    write_config(tmp, hf, bc_section(
        tmp, [{'fasta': host1, 'tax_id': 9606}]))
    first = snakemake(tmp, done)
    check("BoarderCollie mode builds the reference and index",
          first.returncode == 0, first.stdout + first.stderr)
    if first.returncode != 0:
        return
    check("joined reference written to <db_dir>/<host>.fna",
          read(ref) == read(host1))
    check("members.tsv written beside the reference",
          os.path.exists(ref + '.members.tsv'))
    check("index stamped with the reference's md5",
          read(stamp).strip() == host1_md5)
    check("sentinel records the md5", host1_md5 in read(done))
    check("original host FASTA unchanged",
          md5(host1) == host1_md5
          and os.stat(host1).st_mtime_ns == host1_mtime)
    check("BoarderCollie's own indexes are not built",
          not os.path.exists(os.path.join(tmp, 'idx')))

    again = snakemake(tmp, '-n', done)
    check("nothing to redo on an unchanged config",
          again.returncode == 0
          and 'Nothing to be done' in again.stdout + again.stderr,
          again.stdout + again.stderr)

    # Same reference path, other contents: the old index must not be kept.
    old_index = md5(os.path.join(db_dir, 'H.1.bt2'))
    write_config(tmp, hf, bc_section(
        tmp, [{'fasta': host2, 'tax_id': 9606}]))
    second = snakemake(tmp, done)
    check("changed host genomes run again", second.returncode == 0,
          second.stdout + second.stderr)
    if second.returncode != 0:
        return
    check("reference at the same path now holds the new genome",
          read(ref) == read(host2))
    check("index rebuilt for the new reference",
          read(stamp).strip() == md5(host2)
          and md5(os.path.join(db_dir, 'H.1.bt2')) != old_index)
    log = read(os.path.join(tmp, 'output', 'logs', 'qc',
                            'host_bowtie2_build', 'host_bowtie2_build.log'))
    check("rebuild logged", log.startswith('Building'), log)
    check("first host FASTA still unchanged", md5(host1) == host1_md5)

    # A complete index whose stamp matches is reused, not rebuilt.
    os.remove(done)
    index_mtime = os.stat(os.path.join(db_dir, 'H.1.bt2')).st_mtime_ns
    third = snakemake(tmp, done)
    log = read(os.path.join(tmp, 'output', 'logs', 'qc',
                            'host_bowtie2_build', 'host_bowtie2_build.log'))
    check("matching index reused",
          third.returncode == 0 and 'reused' in log
          and os.stat(os.path.join(db_dir, 'H.1.bt2')).st_mtime_ns
          == index_mtime, third.stdout + third.stderr + log)

    shutil.rmtree(tmp)


def test_prebuilt_mode():
    tmp = new_tmp()
    db_dir = os.path.join(tmp, 'idx_prebuilt')
    genome = os.path.join(tmp, 'host1.fna')
    write_config(tmp, {'boardercollie': False, 'db_dir': db_dir,
                       'genome': genome})
    dry = snakemake(tmp, '-n', os.path.join(db_dir, 'host1.bowtie2_build.done'))
    out = dry.stdout + dry.stderr
    check("boardercollie: false plans the index from host_filter.genome",
          dry.returncode == 0 and 'host_bowtie2_build' in out, out)
    check("boardercollie: false loads no BoarderCollie steps",
          'boardercollie_' not in out, out)
    shutil.rmtree(tmp)


def expect_error(name, message, host_filter, boardercollie=None,
                 setup=None):
    tmp = new_tmp()
    if setup:
        setup(tmp)
    hf = {'boardercollie': True, 'boardercollie_dir': BC_DIR,
          'db_dir': os.path.join(tmp, 'hostref'), 'genome': ''}
    hf.update(host_filter(tmp))
    bc = (boardercollie(tmp) if boardercollie else
          bc_section(tmp, [{'fasta': os.path.join(tmp, 'host1.fna'),
                            'tax_id': 9606}]))
    write_config(tmp, hf, bc)
    result = snakemake(tmp, '-n')
    out = result.stdout + result.stderr
    check(name, result.returncode != 0 and message in out, out)
    return tmp


def test_config_errors():
    expect_error("genome set while BoarderCollie is on is refused",
                 "host_filter.genome is set",
                 lambda t: {'genome': os.path.join(t, 'host1.fna')})
    expect_error("missing BoarderCollie clone is refused",
                 "No BoarderCollie workflow",
                 lambda t: {'boardercollie_dir': os.path.join(t, 'nowhere')})
    expect_error("two hosts are refused", "exactly one host",
                 lambda t: {},
                 lambda t: {**bc_section(t, ['GCA_1.1']),
                            'hosts': [{'name': 'A', 'genomes': ['GCA_1.1']},
                                      {'name': 'B', 'genomes': ['GCA_2.1']}]})
    expect_error("reference_out given by hand is refused",
                 "leave reference_out out", lambda t: {},
                 lambda t: bc_section(t, ['GCA_1.1'],
                                      reference_out=os.path.join(t, 'x.fna')))
    expect_error("a host FASTA at the reference path is refused",
                 "one of the host's own genome FASTA files",
                 lambda t: {'db_dir': t},
                 lambda t: {**bc_section(t, []), 'hosts': [{
                     'name': 'host1',
                     'genomes': [{'fasta': os.path.join(t, 'host1.fna'),
                                  'tax_id': 9606}]}]})

    def plant(tmp):
        os.makedirs(os.path.join(tmp, 'hostref'))
        with open(os.path.join(tmp, 'hostref', 'H.fna'), 'w') as handle:
            handle.write(">mine\nACGT\n")
    tmp = expect_error("a file BoarderCollie did not write is refused",
                       "was not written by BoarderCollie", lambda t: {},
                       setup=plant)
    check("and that file is left as it was",
          read(os.path.join(tmp, 'hostref', 'H.fna')) == ">mine\nACGT\n")


if __name__ == '__main__':
    missing = [tool for tool in ('snakemake', 'bowtie2-build')
               if shutil.which(tool) is None]
    if missing or not BC_DIR:
        sys.exit("Needs " + ", ".join(missing + ([] if BC_DIR else
                 ['BOARDERCOLLIE_DIR set to a BoarderCollie clone'])))
    test_boardercollie_mode()
    test_prebuilt_mode()
    test_config_errors()
    print("\n{} failure(s)".format(len(FAILURES)))
    sys.exit(1 if FAILURES else 0)
