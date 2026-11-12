#!/usr/bin/env python3
"""Tests the two ways of giving the host reference (host_reference.smk).

  boardercollie: true   BoarderCollie, loaded as a Snakemake module, joins
                        the host's genomes into
                        output/qc/host_reference/<host>.fna, and the bowtie2
                        index is tied to that file's md5. An older config's
                        host_filter.genome and db_dir are not used.
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


DONE_TARGET = os.path.join('output', 'qc', 'host_reference',
                           'H.bowtie2_build.done')
PARTS = ('1', '2', '3', '4', 'rev.1', 'rev.2')


def shared_indexes(tmp):
    """Folders of BoarderCollie's shared bowtie2 index for host H."""
    base = os.path.join(tmp, 'idx', 'bt2', 'H')
    return sorted(os.listdir(base)) if os.path.isdir(base) else []


def linked_to(db_dir, tmp, key):
    """Every index part in db_dir is the same file as the shared one."""
    shared = os.path.join(tmp, 'idx', 'bt2', 'H', key, 'H')
    return all(os.stat(f"{os.path.join(db_dir, 'H')}.{p}.bt2l").st_ino
               == os.stat(f"{shared}.{p}.bt2l").st_ino for p in PARTS)


def test_boardercollie_mode():
    """The bowtie2 index is BoarderCollie's shared one: built once per host
    key under index_dir, linked into each run folder. A second run folder
    with the same host builds nothing; a changed host gets its own index."""
    tmp = new_tmp()
    db_dir = os.path.join(tmp, 'output', 'qc', 'host_reference')
    ref = os.path.join(db_dir, 'H.fna')
    done = os.path.join(db_dir, 'H.bowtie2_build.done')
    host1 = os.path.join(tmp, 'host1.fna')
    host2 = os.path.join(tmp, 'host2.fna')
    old_dir = os.path.join(tmp, 'old_index')
    os.makedirs(old_dir)
    old_genome = os.path.join(old_dir, 'old_host.fna')
    with open(old_genome, 'w') as handle:
        handle.write(">old\nACGT\n")
    hf = {'boardercollie': True, 'boardercollie_dir': BC_DIR,
          'db_dir': old_dir, 'genome': old_genome}

    host1_md5, host1_mtime = md5(host1), os.stat(host1).st_mtime_ns
    write_config(tmp, hf, bc_section(
        tmp, [{'fasta': host1, 'tax_id': 9606}]))
    first = snakemake(tmp, DONE_TARGET)
    out = first.stdout + first.stderr
    check("BoarderCollie mode builds the reference and index",
          first.returncode == 0, out)
    if first.returncode != 0:
        return
    check("joined reference written to output/qc/host_reference/<host>.fna",
          read(ref) == read(host1))
    check("older config's genome and db_dir left untouched",
          sorted(os.listdir(old_dir)) == ['old_host.fna']
          and read(old_genome) == ">old\nACGT\n")
    check("the host and the shared index are printed",
          'MAGmaker host reference: BoarderCollie, host H' in first.stderr
          and 'not used while' in first.stderr
          and "BoarderCollie's shared bowtie2 index" in first.stderr,
          first.stderr)
    check("members.tsv written beside the reference",
          os.path.exists(ref + '.members.tsv'))
    keys = shared_indexes(tmp)
    check("the index is built once, in BoarderCollie's index_dir",
          len(keys) == 1 and 'rule boardercollie_bowtie2:' in out, keys)
    check("the run folder's index parts are links to the shared index",
          keys and linked_to(db_dir, tmp, keys[0]))
    check("sentinel names the reference and the shared index",
          keys and ref.endswith(read(done).split('\n')[0].split('\t')[1])
          and keys[0] in read(done), read(done))
    check("no index of the run folder's own (no md5 stamp)",
          not os.path.exists(os.path.join(db_dir, 'H.reference.md5')))
    check("original host FASTA unchanged",
          md5(host1) == host1_md5
          and os.stat(host1).st_mtime_ns == host1_mtime)

    again = snakemake(tmp, '-n', DONE_TARGET)
    check("nothing to redo on an unchanged config",
          again.returncode == 0
          and 'Nothing to be done' in again.stdout + again.stderr,
          again.stdout + again.stderr)

    # A second run folder, same host and same shared store: links, no build.
    other = new_tmp()
    write_config(other, hf, bc_section(
        tmp, [{'fasta': host1, 'tax_id': 9606}]))
    second_folder = snakemake(other, DONE_TARGET)
    out = second_folder.stdout + second_folder.stderr
    check("a second run folder with the same host builds no index",
          second_folder.returncode == 0
          and 'rule boardercollie_bowtie2:' not in out
          and 'rule boardercollie_bowtie2_out_H:' in out
          and shared_indexes(tmp) == keys, out)
    check("and links the same shared index",
          keys and linked_to(os.path.join(other, 'output', 'qc',
                                          'host_reference'), tmp, keys[0]))
    shutil.rmtree(other)

    # Other genomes: another host key, its own shared index, relinked here.
    write_config(tmp, hf, bc_section(
        tmp, [{'fasta': host2, 'tax_id': 9606}]))
    changed = snakemake(tmp, DONE_TARGET)
    check("changed host genomes run again", changed.returncode == 0,
          changed.stdout + changed.stderr)
    if changed.returncode != 0:
        return
    new_keys = shared_indexes(tmp)
    new_key = [k for k in new_keys if k not in keys]
    check("reference at the same path now holds the new genome",
          read(ref) == read(host2))
    check("the changed host has its own shared index; the first is kept",
          len(new_keys) == 2 and len(new_key) == 1, new_keys)
    check("the run folder now links the new index, not the old one",
          new_key and linked_to(db_dir, tmp, new_key[0]))
    check("first host FASTA still unchanged", md5(host1) == host1_md5)
    shutil.rmtree(tmp)


def test_own_index_kept():
    """A run folder that built its own index before shared indexes existed
    (stamped <host>.reference.md5) keeps it and the old rule: nothing
    reruns, nothing is built in the shared store."""
    tmp = new_tmp()
    db_dir = os.path.join(tmp, 'output', 'qc', 'host_reference')
    os.makedirs(db_dir)
    with open(os.path.join(db_dir, 'H.reference.md5'), 'w') as handle:
        handle.write("0" * 32 + "\n")
    hf = {'boardercollie': True, 'boardercollie_dir': BC_DIR}
    write_config(tmp, hf, bc_section(
        tmp, [{'fasta': os.path.join(tmp, 'host1.fna'), 'tax_id': 9606}]))
    run = snakemake(tmp, DONE_TARGET)
    out = run.stdout + run.stderr
    log = read(os.path.join(tmp, 'output', 'logs', 'qc', 'host_bowtie2_build',
                            'host_bowtie2_build.log')) if run.returncode == 0 else ''
    check("a folder with its own index keeps the old rule",
          run.returncode == 0 and 'kept as it is' in run.stderr
          and 'rule boardercollie_bowtie2' not in out, out)
    check("that rule (re)builds in the run folder when the stamp does not "
          "match, as before", log.startswith('Building'), log)
    check("nothing is built in the shared store", shared_indexes(tmp) == [])
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
    hf = {'boardercollie': True, 'boardercollie_dir': BC_DIR}
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
    def empty_dir(tmp):
        os.makedirs(os.path.join(tmp, 'not_a_clone'))
    expect_error("a folder that is not a BoarderCollie clone is refused",
                 "No BoarderCollie workflow",
                 lambda t: {'boardercollie_dir': os.path.join(t, 'not_a_clone')},
                 setup=empty_dir)
    expect_error("two hosts are refused", "exactly one host",
                 lambda t: {},
                 lambda t: {**bc_section(t, ['GCA_1.1']),
                            'hosts': [{'name': 'A', 'genomes': ['GCA_1.1']},
                                      {'name': 'B', 'genomes': ['GCA_2.1']}]})
    expect_error("reference_out given by hand is refused",
                 "leave reference_out out", lambda t: {},
                 lambda t: bc_section(t, ['GCA_1.1'],
                                      reference_out=os.path.join(t, 'x.fna')))
    expect_error("bowtie2_out given by hand is refused",
                 "leave bowtie2_out out", lambda t: {},
                 lambda t: bc_section(t, ['GCA_1.1'],
                                      bowtie2_out=os.path.join(t, 'x')))
    ref_dir = lambda t: os.path.join(t, 'output', 'qc', 'host_reference')

    def own_fasta(tmp):
        os.makedirs(ref_dir(tmp))
        shutil.copy(os.path.join(tmp, 'host1.fna'), ref_dir(tmp))
    tmp = expect_error("a host FASTA at the reference path is refused",
                       "one of the host's own genome FASTA files",
                       lambda t: {},
                       lambda t: {**bc_section(t, []), 'hosts': [{
                           'name': 'host1',
                           'genomes': [{'fasta': os.path.join(
                               ref_dir(t), 'host1.fna'), 'tax_id': 9606}]}]},
                       setup=own_fasta)
    check("and that FASTA is left as it was",
          read(os.path.join(ref_dir(tmp), 'host1.fna'))
          == read(os.path.join(tmp, 'host1.fna')))

    def plant(tmp):
        os.makedirs(ref_dir(tmp))
        with open(os.path.join(ref_dir(tmp), 'H.fna'), 'w') as handle:
            handle.write(">mine\nACGT\n")
    tmp = expect_error("a file BoarderCollie did not write is refused",
                       "was not written by BoarderCollie", lambda t: {},
                       setup=plant)
    check("and that file is left as it was",
          read(os.path.join(ref_dir(tmp), 'H.fna')) == ">mine\nACGT\n")


def git(*args, cwd=None):
    return subprocess.run(['git', *args], cwd=cwd, capture_output=True,
                          text=True, check=True).stdout.strip()


def test_automatic_clone():
    """MAGmaker clones BoarderCollie itself, at the version it pins, and
    moves its own clone when the pin changes; a clone it did not make is
    never changed. Cloned from the local BoarderCollie clone (BC_DIR) so the
    test needs no network."""
    tmp = new_tmp()
    clone = os.path.join(tmp, 'external', 'BoarderCollie')
    head = git('rev-parse', 'HEAD', cwd=BC_DIR)
    # An older commit that already has bowtie2_out (MAGmaker sets it): the
    # feature commit merged as HEAD's second parent.
    older = git('rev-parse', 'HEAD^2', cwd=BC_DIR)
    bc = bc_section(tmp, [{'fasta': os.path.join(tmp, 'host1.fna'),
                           'tax_id': 9606}])

    def dry(version):
        write_config(tmp, {'boardercollie': True, 'boardercollie_dir': clone,
                           'boardercollie_url': BC_DIR,
                           'boardercollie_version': version}, bc)
        return snakemake(tmp, '-n')

    first = dry(older)
    check("missing BoarderCollie is cloned at the pinned version",
          first.returncode == 0
          and git('rev-parse', 'HEAD', cwd=clone) == older
          and os.path.exists(os.path.join(clone, '.git',
                                          'magmaker-managed')),
          first.stdout + first.stderr)
    check("no partial clone left behind",
          os.listdir(os.path.dirname(clone)) == ['BoarderCollie'])
    check("the clone shows no local changes",
          git('status', '--porcelain', cwd=clone) == '')

    second = dry(head)
    check("MAGmaker's own clone moves to a newly pinned version",
          second.returncode == 0
          and git('rev-parse', 'HEAD', cwd=clone) == head
          and 'moving BoarderCollie' in second.stderr,
          second.stdout + second.stderr)

    with open(os.path.join(clone, 'README.md'), 'a') as handle:
        handle.write("local edit\n")
    third = dry(older)
    check("a clone with local changes is not moved",
          third.returncode != 0 and 'has local changes' in third.stderr
          and git('rev-parse', 'HEAD', cwd=clone) == head,
          third.stdout + third.stderr)

    bc_head = git('rev-parse', 'HEAD', cwd=BC_DIR)
    write_config(tmp, {'boardercollie': True, 'boardercollie_dir': BC_DIR,
                       'boardercollie_version': older}, bc)
    fourth = snakemake(tmp, '-n')
    check("a clone MAGmaker did not make is used as it is",
          fourth.returncode == 0
          and git('rev-parse', 'HEAD', cwd=BC_DIR) == bc_head,
          fourth.stdout + fourth.stderr)

    write_config(tmp, {'boardercollie': True,
                       'boardercollie_dir': os.path.join(tmp, 'other'),
                       'boardercollie_url': os.path.join(tmp, 'no_repo')}, bc)
    fifth = snakemake(tmp, '-n')
    check("a failed clone stops the run with git's message",
          fifth.returncode != 0 and 'git clone' in fifth.stderr
          and not os.path.exists(os.path.join(tmp, 'other')),
          fifth.stdout + fifth.stderr)
    shutil.rmtree(tmp)


def test_mitochondrion_passes_through():
    """A genome entry naming its species' mitochondrion (BoarderCollie's
    mitochondrion key) is planned through MAGmaker: fetched, then joined
    right after its genome."""
    tmp = new_tmp()
    host1 = os.path.join(tmp, 'host1.fna')
    hf = {'boardercollie': True, 'boardercollie_dir': BC_DIR}
    write_config(tmp, hf, bc_section(
        tmp, [{'fasta': host1, 'tax_id': 9606,
               'mitochondrion': 'NC_012920.1'}]))
    dry = snakemake(tmp, '-n', '-p', os.path.join(
        'output', 'qc', 'host_reference', 'H.bowtie2_build.done'))
    out = dry.stdout + dry.stderr
    check("a mitochondrion entry plans", dry.returncode == 0, out[-2000:])
    check("its fetch step is planned",
          'rule boardercollie_fetch_mitochondrion:' in out, out[-2000:])
    check("it is joined after its genome",
          '--accessions host1 NC_012920.1 ' in out, out[-2000:])


if __name__ == '__main__':
    missing = [tool for tool in ('snakemake', 'bowtie2-build')
               if shutil.which(tool) is None]
    if missing or not BC_DIR:
        sys.exit("Needs " + ", ".join(missing + ([] if BC_DIR else
                 ['BOARDERCOLLIE_DIR set to a BoarderCollie clone'])))
    test_boardercollie_mode()
    test_own_index_kept()
    test_prebuilt_mode()
    test_config_errors()
    test_automatic_clone()
    test_mitochondrion_passes_through()
    print("\n{} failure(s)".format(len(FAILURES)))
    sys.exit(1 if FAILURES else 0)
