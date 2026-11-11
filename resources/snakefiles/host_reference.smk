import copy
import os
import shutil
import subprocess
import sys
from os.path import abspath, basename, exists, join, normpath, splitext

# Host reference: the FASTA host_filter maps reads against, and its bowtie2
# index. Two ways, chosen by host_filter.boardercollie:
#
#   true  (default)  The host's genomes are listed in the boardercollie
#                    section. BoarderCollie, loaded as a Snakemake module,
#                    downloads, screens, cleans and joins them. The joined
#                    reference, and MAGmaker's index of it, go to
#                    output/qc/host_reference/. host_filter.genome and
#                    host_filter.db_dir are not used, so older project
#                    configs that still set them run with BoarderCollie too.
#   false            host_filter.genome is an existing FASTA, used as it is
#                    with its index in db_dir, exactly as before BoarderCollie.
#
# Defines host_reference (the FASTA), host_base (the bowtie2 index stem) and
# rule host_bowtie2_build.

# BoarderCollie mode: the joined reference and its index are per project.
HOST_REFERENCE_DIR = "output/qc/host_reference"

# BoarderCollie is fetched by MAGmaker itself: cloned on the first run into
# resources/external/BoarderCollie (git-ignored) at the version below, and
# moved to that version whenever MAGmaker pins a newer one. To use a newer
# BoarderCollie, change BOARDERCOLLIE_VERSION to its commit and test.
BOARDERCOLLIE_URL = "https://github.com/SprockettLab/BoarderCollie.git"
BOARDERCOLLIE_VERSION = "a393e862693c0e367fb86c5db3524053d9c552df"
BOARDERCOLLIE_DEFAULT_DIR = normpath(join(workflow.current_basedir, '..',
                                          'external', 'BoarderCollie'))
# Kept inside .git, so it never shows as a change in the clone: marks a
# clone MAGmaker made and may move to another version. A clone without it
# (one named in host_filter.boardercollie_dir) is used as it is.
MANAGED_MARK = join('.git', 'magmaker-managed')


def _git(*args, cwd=None):
    """Run git; raise ValueError with its message on failure. A clone that
    needs a login fails at once instead of waiting for a password."""
    env = {**os.environ, 'GIT_TERMINAL_PROMPT': '0'}
    result = subprocess.run(['git', *args], cwd=cwd, env=env,
                            capture_output=True, text=True)
    if result.returncode != 0:
        raise ValueError("git " + " ".join(args) + " failed:\n"
                         + (result.stderr or result.stdout).strip())
    return result.stdout.strip()


def ensure_boardercollie(clone, url, version):
    """in : clone folder, repository URL, commit
       out: nothing; afterwards a MAGmaker-made clone is at that commit

    A missing clone is made in a folder beside it and renamed into place,
    so a run that is interrupted, or two runs starting at once, never leave
    half a clone where the next run would take it for a whole one. A clone
    without MAGmaker's mark is never changed."""
    if exists(clone):
        if not exists(join(clone, MANAGED_MARK)):
            return
        if _git('rev-parse', 'HEAD', cwd=clone) == version:
            return
        if _git('status', '--porcelain', '--untracked-files=no', cwd=clone):
            raise ValueError(
                f"{clone} has local changes, so MAGmaker will not move it to "
                f"BoarderCollie {version}. Remove the folder and run again; "
                "it is cloned again.")
        print(f"MAGmaker: moving BoarderCollie in {clone} to {version}",
              file=sys.stderr)
        _git('fetch', '--quiet', 'origin', cwd=clone)
        _git('checkout', '--quiet', '--detach', version, cwd=clone)
        return
    print(f"MAGmaker: cloning BoarderCollie {version} into {clone}",
          file=sys.stderr)
    os.makedirs(os.path.dirname(clone), exist_ok=True)
    partial = f"{clone}.partial.{os.getpid()}"
    try:
        _git('clone', '--quiet', url, partial)
        _git('checkout', '--quiet', '--detach', version, cwd=partial)
        open(join(partial, MANAGED_MARK), 'w').close()
        try:
            os.rename(partial, clone)
        except OSError:
            # Another run made the clone first; use that one.
            if not exists(join(clone, MANAGED_MARK)):
                raise
    finally:
        shutil.rmtree(partial, ignore_errors=True)


def _host_fail(message):
    """Abort with a message readable without a Python traceback."""
    sys.exit("\nMAGmaker host reference error\n" + "-" * 60 + "\n"
             + message + "\n")


_hf = config['host_filter']
USE_BOARDERCOLLIE = _hf.get('boardercollie', True)
if not isinstance(USE_BOARDERCOLLIE, bool):
    _host_fail("host_filter.boardercollie must be true or false, got "
               f"{USE_BOARDERCOLLIE!r}.")


def boardercollie_setup(hf, bc_config):
    """in : the host_filter and boardercollie config sections
       out: (reference path, index stem, Snakefile path, module config)

    The module config is a copy with reference_out filled in; the caller's
    config is not changed."""
    clone = hf.get('boardercollie_dir') or BOARDERCOLLIE_DEFAULT_DIR
    ensure_boardercollie(clone,
                         hf.get('boardercollie_url') or BOARDERCOLLIE_URL,
                         hf.get('boardercollie_version')
                         or BOARDERCOLLIE_VERSION)
    snakefile = join(clone, 'workflow', 'Snakefile')
    if not exists(snakefile):
        raise ValueError(
            f"No BoarderCollie workflow at {snakefile}.\n"
            "host_filter.boardercollie_dir must be a clone of "
            "github.com/SprockettLab/BoarderCollie, or left out so that "
            "MAGmaker fetches BoarderCollie itself.")

    bc = copy.deepcopy(bc_config or {})
    hosts = bc.get('hosts') or []
    if len(hosts) != 1:
        raise ValueError(
            "boardercollie.hosts must list exactly one host, the one reads "
            f"are filtered against; it lists {len(hosts)}.")
    host = hosts[0]
    name = str(host.get('name') or '')
    if not name or '/' in name:
        raise ValueError("boardercollie.hosts: the host needs a name without "
                         f"'/', got {name!r}.")
    if host.get('reference_out'):
        raise ValueError(
            "boardercollie.hosts: leave reference_out out. MAGmaker sets it "
            f"to {HOST_REFERENCE_DIR}/<host name>.fna.")

    reference = normpath(join(HOST_REFERENCE_DIR, name + '.fna'))

    # The original FASTA is never changed. BoarderCollie writes the joined
    # reference with ln -f (replacing whatever is at the path), and
    # Snakemake deletes a job's outputs before it runs, so the path must not
    # be one of the host's own FASTA files, nor any file BoarderCollie did
    # not write (a reference it wrote always has members.tsv beside it).
    for genome in host.get('genomes') or []:
        if isinstance(genome, dict) and genome.get('fasta') and \
                abspath(genome['fasta']) == abspath(reference):
            raise ValueError(
                f"The joined reference would be written to {reference}, "
                "which is one of the host's own genome FASTA files. Move "
                "that FASTA elsewhere.")
    if exists(reference) and not exists(reference + '.members.tsv'):
        raise ValueError(
            f"{reference} exists but was not written by BoarderCollie (no "
            f"{basename(reference)}.members.tsv beside it). MAGmaker will not "
            "overwrite it. Move it elsewhere, or use it as it is with  "
            "host_filter.boardercollie: false  and  host_filter.genome: "
            + reference)

    host['reference_out'] = reference
    return (reference, normpath(join(HOST_REFERENCE_DIR, name)), snakefile,
            bc)


def boardercollie_summary(hf, bc):
    """Which host reference this run uses, printed when the workflow is read."""
    host = bc['hosts'][0]
    genomes = [g if isinstance(g, str) else
               (g.get('accession') or g.get('fasta'))
               for g in host.get('genomes') or []]
    lines = ["MAGmaker host reference: BoarderCollie, host " + host['name']
             + ", genomes " + ", ".join(map(str, genomes))
             + ", written to " + host['reference_out']]
    ignored = [k for k in ('genome', 'db_dir') if hf.get(k)]
    if ignored:
        lines.append("  (host_filter." + " and host_filter.".join(ignored)
                     + " not used while host_filter.boardercollie is true)")
    return "\n".join(lines)


if USE_BOARDERCOLLIE:
    try:
        (host_reference, host_base, BOARDERCOLLIE_SNAKEFILE,
         BOARDERCOLLIE_CONFIG) = boardercollie_setup(
            _hf, config.get('boardercollie'))
    except ValueError as err:
        _host_fail(str(err))
    print(boardercollie_summary(_hf, BOARDERCOLLIE_CONFIG), file=sys.stderr)
    include: "host_index_boardercollie.smk"
else:
    if not _hf.get('genome') or not _hf.get('db_dir'):
        _host_fail("host_filter.boardercollie is false, so host_filter.genome "
                   "must name an existing host FASTA and host_filter.db_dir "
                   "the directory its bowtie2 index is kept in.")
    host_reference = _hf['genome']
    host_base = join(_hf['db_dir'], splitext(basename(_hf['genome']))[0])
    include: "host_index_prebuilt.smk"
