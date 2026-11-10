import copy
import os
import sys
from os.path import abspath, basename, exists, join, normpath, splitext

# Host reference: the FASTA host_filter maps reads against, and its bowtie2
# index. Two ways, chosen by host_filter.boardercollie:
#
#   true  (default)  The host's genomes are listed in the boardercollie
#                    section. BoarderCollie, loaded as a Snakemake module,
#                    downloads, screens, cleans and joins them, and writes
#                    the joined reference to <db_dir>/<host name>.fna.
#                    host_filter.genome must be empty.
#   false            host_filter.genome is an existing FASTA, used as it is
#                    with its index in db_dir, exactly as before BoarderCollie.
#
# Defines host_reference (the FASTA), host_base (the bowtie2 index stem) and
# rule host_bowtie2_build.


def _host_fail(message):
    """Abort with a message readable without a Python traceback."""
    sys.exit("\nMAGmaker host reference error\n" + "-" * 60 + "\n"
             + message + "\n")


_hf = config['host_filter']
USE_BOARDERCOLLIE = _hf.get('boardercollie', True)
if not isinstance(USE_BOARDERCOLLIE, bool):
    _host_fail("host_filter.boardercollie must be true or false, got "
               f"{USE_BOARDERCOLLIE!r}.")
if not _hf.get('db_dir'):
    _host_fail("host_filter.db_dir is not set. It is the directory the "
               "bowtie2 index is kept in.")


def boardercollie_setup(hf, bc_config):
    """in : the host_filter and boardercollie config sections
       out: (reference path, index stem, Snakefile path, module config)

    Plain Python, so the checks can be tested without Snakemake. The
    module config is a copy with reference_out filled in; the caller's
    config is not changed."""
    if hf.get('genome'):
        raise ValueError(
            "host_filter.genome is set (" + str(hf['genome']) + ") while "
            "host_filter.boardercollie is true.\n\n"
            "To use that FASTA as it is, set  host_filter.boardercollie: false.\n"
            "To have BoarderCollie build the host reference, leave "
            "host_filter.genome empty and list the host's genomes under "
            "boardercollie.hosts.")
    clone = hf.get('boardercollie_dir')
    if not clone:
        raise ValueError(
            "host_filter.boardercollie_dir is not set. It is the path of a "
            "local clone of BoarderCollie.")
    snakefile = join(clone, 'workflow', 'Snakefile')
    if not exists(snakefile):
        raise ValueError(
            f"No BoarderCollie workflow at {snakefile}.\n"
            "host_filter.boardercollie_dir must be a clone of "
            "github.com/SprockettLab/BoarderCollie.")

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
            "to <host_filter.db_dir>/<host name>.fna.")

    reference = normpath(join(hf['db_dir'], name + '.fna'))

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
                "which is one of the host's own genome FASTA files. Set "
                "host_filter.db_dir to another directory.")
    if exists(reference) and not exists(reference + '.members.tsv'):
        raise ValueError(
            f"{reference} exists but was not written by BoarderCollie (no "
            f"{basename(reference)}.members.tsv beside it). MAGmaker will not "
            "overwrite it. Set host_filter.db_dir to another directory, or "
            "use the file as it is with  host_filter.boardercollie: false  "
            "and  host_filter.genome: " + reference)

    host['reference_out'] = reference
    return reference, normpath(join(hf['db_dir'], name)), snakefile, bc


if USE_BOARDERCOLLIE:
    try:
        (host_reference, host_base, BOARDERCOLLIE_SNAKEFILE,
         BOARDERCOLLIE_CONFIG) = boardercollie_setup(
            _hf, config.get('boardercollie'))
    except ValueError as err:
        _host_fail(str(err))
    include: "host_index_boardercollie.smk"
else:
    if not _hf.get('genome'):
        _host_fail("host_filter.boardercollie is false, so host_filter.genome "
                   "must name an existing host FASTA.")
    host_reference = _hf['genome']
    host_base = join(_hf['db_dir'], splitext(basename(_hf['genome']))[0])
    include: "host_index_prebuilt.smk"
