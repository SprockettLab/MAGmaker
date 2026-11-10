import copy
import os
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
