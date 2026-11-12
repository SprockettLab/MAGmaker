import pandas as pd
from os.path import join

configfile: "config.yaml"

# Loads the metadata table and defines metadata_table, samples, seqruns,
# reads, get_read() and seqruns_for(). Shared with Snakefile so the two
# stages cannot disagree about the sample sheet.
include: "resources/snakefiles/common.smk"

# Sample and sequencing-run identifiers never contain a path separator, but
# Snakemake wildcards match across "/" unless they are constrained. Without
# this, output/qc/fastp/se/<sample>.<run>.fastp.json is matchable by fastp_pe
# by binding sample="se/<sample>", so the paired rule is selected for a
# single-end report and dies in its input function with a KeyError naming a
# sample that was never in the metadata. Constraining the two wildcards the
# se/ subdirectory can be absorbed into sends each report to its own rule.
wildcard_constraints:
    sample=r"[^/]+",
    seqrun=r"[^/]+"

binning_fp = config['binning']

binning_df = pd.read_csv(binning_fp,
                         header=0,
                         index_col=0,
                         sep='\t',
                         na_filter=False)

def parse_groups(group_series):
    groups = {}
    for sample, grps in group_series.items():
        if not grps:
            continue
        grp_list = grps.split(',')
        for grp in grp_list:
            if grp not in groups:
                groups[grp] = [sample]
            else:
                groups[grp].append(sample)
    return(groups)

def make_pairings(read_grp, ctg_grp):
    if read_grp.keys() != ctg_grp.keys():
        raise ValueError('Not all keys in both from and to groups!')

    pairings = []
    contig_pairings = {}
    for grp in read_grp.keys():
        r = read_grp[grp]
        c = ctg_grp[grp]

        for i in r:
            for  j in c:
                pairings.append((i, j))
                if j not in contig_pairings:
                    contig_pairings[j] = [i]
                else:
                    contig_pairings[j].append(i)

    return(pairings, contig_pairings)

contig_groups = parse_groups(binning_df['Contig_Groups'])
read_groups = parse_groups(binning_df['Read_Groups'])
pairings, contig_pairings = make_pairings(read_groups, contig_groups)

print('Contig samples: %s' % contig_groups)
print('Read samples: %s' % read_groups)
print('Pairings: %s' % pairings)
print('Contig Pairings: %s' % contig_pairings)

def get_contigs(sample, binning_df):
    return(binning_df.loc[sample, 'Contigs'])


# A binning row is a contig set. Short-read rows are named after their sample;
# long-read rows are <sample>__<assembler> (generate_binning_config). The
# optional Source_Sample and Assembler columns say which sample and assembler a row
# is; a table without them (older, or hand-written) keeps the old meaning:
# the row is the sample, assembled by the first configured assembler that can.
def binning_column(row, column):
    if column in binning_df.columns:
        value = str(binning_df.loc[row, column]).strip()
        if value:
            return value
    return None


def sample_of(row):
    """The sample whose reads a binning row's contigs were assembled from."""
    return binning_column(row, 'Source_Sample') or row


def assembler_of(row):
    """The assembler that made a binning row's contigs."""
    named = binning_column(row, 'Assembler')
    if named:
        return named
    for assembler in config['assemblers']:
        if sample_of(row) in set(samples_for_assembler(assembler)):
            return assembler
    return 'unknown'


def is_long_read_row(row):
    return assembler_of(row) in LONG_READ_ASSEMBLER_NAMES


def row_contigs(row):
    """The contigs FASTA of a binning row, as the table names it."""
    return get_contigs(row, binning_df)

include: "resources/snakefiles/qc.smk"
include: "resources/snakefiles/long_reads.smk"
include: "resources/snakefiles/assemble.smk"
include: "resources/snakefiles/mapping.smk"
include: "resources/snakefiles/binning.smk"
include: "resources/snakefiles/selected_bins.smk"
include: "resources/snakefiles/cmseq.smk"
include: "resources/snakefiles/mag_qc.smk"
include: "resources/snakefiles/virus.smk"
include: "resources/snakefiles/long_read_binning.smk"


rule select_bins:
    input:
        lambda wildcards: expand("output/selected_bins/{mapper}/" + SELECTED_FASTAS + "/{contig_sample}/.done",
                                 mapper=config['mappers'],
                                 contig_sample=contig_pairings.keys())

rule bin_all:
    # Only the binners listed in `binners:`. That key existed and was read
    # by nothing: this target asked for all three unconditionally, and both
    # consolidation tools took their bin sets from a hardcoded list, so
    # removing a binner from the config changed nothing at all.
    input:
        lambda wildcards: binner_bin_dirs(list(contig_pairings.keys()))

rule map_all:
    input:
        expand("output/mapping/{mapper}/sorted_bams/{pairing[0]}_Mapped_To_{pairing[1]}.bam",
               mapper=config['mappers'],
               pairing=pairings)
