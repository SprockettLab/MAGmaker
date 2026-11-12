#!/usr/bin/env python3
"""Copy a FASTA, making repeated sequence IDs unique, and report the counts.

CheckV refuses a FASTA with a repeated ID ("Please remove duplicated sequence
IDs from the input FASTA file"). VIBRANT produces them whenever one contig
yields several phage regions, which only shows up on large assemblies (WF30
saw it on 2 of 4 test samples, on the two largest).

A repeat is RENAMED to <id>__dup<N>, never dropped: each one is a distinct
predicted region, and dropping one would lose a call silently. The summary's
name mapping strips __dup<N> again, so the source contig is unchanged.

Writes <out>.counts.tsv beside the output with sequences_in, unique_ids_in and
sequences_out, so a later table can be checked against what went in.

Usage: unique_fasta_ids.py <in.fna> <out.fna>
An input that does not exist, or is empty, gives an empty output and counts of
zero: a predictor that called nothing is not an error.
"""

import os
import sys


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, dest = sys.argv[1], sys.argv[2]
    os.makedirs(os.path.dirname(os.path.abspath(dest)) or '.', exist_ok=True)

    seen, n_in, n_out = {}, 0, 0
    with open(dest, 'w') as out:
        if os.path.exists(src) and os.path.getsize(src) > 0:
            with open(src) as handle:
                for line in handle:
                    if not line.startswith('>'):
                        out.write(line)
                        continue
                    n_in += 1
                    n_out += 1
                    name = line[1:].split()[0] if line[1:].split() else 'unnamed'
                    seen[name] = seen.get(name, 0) + 1
                    if seen[name] > 1:
                        out.write('>{0}__dup{1}\n'.format(name, seen[name]))
                    else:
                        out.write(line)

    with open(dest + '.counts.tsv', 'w') as handle:
        handle.write('sequences_in\tunique_ids_in\tsequences_out\n')
        handle.write('{0}\t{1}\t{2}\n'.format(n_in, len(seen), n_out))

    if n_in != n_out:
        sys.exit('ERROR: {0} sequences in, {1} written'.format(n_in, n_out))
    print('{0} sequences, {1} unique ids, {2} written to {3}'.format(
        n_in, len(seen), n_out, dest))


if __name__ == '__main__':
    main()
