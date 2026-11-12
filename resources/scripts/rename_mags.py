"""
Reads mag_summary.tsv and writes each MAG into renamed_mags/ as <New_Name>.fa.
Edit New_Name in the table and re-run this rule to update file names.

Every contig is renamed <MAG_ID>__<original contig name>, so contig names are
unique across all MAGs. Assemblers reuse the same names in every sample
(MEGAHIT's k141_1, Flye's contig_1), so without this two MAGs from different
samples can share contig names, which breaks any step that pools MAGs into one
reference (a dRep catalog, read mapping against the catalog). MAG_ID is used
rather than New_Name because it is short, unique by construction, and does not
change when New_Name is edited. Only the copies are renamed: the MAG FASTAs
under selected_bins/ are not touched. contig_map.tsv records every new name
with the original name and source file.
"""
import os
import sys
import shutil
import pandas as pd

log_path = str(snakemake.log[0])
sys.stderr = open(log_path, 'w')
sys.stdout = sys.stderr

summary_path = str(snakemake.input.summary)
renamed_dir = snakemake.params.renamed_dir

df = pd.read_csv(summary_path, sep='\t')

required = {'MAG_ID', 'New_Name', 'Original_Path'}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"mag_summary.tsv is missing required columns: {missing}")

# Clear destination so stale files from previous runs don't accumulate
if os.path.isdir(renamed_dir):
    shutil.rmtree(renamed_dir)
os.makedirs(renamed_dir)

seen = {}          # new contig name -> MAG_ID that holds it
contig_rows = []
n_renamed = 0
for _, row in df.iterrows():
    # Rows with MAG_ID NONE record a sample that produced no MAGs. They
    # carry no FASTA by construction, so they are skipped without the
    # missing-source warning, which would otherwise fire once per such
    # sample and read like a fault.
    mag_id = str(row['MAG_ID'])
    if mag_id == 'NONE':
        continue
    src = str(row['Original_Path'])
    new_name = str(row['New_Name'])
    dst = os.path.join(renamed_dir, f"{new_name}.fa")
    if not os.path.exists(src):
        print(f"Warning: source FASTA not found, skipping: {src}")
        continue
    if os.path.exists(dst):
        raise ValueError(f"two MAGs would both be written to {dst}; New_Name must be unique")
    with open(src) as fin, open(dst, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                original, _, description = line[1:].rstrip('\n').partition(' ')
                new = f"{mag_id}__{original}"
                if new in seen:
                    raise ValueError(f"contig name {new} occurs twice (MAGs {seen[new]} and {mag_id})")
                seen[new] = mag_id
                contig_rows.append({'MAG_ID': mag_id, 'New_Name': new_name, 'contig': new,
                                    'original_contig': original, 'original_path': src})
                fout.write(f">{new} {description}\n" if description else f">{new}\n")
            else:
                fout.write(line)
    n_renamed += 1

pd.DataFrame(contig_rows, columns=['MAG_ID', 'New_Name', 'contig', 'original_contig', 'original_path']) \
    .to_csv(os.path.join(renamed_dir, 'contig_map.tsv'), sep='\t', index=False)
print(f"Renamed {n_renamed} MAGs into {renamed_dir}; {len(contig_rows)} contigs renamed "
      f"<MAG_ID>__<original name>, listed in contig_map.tsv")
