"""
Reads mag_summary.tsv and copies MAG FASTAs into renamed_mags/ using the
new_name column. Edit new_name in the table and re-run this rule to update names.
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

required = {'new_name', 'original_path'}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"mag_summary.tsv is missing required columns: {missing}")

# Clear destination so stale files from previous runs don't accumulate
if os.path.isdir(renamed_dir):
    shutil.rmtree(renamed_dir)
os.makedirs(renamed_dir)

for _, row in df.iterrows():
    src = str(row['original_path'])
    new_name = str(row['new_name'])
    dst = os.path.join(renamed_dir, f"{new_name}.fa")
    if not os.path.exists(src):
        print(f"Warning: source FASTA not found, skipping: {src}")
        continue
    shutil.copy2(src, dst)

print(f"Renamed {len(df)} MAGs into {renamed_dir}")
