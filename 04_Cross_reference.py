# Cross‑reference with your trait index file  
import pandas as pd

# Load trait index
df = pd.read_csv("/home/dyu/DanYu_BVSc_PhD/yudan_sumstats/gwas_targets.index.tsv", sep="\t")

# Load valid GCST IDs
with open("/home/dyu/DanYu_BVSc_PhD/yudan_sumstats/gcst_with_data.txt") as f:
    valid_ids = set(line.strip() for line in f)

# Keep only rows with valid GCST IDs
df_filtered = df[df["GCST_ID"].isin(valid_ids)]

# Save filtered file
df_filtered.to_csv("/home/dyu/DanYu_BVSc_PhD/yudan_sumstats/gwas_targets.filtered.tsv", sep="\t", index=False)

print(f"Original index: {len(df)} rows")
print(f"Valid GCST IDs: {len(valid_ids)}")
print(f"Filtered index: {len(df_filtered)} rows")
