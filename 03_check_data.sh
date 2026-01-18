# Check which GCST folders actually contain usable files  
#!/bin/bash

# Absolute paths
RAW_DIR="/home/dyu/DanYu_BVSc_PhD/yudan_sumstats/raw_sumstats"
OUT_FILE="/home/dyu/DanYu_BVSc_PhD/yudan_sumstats/gcst_with_data.txt"

# Loop through GCST folders
for d in "$RAW_DIR"/GCST*; do
  if ls "$d"/*.{tsv.gz,csv,txt.gz} >/dev/null 2>&1; then
    echo "$(basename "$d")"
  fi
done > "$OUT_FILE"
