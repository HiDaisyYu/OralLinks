#!/bin/bash

INPUT_FILE="gwas_targets.tsv"
CLEAN_FILE="gwas_targets.cleaned.tsv"
DOWNLOAD_DIR="./raw_sumstats"
LOG_FILE="download.log"
SUCCESS_LOG="success.log"
PARALLEL_JOBS=50   # tune based on bandwidth

mkdir -p "$DOWNLOAD_DIR"
: > "$LOG_FILE"
: > "$SUCCESS_LOG"

# --- Clean the input file ---
# 1. Keep only rows with valid GCST accession in column 3
# 2. Remove quotes from Trait
# 3. Replace commas in Sample_Size with spaces
# 4. Collapse multiple spaces
awk -F'\t' 'NR==1{print;next}{
  if ($3 ~ /^GCST[0-9]+$/) {
    gsub(/"/,"",$2);          # remove quotes
    gsub(/,/," ",$4);         # replace commas with spaces
    gsub(/  +/," ",$4);       # collapse multiple spaces
    print $1"\t"$2"\t"$3"\t"$4
  }
}' "$INPUT_FILE" > "$CLEAN_FILE"

download_one() {
    category="$1"
    trait="$2"
    accession="$3"
    sample_size="$4"

    target_dir="$DOWNLOAD_DIR/$accession"
    mkdir -p "$target_dir"

    # Skip if already downloaded
    if ls "$target_dir/${accession}"*.tsv.gz >/dev/null 2>&1 || \
       ls "$target_dir/${accession}"*.csv >/dev/null 2>&1; then
        echo "[SKIP] $accession already exists."
        return
    fi

    # Validate accession
    if [[ ! "$accession" =~ ^GCST[0-9]+$ ]]; then
        echo "Warning: Invalid ID format $accession"
        return
    fi

    num_id=${accession#GCST}
    num_id=$(echo "$num_id" | sed 's/^0*//')

    lower=$(( ( (num_id - 1) / 1000 ) * 1000 + 1 ))
    upper=$(( lower + 999 ))
    range_str=$(printf "GCST%08d-GCST%08d" $lower $upper)
    base_url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/$range_str/$accession"

    echo "⬇️  [$category] Downloading $accession ($trait)..."

    # Try harmonised first
    harm_list=$(wget -q -O - "$base_url/harmonised/" | grep -oE '[0-9A-Za-z._-]+\.h\.tsv\.gz')
    if [ -n "$harm_list" ]; then
        for harm_file in $harm_list; do
            wget -q --show-progress -c "$base_url/harmonised/$harm_file" -O "$target_dir/$harm_file"
            echo "   ✅ Success (Harmonised): $accession ($harm_file)" >> "$LOG_FILE"
            echo -e "$accession\t$trait\t$category\t$harm_file" >> "$SUCCESS_LOG"
        done
        return
    fi

    # Then try raw .tsv.gz
    raw_list=$(wget -q -O - "$base_url/" | grep -oE '[0-9A-Za-z._-]+\.tsv\.gz')
    if [ -n "$raw_list" ]; then
        for raw_file in $raw_list; do
            wget -q --show-progress -c "$base_url/$raw_file" -O "$target_dir/$raw_file"
            echo "   ✅ Success (Raw): $accession ($raw_file)" >> "$LOG_FILE"
            echo -e "$accession\t$trait\t$category\t$raw_file" >> "$SUCCESS_LOG"
        done
        return
    fi

    # Finally try raw .csv
    csv_list=$(wget -q -O - "$base_url/" | grep -oE '[0-9A-Za-z._-]+\.csv')
    if [ -n "$csv_list" ]; then
        for csv_file in $csv_list; do
            wget -q --show-progress -c "$base_url/$csv_file" -O "$target_dir/$csv_file"
            echo "   ✅ Success (CSV): $accession ($csv_file)" >> "$LOG_FILE"
            echo -e "$accession\t$trait\t$category\t$csv_file" >> "$SUCCESS_LOG"
        done
        return
    fi

    echo "   ❌ Failed: No data found for $accession" >> "$LOG_FILE"
}

export -f download_one
export DOWNLOAD_DIR LOG_FILE SUCCESS_LOG

# Run in parallel on the cleaned file
tail -n +2 "$CLEAN_FILE" | parallel --colsep '\t' -j $PARALLEL_JOBS download_one {1} {2} {3} {4}

# Clean up zero-byte junk files
find "$DOWNLOAD_DIR" -type f -size 0 -delete

# make a index file: system, disease, GCST ID. 
awk -F'\t' 'NR==1{print "System\tTrait\tGCST_ID";next}{print $1"\t"$2"\t"$3}' gwas_targets.clean.tsv > gwas_targets.index.tsv
echo "🎉 All jobs finished!"
