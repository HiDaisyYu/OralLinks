#!/bin/bash

INPUT_FILE="gwas_targets.tsv"
DOWNLOAD_DIR="./raw_sumstats"
LOG_FILE="download.log"
PARALLEL_JOBS=100   # adjust based on bandwidth

mkdir -p "$DOWNLOAD_DIR"

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: $INPUT_FILE not found. Run the python filter script first."
    exit 1
fi

download_one() {
    category="$1"
    trait="$2"
    accession="$3"
    sample_size="$4"

    target_dir="$DOWNLOAD_DIR/$category"
    mkdir -p "$target_dir"

    if ls "$target_dir/${accession}"*.h.tsv.gz >/dev/null 2>&1; then
        echo "[SKIP] $accession already exists."
        return
    fi

    num_id=$(echo "$accession" | sed 's/^GCST//' | sed 's/^0*//')
    if ! [[ "$num_id" =~ ^[0-9]+$ ]]; then
        echo "Warning: Invalid ID format $accession"
        return
    fi

    lower=$(( ( (num_id - 1) / 1000 ) * 1000 + 1 ))
    upper=$(( lower + 999 ))
    range_str=$(printf "GCST%06d-GCST%06d" $lower $upper)
    base_url="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/$range_str/$accession"

    echo "⬇️  [$category] Downloading $accession ($trait)..."

    harm_list=$(wget -q -O - "$base_url/harmonised/" | grep -oE '[0-9A-Za-z._-]+\.h\.tsv\.gz')
    if [ -n "$harm_list" ]; then
        for harm_file in $harm_list; do
            wget -q --show-progress -c "$base_url/harmonised/$harm_file" -O "$target_dir/$harm_file"
            echo "   ✅ Success (Harmonised): $accession ($harm_file)" >> "$LOG_FILE"
        done
    else
        raw_list=$(wget -q -O - "$base_url/" | grep -oE '[0-9A-Za-z._-]+\.tsv\.gz')
        if [ -n "$raw_list" ]; then
            for raw_file in $raw_list; do
                wget -q --show-progress -c "$base_url/$raw_file" -O "$target_dir/$raw_file"
                echo "   ✅ Success (Raw format): $accession ($raw_file)" >> "$LOG_FILE"
            done
        else
            echo "   ❌ Failed: No data found for $accession" >> "$LOG_FILE"
        fi
    fi
}

export -f download_one
export DOWNLOAD_DIR LOG_FILE

# Run in parallel, splitting on tabs
tail -n +2 "$INPUT_FILE" | parallel --colsep '\t' -j $PARALLEL_JOBS download_one {1} {2} {3} {4}

# Clean up zero-byte junk files
find "$DOWNLOAD_DIR" -type f -size 0 -delete

echo "🎉 All jobs finished!"
