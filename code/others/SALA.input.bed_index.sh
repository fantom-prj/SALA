#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <results_path> <tag>"
    exit 1
fi

# Assign input parameters
results_path="$1"
tag="$2"

# Define file paths

file1=$3
file2=$4
file3=$5
file4=$6
file5="${file4}.gz"

tabix_bin="tabix"
bgzip_bin="bgzip"

# Function to compress and index BED files
compress_and_index() {
    local filepath="$1"
    
    if [[ -f "$filepath" ]]; then
        local outpath="${filepath%.gz}.bgz"  # Change .gz to .bgz for output filename
        zcat "$filepath" | "$bgzip_bin" > "$outpath"
        "$tabix_bin" -p bed "$outpath"
        echo "Processed: $outpath"
    else
        echo "Error: File $filepath not found!" >&2
    fi
}

gzip -c "$file4" > "$file5"

compress_and_index "$file1"
compress_and_index "$file2"
compress_and_index "$file3"
compress_and_index "$file5"

echo "All files processed successfully!"
