#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <results_path> <tag>"
    exit 1
fi

# Assign input parameters
results_path="$1"
tag="$2"

# Define file paths

file1="${results_path}/annotate/${tag}/bed/${tag}.cluster.coord.bed.gz"
file2="${results_path}/aggregate/${tag}/bed/${tag}.aggregate.collapse.ctss.bed.gz"
file3="${results_path}/ctes_cluster/${tag}/bed/${tag}.tssCluster.bed.gz"
file4="${results_path}/ctes_signal/${tag}/${tag}.end3.bed"
file5="${results_path}/ctes_signal/${tag}/${tag}.end5.bed"
file6="${file4}.gz"
file7="${file5}.gz"

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

# Process each file
gzip -c "$file4" > "$file6"
gzip -c "$file5" > "$file7"

compress_and_index "$file1"
compress_and_index "$file2"
compress_and_index "$file3"
compress_and_index "$file6"
compress_and_index "$file7"

echo "All files processed successfully!"
