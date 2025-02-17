#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <scafe_path> <n3_cluster_path> <tag> <resource_directory>"
    exit 1
fi

# Assign input parameters
scafe_path="$1"
n3_cluster_path="$2"
tag="$3"
resource_directory="$4"

# Define file paths
file1="${scafe_path}/aggregate/run_full/out/annotate/${tag}/bed/${tag}.cluster.coord.bed.gz"
file2="${scafe_path}/aggregate/run_full/out/aggregate/${tag}/bed/${tag}.aggregate.collapse.ctss.bed.gz"
file3="${n3_cluster_path}/scafe/cluster/${tag}.CTES.s3_n3_c5/bed/${tag}.CTES.s3_n3_c5.tssCluster.bed.gz"
file4="${n3_cluster_path}/end3_bed_bigwig/${tag}.end3.bed"
file5="${n3_cluster_path}/end3_bed_bigwig/${tag}.end5.bed"
file6="${n3_cluster_path}/end3_bed_bigwig/${tag}.end3.bed.gz"
file7="${n3_cluster_path}/end3_bed_bigwig/${tag}.end5.bed.gz"

tabix_bin="${resource_directory}/bin/tabix/tabix"
bgzip_bin="${resource_directory}/bin/bgzip/bgzip"

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
