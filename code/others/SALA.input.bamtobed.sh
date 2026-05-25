#!/bin/bash

# Fail if any command in a pipe fails (e.g. sort killed by OOM)
set -euo pipefail

# Check if the input file is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <bam_list.txt> <out_dir>"
    exit 1
fi

bam_list="$1"  # The text file containing BAM file paths
output_dir="$2"

output_dir2="$2/bed"
mkdir -p "$output_dir"
mkdir -p "$output_dir2"

all_bed="$output_dir/combined.bed.bgz"
bgzip="bgzip"
tabix="tabix"
bedtools="bedtools"
samtools="samtools"

# Temporary file to store all BED entries before sorting
tmp_all_bed=$(mktemp)

# Read the BAM list line by line
while IFS=$'\t' read -r index_id sample_id bam_path; do
    echo "Processing $sample_id with BAM: $bam_path"

    # Define output paths
    bed_output="$output_dir2/${sample_id}.bed.bgz"

    # Step 1: Convert BAM to BED12 and compress
    "$samtools" view -bh "$bam_path" | \
        "$bedtools" bamtobed -i stdin -bed12 | \
        sort -k1,1 -k2,2n | \
        "$bgzip" -c > "$bed_output"

    # Step 2: Modify column 4 and recompress
    # Add library ID prefix
    # Add sequential suffix for secondary alignments
    tmp_output="$output_dir2/${sample_id}_output_tmp.bed.bgz"

    "$bgzip" -d -c "$bed_output" | \
        awk -v id="$index_id" 'BEGIN{OFS="\t"} {$4=id"_"$4; print}' | \
        awk 'BEGIN{OFS="\t"} {count[$4]++; $4=$4"_"count[$4]; print}' | \
        "$bgzip" -c > "$tmp_output" && mv "$tmp_output" "$bed_output"

    # Step 3: Index the BED file
    "$tabix" -p bed "$bed_output"

    # Append data to temporary file for merging
    bgzip -d -c "$bed_output" >> "$tmp_all_bed"

    echo "Completed processing for $sample_id"
done < "$bam_list"

# Step 4: Sort and compress all combined BED entries
sort --parallel=4 -k1,1 -k2,2n "$tmp_all_bed" | "$bgzip" -c > "$all_bed"

# Step 5: Index the merged all.bed.bgz file
"$tabix" -p bed "$all_bed"

# Remove temporary file
rm "$tmp_all_bed"

echo "All samples processed! Merged BED saved as $all_bed"
