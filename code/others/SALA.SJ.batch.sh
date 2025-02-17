#!/bin/bash

# Path to your sample list (Tab-separated: Sample ID, BAM file path)
SAMPLE_LIST="./data/bam.list_b4TC.txt"

# Loop through each line in the sample list
while IFS=$'\t' read -r SAMPLE_ID BAM_PATH; do
    # Define output paths and parameters
    OUT_PREFIX="${SAMPLE_ID}"
    OUT_DIR="./demo_output_local/input/junction_extractor/output"
    
    # Run the junction_extractor Perl script with the specified parameters
    ./code/SALA/junction_extractor/junction_extractor_v0.1.pl \
    --in_bam="$BAM_PATH" \
    --chrom_size_path=./resources/chrom.sizes.tsv \
    --chrom_fasta_path=./resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz \
    --out_prefix="$OUT_PREFIX" \
    --out_dir="$OUT_DIR" \
    --max_thread=1 \
    --min_nt_qual=10 \
    --min_MAPQ=20 \
    --samtools_bin=./resources/bin/samtools/samtools \
    --bedtools_bin=./resources/bin/bedtools/bedtools \
    --tabix_bin=./resources/bin/tabix/tabix \
    --bgzip_bin=./resources/bin/bgzip/bgzip
done < "$SAMPLE_LIST"
