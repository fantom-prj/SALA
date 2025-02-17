#!/bin/bash

# Get the directory of the script and move there
cwd=$(dirname "$0")
cd "$cwd"

# Check if the correct number of arguments is provided
if [ "$#" -ne 8 ]; then
    echo "Usage: $0 <transcript_bed> <outDir> <outPrefix> <genome> <chrom_size_path> <sala_path> <scafe_path> <bedGraphToBigWig_bin>"
    exit 1
fi

# Allow first 8 parameters to be passed as input arguments
transcript_bed=${1:-'./demo_output_local/input/bam_to_bed/Neuron_series_demo.bed.bgz'}
outDir=${2:-'./demo_output_local/input/CTES_clusters'}
outPrefix=${3:-'Neuron_series_demo'}
genome=${4:-'hg38.gencode_v39'}
chrom_size_path=${5:-'./resources/hg38.chrom.sizes.sorted.txt'}
sala_path=${6:-'./code/SALA/step2.call.pl'}
scafe_path=${7:-'./code/SCAFEv1.0.1'}
bedGraphToBigWig_bin=${8:-'./resources/bin/bedGraphToBigWig/bedGraphToBigWig'}


CTES_path="$outDir/$outPrefix.end3.bed"
outDir1="$outDir/end3_bed_bigwig"
outDir2="$outDir/scafe/cluster"
outDir3="$outDir/scafe/ctss_to_bigwig"
scafe_cluster="$scafe_path/scripts/scafe.tool.cm.cluster"
scafe_ctss_to_bigwig="$scafe_path/scripts/scafe.tool.cm.ctss_to_bigwig"
transcript_bed_to_end_bed_bigwig="$sala_path/3n_cluster/transcript_bed_to_end_bed_bigwig.pl"

# Step 1: Run transcript_bed_to_end_bed_bigwig
perl "$transcript_bed_to_end_bed_bigwig" "$transcript_bed" "$chrom_size_path" "$outPrefix" "$outDir1" "$bedGraphToBigWig_bin"
cp "$0" "$outDir1"

# Step 2: Run clustering with different parameters
$scafe_cluster \
--overwrite=yes \
--cluster_ctss_bed_path="$CTES_path" \
--count_ctss_bed_path="$CTES_path" \
--min_summit_count=2 \
--min_nt_count=2 \
--min_cluster_count=3 \
--outputPrefix="$outPrefix.CTES.s2_n3_c3" \
--outDir="$outDir2"

$scafe_cluster \
--overwrite=yes \
--cluster_ctss_bed_path="$CTES_path" \
--count_ctss_bed_path="$CTES_path" \
--min_summit_count=3 \
--min_nt_count=3 \
--min_cluster_count=5 \
--outputPrefix="$outPrefix.CTES.s3_n3_c5" \
--outDir="$outDir2"

$scafe_cluster \
--overwrite=yes \
--cluster_ctss_bed_path="$CTES_path" \
--count_ctss_bed_path="$CTES_path" \
--min_summit_count=5 \
--min_nt_count=5 \
--min_cluster_count=10 \
--outputPrefix="$outPrefix.CTES.s5_n5_c10" \
--outDir="$outDir2"

# Step 3: Convert CTSS to BigWig
$scafe_ctss_to_bigwig \
--ctss_bed_path="$CTES_path" \
--genome="$genome" \
--outputPrefix="$outPrefix.CTES" \
--outDir="$outDir3"