#!/bin/bash

#Usage: <sala_path> <scafe_path> <resources_path> <out_prefix> <cpat_path[optional]>

# Check if the correct number of arguments is provided
if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <sala_path> <scafe_path> <resources_path> <out_prefix> [cpat_path (optional)]"
    exit 1
fi

# Get command-line parameters (or use default values)
sala_path=${1:-'./demo_output_local/sala/transcript/Neuron_series_demo'}
scafe_path=${2:-'./demo_output_local/scafe'}
resources_path=${3:-'./resources'}
out_prefix=${4:-'Neuron_series_demo'}
cpat_path=${5:-''}  # If empty, CPAT will be skipped

# Fixed paths (not passed as parameters)
genome_fasta="$resources_path/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
hexamer_file="$resources_path/CPAT/Human_Hexamer.tsv"
logit_model="$resources_path/CPAT/Human_logitModel.RData"
cpat_log="$cpat_path/CPAT_run_info.log"

# Step 1: Convert BED12 to BED6
bed12ToBed6 -i "$sala_path/bed/$out_prefix.model.bed.bgz" | gzip > "$sala_path/bed/$out_prefix.model.bed6.bed.gz"

# Step 2: Detect GENCODE 3n end from n3 cluster 
bedtools intersect -s -wa -wb -a "$sala_path/bed/$out_prefix.end3.bed.bgz" -b "$resources_path/GENCODE_info/gencode.3n.bed.gz" | gzip > "$sala_path/bed/$out_prefix.end3.cluster.region.gencode3n.bed.gz"

# Step 3: Detect GENCODE 5n and SCAFE cluster from n5 cluster
bedtools intersect -s -wa -wb -a "$sala_path/bed/$out_prefix.end5.bed.bgz" -b "$scafe_path/aggregate/run_full/out/annotate/$out_prefix/bed/$out_prefix.cluster.coord.bed.gz" | gzip > "$sala_path/bed/$out_prefix.end5.cluster.region.cluster.bed.gz"

bedtools intersect -s -wa -wb -a "$sala_path/bed/$out_prefix.end5.bed.bgz" -b "$resources_path/GENCODE_info/gencode.5n.bed.gz" | gzip > "$sala_path/bed/$out_prefix.end5.cluster.region.gencode5n.bed.gz"

# Step 4: Extract FASTA sequences
bedtools getfasta -s -nameOnly -split -fi "$genome_fasta" -bed "$sala_path/bed/$out_prefix.model.bed.bgz" > "$sala_path/bed/$out_prefix.model.fasta"

# Step 5: Run CPAT only if cpat_path is provided
if [[ -n "$cpat_path" ]]; then
    echo "Running CPAT..."
    mkdir -p "$cpat_path"  # Ensure CPAT output directory exists
    cpat -x "$hexamer_file" -d "$logit_model" --top-orf=5 -g "$sala_path/bed/$out_prefix.model.fasta" -o "$cpat_path/output" --log-file "$cpat_log" >/dev/null 2>&1
else
    echo "Skipping CPAT (No output path provided)"
fi

# Step 6: Annotate promoter-type
bedtools closest -a "$scafe_path/aggregate/run_full/out/annotate/$out_prefix/bed/$out_prefix.CRE.coord.bed.gz"  -b "$resources_path/CRE.from.SCREEN.encodev3/GRCh38-ELS.all.enhancer.sort.bed" -D a> "$scafe_path/aggregate/run_full/out/annotate/$out_prefix/bed/$out_prefix.CRE.coord.e.bed"
bedtools closest -a "$scafe_path/aggregate/run_full/out/annotate/$out_prefix/bed/$out_prefix.CRE.coord.bed.gz"  -b "$resources_path/CRE.from.SCREEN.encodev3/GRCh38-PLS.all.promoter.sort.bed" -D a> "$scafe_path/aggregate/run_full/out/annotate/$out_prefix/bed/$out_prefix.CRE.coord.p.bed"
bedtools closest -a "$scafe_path/aggregate/run_full/out/annotate/$out_prefix/bed/$out_prefix.CRE.coord.bed.gz"  -b "$resources_path/CRE.from.SCREEN.encodev3/GRCh38-CTCF.sort.bed" -D a> "$scafe_path/aggregate/run_full/out/annotate/$out_prefix/bed/$out_prefix.CRE.coord.ctcf.bed"

bedtools closest -a "$scafe_path/aggregate/run_full/out/annotate/$out_prefix/bed/$out_prefix.CRE.coord.bed.gz" -b "$resources_path/F5_enhancer/hg38_robust_enhancers.sort.bed" -D a> "$scafe_path/aggregate/run_full/out/annotate/$out_prefix/bed/$out_prefix.CRE.coord.andersson.robust.e.bed"
bedtools closest -a "$scafe_path/aggregate/run_full/out/annotate/$out_prefix/bed/$out_prefix.CRE.coord.bed.gz" -b "$resources_path/F5_enhancer/hg38_permissive_enhancers.sort.bed" -D a> "$scafe_path/aggregate/run_full/out/annotate/$out_prefix/bed/$out_prefix.CRE.coord.andersson.permissive.e.bed"
