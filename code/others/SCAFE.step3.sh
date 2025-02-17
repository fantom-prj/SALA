#!/bin/bash

if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <baseDir> <in_ctss_list> <genome> <ctss_scope_bed_path> <countRegion_bed_path> <SCAFE_code_dir> <outputPrefix>"
    exit 1
fi

# Input parameters
baseDir=$1
in_ctss_list=$2
genome=$3
ctss_scope_bed_path=$4
countRegion_bed_path=$5
scriptDir=$6
tag=$7

#baseDir="./demo_output_local/scafe/count"
#in_ctss_list="step2.CTSS.list.txt"
#genome="hg38.gencode_v39"
#ctss_scope_bed_path="./demo_output_local/scafe/aggregate/run_full/out/annotate/Neuron_series_demo/bed/Neuron_series_demo.cluster.coord.bed.gz"
#countRegion_bed_path="./demo_output_local/scafe/aggregate/run_full/out/annotate/Neuron_series_demo/bed/Neuron_series_demo.CRE.coord.bed.gz" 
#scriptDir="./code/SCAFEv1.0.1/scripts"


outDir="${baseDir}/per_lib/"
matrixDir="${baseDir}/count_matrix"
run_scafe_script_path="${scriptDir}/scafe.tool.bk.count"


# Create directories
mkdir -pm 755 "$outDir"
mkdir -pm 755 "$matrixDir"

# Initialize associative arrays (requires Bash 4 or later)
declare -A lib_info
declare -A data_hsh

# Read the CTSS list file
while IFS=$'\t' read -r libID ctss_bed_path ung_ctss; do
    if [[ ! -s "$ctss_bed_path" ]]; then
        echo "Error: ctss_bed_path $ctss_bed_path does not exist." >&2
        exit 1
    fi
    lib_info["$libID"]="$ctss_bed_path"
done < "$in_ctss_list"

# Process each library
for libID in "${!lib_info[@]}"; do
    ctss_bed_path="${lib_info[$libID]}"
    outputPrefix="$libID"
    count_log_path="${outDir}/${outputPrefix}/log/count.log.tsv"
    lib_info["${libID}_log"]="$count_log_path"

    # Construct the SCAFE command
    cmd="$run_scafe_script_path \
        --countRegion_bed_path=$countRegion_bed_path \
        --ctss_bed_path=$ctss_bed_path \
        --ctss_scope_bed_path=$ctss_scope_bed_path \
        --ctss_scope_slop_bp=0 \
        --genome=$genome \
        --outputPrefix=$outputPrefix \
        --outDir=$outDir"

    # Run the command if the log file doesn't exist
    if [[ ! -s "$count_log_path" ]]; then
        eval "$cmd"
    fi

    # Read the count log
    if [[ -s "$count_log_path" ]]; then
        tail -n +2 "$count_log_path" | while IFS=$'\t' read -r CREID count; do
            data_hsh["$CREID,$libID"]="$count"
        done
    fi
done

# Generate the count matrix
output_file="${matrixDir}/${tag}.count.txt"

# Extract unique library IDs (columns for the matrix)
library_ids=()
for libID in "${!lib_info[@]}"; do
    # Exclude any keys ending with "_log"
    if [[ ! "$libID" =~ _log$ ]]; then
        library_ids+=("$libID")
    fi
done

{
    # Print the header row
    echo -e "CREID\t$(printf "%s\t" "${library_ids[@]}" | sed 's/\t$//')"

    # Extract unique CRE IDs (rows for the matrix)
    all_CREIDs=$(printf "%s\n" "${!data_hsh[@]}" | cut -d',' -f1 | sort -u)

    # Print the data rows
    for CREID in $all_CREIDs; do
        echo -n "$CREID"
        for libID in "${library_ids[@]}"; do
            count="${data_hsh["$CREID,$libID"]}"
            echo -ne "\t${count:-0}"  # Use 0 if count is not available
        done
        echo
    done
} > "$output_file"

echo "Count matrix saved to: $output_file"