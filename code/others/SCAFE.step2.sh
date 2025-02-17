#!/bin/bash

# Check if enough arguments are provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <lib_list_path> <genome> <baseDir> <outputPrefix> <SCAFE_code_dir>"
    exit 1
fi

# Input parameters
lib_list_path=$1
genome=$2
baseDir=$3
outputPrefix=$4
scriptDir=$5

#lib_list_path='step2.CTSS.list.txt'
#genome="hg38.gencode_v39"
#baseDir='./demo_output_local/scafe/aggregate/run_full'
#outputPrefix="Neuron_series_demo"
#scriptDir="./code/SCAFEv1.0.1/scripts"

$scriptDir/scafe.tool.cm.aggregate \
--lib_list_path=$lib_list_path \
--max_thread=5 \
--genome=$genome \
--outputPrefix=$outputPrefix \
--outDir=$baseDir/out/aggregate

$scriptDir/scafe.tool.cm.cluster \
--overwrite=yes \
--cluster_ctss_bed_path=$baseDir/out/aggregate/$outputPrefix/bed/$outputPrefix.aggregate.collapse.ctss.bed.gz \
--count_ctss_bed_path=$baseDir/out/aggregate/$outputPrefix/bed/$outputPrefix.aggregate.unencoded_G.collapse.ctss.bed.gz \
--min_summit_count=0 \
--min_cluster_count=1 \
--outputPrefix=$outputPrefix \
--outDir=$baseDir/out/cluster

$scriptDir/scafe.tool.cm.ctss_to_bigwig \
--genome=$genome \
--ctss_bed_path=$baseDir/out/aggregate/$outputPrefix/bed/$outputPrefix.aggregate.collapse.ctss.bed.gz \
--outputPrefix=$outputPrefix.all \
--outDir=$baseDir/out/ctss_to_bigwig

$scriptDir/scafe.tool.cm.ctss_to_bigwig \
--genome=$genome \
--ctss_bed_path=$baseDir/out/aggregate/$outputPrefix/bed/$outputPrefix.aggregate.unencoded_G.collapse.ctss.bed.gz \
--outputPrefix=$outputPrefix.ung \
--outDir=$baseDir/out/ctss_to_bigwig

$scriptDir/scafe.tool.cm.filter \
--overwrite=yes \
--ctss_bed_path=$baseDir/out/aggregate/$outputPrefix/bed/$outputPrefix.aggregate.collapse.ctss.bed.gz \
--ung_ctss_bed_path=$baseDir/out/aggregate/$outputPrefix/bed/$outputPrefix.aggregate.unencoded_G.collapse.ctss.bed.gz \
--tssCluster_bed_path=$baseDir/out/cluster/$outputPrefix/bed/$outputPrefix.tssCluster.bed.gz \
--genome=$genome \
--outputPrefix=$outputPrefix \
--outDir=$baseDir/out/filter \

$scriptDir/scafe.tool.cm.annotate \
--overwrite=yes \
--tssCluster_bed_path=$baseDir/out/filter/$outputPrefix/bed/$outputPrefix.tssCluster.default.filtered.bed.gz \
--tssCluster_info_path=$baseDir/out/filter/$outputPrefix/log/$outputPrefix.tssCluster.log.tsv \
--min_CRE_count=3 \
--genome=$genome \
--outputPrefix=$outputPrefix \
--outDir=$baseDir/out/annotate

$scriptDir/scafe.tool.cm.directionality \
--overwrite=yes \
--CRE_bed_path=$baseDir/out/annotate/$outputPrefix/bed/$outputPrefix.CRE.coord.bed.gz \
--CRE_info_path=$baseDir/out/annotate/$outputPrefix/log/$outputPrefix.CRE.info.tsv.gz \
--ctss_bed_path=$baseDir/out/aggregate/$outputPrefix/bed/$outputPrefix.aggregate.collapse.ctss.bed.gz \
--outputPrefix=$outputPrefix \
--outDir=$baseDir/out/directionality/
