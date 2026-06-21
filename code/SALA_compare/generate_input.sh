#!/bin/bash
project_root=$(pwd)

outpath="$project_root/code/SALA_compare/CAT/input"
mkdir -p "$outpath"
cd "$outpath"

base_tag='Neuron_series_demo'
ref_tag='CAT'

transcript_bed_to_end_bed_bigwig="$project_root/code/SALA_compare/transcript_bed_to_end_bed_bigwig.pl"

assemble_bed_gz="$project_root/demo_output_local/sala/transcript/Neuron_series_demo/log/table4.All_Ref.bed12.bed.gz"

assemble_info_tsv_gz="$project_root/demo_output_local/sala/gene/table_filtered_gene/Neuron_series_demo/log/Neuron_series_demo.model.info.tsv.gz"

ref_bed_gz="$project_root/code/SALA_compare/CAT/lv2_permissive.trnscpt.hg38.bed.gz"

base_rng_end5_gz="$project_root/demo_output_local/sala/transcript/Neuron_series_demo/bed/Neuron_series_demo.end5.bed.bgz"

base_rng_end3_gz="$project_root/demo_output_local/sala/transcript/Neuron_series_demo/bed/Neuron_series_demo.end3.bed.bgz"

merge_d_end5=75
merge_d_end3=150
base_rng_end5=./$base_tag.end5.bed
base_rng_end3=./$base_tag.end3.bed

cd $outpath

gzip -dc $base_rng_end5_gz | cut -f 1-6 >$base_rng_end5
gzip -dc $base_rng_end3_gz | cut -f 1-6 >$base_rng_end3

zcat $assemble_info_tsv_gz | tail -n +2 | awk -F'\t' '{print $1"\t"$7"\t"$9"\t"$9"\t"$8}' >./$base_tag.transcript_to_gene.tsv
zcat $assemble_bed_gz | bgzip -c >./$base_tag.bed.bgz
tabix -p bed ./$base_tag.bed.bgz

zcat $ref_bed_gz | bgzip -c >./$ref_tag.bed.bgz
tabix -p bed ./$ref_tag.bed.bgz

zcat ./$ref_tag.bed.bgz ./$base_tag.bed.bgz | sort -k1,1 -k2,2n | bgzip -c >./$ref_tag.$base_tag.transcript.bed.bgz
tabix -p bed ./$ref_tag.$base_tag.transcript.bed.bgz

#$transcript_bed, $merge_d_end5, $merge_d_end3, $base_rng_end5, $base_rng_end3, $outPrefix, $outDir
perl $transcript_bed_to_end_bed_bigwig ./$ref_tag.$base_tag.transcript.bed.bgz $merge_d_end5 $merge_d_end3 $base_rng_end5 $base_rng_end3 $ref_tag.$base_tag $outpath
