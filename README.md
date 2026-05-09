# SALA
<div style="text-align:center"><img src="figs/SALA.png" width="500"></div>

The transcript Start-site Aware Long-read Assembler (SALA) is developed for de novo assembling long-read into transcript and gene models, considering support from confident transcription start site. SALA incorporates confident TSS clusters de novo identified from the long-read data or pre-defined confident TSS clusters.

## Table of contents
* [Dependencies](#depend)
* [Installation](#installation)
* [Running SALA](#how_to_run)
  * [Assembling into transcript models](#transcript_model)
  * [Grouping into initial gene models](#gene_group)
  * [Count matrix per transcript model](#transcript_count)
  * [Annotating the transcript model](#transcript_filter)
  * [Getting final gene model](#final_gene_model)
  * [Final log table and gtf output](#final_table_gtf)
* [Working with the SALA results](#SALA_result)
* [Citing SALA](#SALA_cite)
* [Contribution](#contribution)
* [References](#ref)


# <a name="depend"></a>Dependencies

## SALA requires the following tools to run
* Perl (tested with v5.26.2, installed from https://www.perl.org/get.html)
* R (tested with v4.2, installed from: https://cran.r-project.org/)
* SCAFE (Moody et al. 2022) (tested with v1.0.1, located at ./code/SCAFEv1.0.1/scripts)
* paraclu (Frith et al. 2008) (located at ./resources/bin/paraclu/paraclu)
* samtools (Danecek et al. 2021) (tested with v1.11 , located at ./resources/bin/samtools/samtools)
* bedtools (Quinlan and Hall 2010) (tested with v2.30.0 , located at ./resources/bin/bedtools/bedtools)
* tabix (Li 2011) (tested with v1.15.1 , located at ./resources/bin/tabix/tabix)
* bgzip (Li 2011) (tested with v1.15.1 , located at ./resources/bin/bgzip/bgzip)
* bedGraphToBigWig (Yao et al. 2017) (tested with version 2.8, located at ./resources/bin/bedGraphToBigWig/bedGraphToBigWig)

## The followings tools are recommended to install
* TranscriptClean (Wyman and Mortazavi 2019) (tested with v2.0.3, installed from https://github.com/mortazavilab/TranscriptClean)
* bambu (Chen et al. 2023) (tested with v3.2.4, installed from https://github.com/GoekeLab/bambu)
* bedparse (Leonardi 2019) (tested with v0.2.3, installed from https://github.com/tleonardi/bedparse). 

# <a name="installation"></a>Installation
To obtain SALA:
```
#--- make a directory to install SALA
mkdir -pm 755 /my/path/to/install/
cd /my/path/to/install/

#--- Obtain SALA from github
git clone https://github.com/fantom-prj/SALA
cd SALA

#--- export SALA scripts dir to PATH for system-wide call of SALA commands 
echo "export PATH=\$PATH:$(pwd)/code/SCAFEv1.0.1/scripts" >>~/.bashrc
echo "export PATH=\$PATH:$(pwd)/code/SALA" >>~/.bashrc
echo "export PATH=\$PATH:$(pwd)/code/others" >>~/.bashrc
source ~/.bashrc

#--- making sure the scripts and binaries are executable
chmod 755 -R ./code/
chmod 755 -R ./resources/bin/
```

This package itself does not require installation. Essential binary files for Linux platform are included in ./resources/bin (for SALA) and ./code/SCAFEv1.0.1/resources/bin (for SCAFE). If other platform is used, the binary files need to be replaced by the ones from your system. Alternative bin set for Mac OS can be downloaded [here](https://drive.google.com/drive/folders/13SqjSH0eGSH5xnKH4RE43yXtIinvCZet?usp=drive_link). Please replace the downloaded bin folder with the bin folder for SALA and SCAFE.

# <a name="how_to_run"></a>How to run
Please refer to the [wiki page](https://github.com/fantom-prj/SALA/wiki) to run a demo

## <a name="transcript_model"></a>Assembling into transcript models
This tool assigns long-read sequencing data (as query) to a set of reference transcripts (e.g. GENCODE) using a 5' end centric approach. Several reference transcript annotation sets are available: [GENCODE_V39](https://drive.google.com/file/d/14Wd4vab7_8LZi3Z_qyXYbYnTBEoq6Rau/view?usp=drive_link), [GENCODE_V47](https://drive.google.com/file/d/1n2tl-ejmeQwIwgBRDpflyuSNukrAWTtw/view?usp=drive_link), [GENCODE_VM25](https://drive.google.com/file/d/1kI2pL8ZwIQ7b-4D-8O5ZlNUUKYzMAP2J/view?usp=drive_link), [GENCODE_VM36](https://drive.google.com/file/d/1vURkydNDW27Qlxo6hSC-lMAVs9KZjgWJ/view?usp=drive_link). This code will take a set of user-defined confident 5' end clusters (or de novo defined by SCAFE clustering) and 3' end clusters (or de novo defined by clustering) and assign the query reads to the reference transcripts with the following step:

1. A query read is classified as complete if both of its 5' and 3' end overlap a confident 5' and 3' end cluster, otherwise as incomplete.
2. An incomplete read without a confident 5' cluster will be flagged.
3. A complete query read will be assigned to a reference transcript if it shares the same 1) 5' end cluster, 2) 3' end cluster and 3) internal splicing structures (i.e. same splicing junctions or both unspliced).
4. An incomplete query read will be assigned to a reference transcript if it shares the same 5'end cluster and a partial internal splicing structure (i.e. contains part if the reference transcript splicing junctions or unspliced but overlap with reference transcript 1st exon).
5. All unassigned query reads will be flagged as novel.
6. Novel complete query reads with the same 1) 5'end cluster, 2) 3'end cluster and 3) internal splicing structure (i.e. same splicing junctions or both unspliced) will be collapsed as a novel transcript model will new ID assigned.
7. Novel incomplete query reads will be assigned to the novel transcript models if it shares the same 5'end cluster and a partial internal splicing structure (i.e. contains part of the reference transcript splicing junctions or unspliced but overlap with novel transcript models 1st exon).
8. All remaining unassigned novel incomplete query reads will be grouped by their 1) 5'end clusters and 2) internal splicing structure (i.e. same splicing junctions or unspliced) and each group will be collapsed as a novel transcript model with a new ID.
9. The 5' end all transcript models will be adjusted to the summit of the 5'end clusters (de novo or user defined).
10. The 3' end complete transcript models will be adjusted to the summit of the 3'end clusters (de novo or user defined).
11. The 3' end incomplete transcript models will be adjusted to the furthest 3'end of its query reads.

```
Usage: end5_guided_assembler_v1.1.pl [options] --qry_bed_bgz --ref_bed_bgz --out_dir
   
   --qry_bed_bgz                <required> [path]    bed 12 of the long-reads, 4th column must be read ID and in bgz format, 
                                                     for multiple query bed, user can supply a list of path in plain text format, one line one path
   --ref_bed_bgz                <required> [path]    bed 12 of the reference transcript models, 4th column must be transcript ID and in bgz format
   --out_dir                    <required> [path]    output directory
   --chrom_size_path            <required> [path]    a txt file contains the chromsome size in format of chrom\tsize
   --chrom_fasta_path           <required> [path]    genome fasta file
   --conf_end5_bed_bgz          <required> [path]    a bed bgz 12 file contains the 5'end clusters, summit must be provide in the thick end column
   --conf_end3_bed_bgz          <required> [path]    a bed bgz 12 file contains the 3'end clusters, summit must be provide in the thick end column
   --signal_end5_bed_bgz        <required> [path]    a single nucleotide piled up end5 signal bed (ctss bed file) used to define conf_end5_bed_bgz
   --signal_end3_bed_bgz        <required> [path]    a single nucleotide piled up end3 signal bed (ctes bed file) used to define conf_end3_bed_bgz
   --out_prefix                 (optional) [string]  output files prefix, if not defined, qry_bed_bgz filename will be used
   --novel_model_prefix         (optional) [string]  prefix of the novel transcript models [default=ONTC]
   --min_qry_score              (optional) [integer] the minimum score in the query bed file (assumes MAPQ) to be taken for assembly [default=10]
   --conf_end3_merge_flank      (optional) [integer] the flanking distance (on each side) of the 3'end clusters used to merge as a end3 region.
                                                     Use '-1' to turn off. [default=50]
   --conf_end5_merge_flank      (optional) [integer] the flanking distance (on each side) of the 5'end clusters used to merge as a end5 region.
                                                     Use '-1' to turn off. [default=50]
   --conf_end3_add_ref          (optional) [yes/no]  to add reference 3'end into the user defined confident 3'end clusters or not. if yes, the ref 3'end 
                                                     will bed extended by conf_end3_merge_flank nt and merged with confident 3'end clusters
   --min_exon_length            (optional) [integer] minimum length of an exon in a transcript to be considered as valid. If a transcript contains
                                                     an exon shorter than min_exon_length, the transcript will be discarded [default=1]
   --min_transcript_length      (optional) [integer] minimum length of a transcript (including intron) to be considered as valid. If a transcript 
                                                     is shorter than min_transcript_length, the transcript will be discarded [default=50]
   --filter_conf_end5           (optional) [yes/no]  to filter out query reads that is out the original ranges in the conf_end5_bed_bgz
                                                     will bed extended by conf_end3_merge_flank nt and merged with confident 3'end clusters
   --trnscpt_set_end_priority   (optional) [string]  Priority of methods to determine the ends of transcript set? 
                                                     1) based on "summit" : the signal summit in confident end3/end5 clusters, in signal_end*_bed_bgz
                                                     2) based on "commonest" : the observed position that is the most frequent in transcripts of the set
                                                     3) based on "longest": the observed position that is the furtherest in transcripts of the set
                                                     for (1) and (2), there is chances of causing conflicts in the transcript set ranges (e.g. 3'end is
                                                     more downstream than the 5'end in the transcript set). (3) is guranteed to be conflict free.
                                                     use a colon (:) delimited string to indiciate priority e.g. "summit:commonest:longest"
                                                     [default=summit:commonest:longest]
   --enforce_qry_original_end   (optional) [yes/no]  Overrides 
   --print_trnscrptID           (optional) [yes/no]  Print out the transcript ID or not 
   --doubtful_end_merge_dist    (optional) [integer] Distance to merge incomplete ends as groups [default=100]
   --doubtful_end_avoid_summit  (optional) [yes/no]  Overrides --trnscpt_set_end_priority from using "summit" 
   --retain_no_qry_ref_bound_set(optional) [yes/no]  report the bound set or not if the bound set is not detected from the query reads
   --min_summit_dist_split      (optional) [integer] When splitting an end cluster into two, the minimum distance between two summits [default=50]
   --min_size_split             (optional) [integer] When splitting an end cluster into two, the minimum size of the cluster [default=100]
   --min_frac_split             (optional) [integer] When splitting an end cluster into two, the minimum fraction of signal from the two summits [default=0.2]
   --max_thread                 (optional) [integer] number of threads to be used [default=5]
   --bedtools_bin               (optional) [path]    path to the binary of bedtools, if not provided, "bedtools" will be called
   --tabix_bin                  (optional) [path]    path to the binary of tabix, if not provided, "tabix" will be called
   --bgzip_bin                  (optional) [path]    path to the binary of bgzip, if not provided, "bgzip" will be called
   --conf_junction_bed          (optional) [path]    path to confident splicing junction bed

```
Test run:
```
cd /your/SALA/directory/resources
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
gunzip GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz

cd /your/SALA/directory/

./code/SALA/end5_guided_assembler.pl \
--qry_bed_bgz=./test_input/input/bam_to_bed/combined.bed.bgz \
--ref_bed_bgz=./resources/GENCODE_V39/transcript.bed.bgz \
--out_dir=./test_input/sala/transcript \
--chrom_size_path=./resources/chrom.sizes_major.tsv \
--chrom_fasta_path=./resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
--conf_end5_bed_bgz=./test_input/scafe/aggregate/run_full/out/annotate/Neuron_series_demo/bed/Neuron_series_demo.cluster.coord.bed.bgz \
--conf_end3_bed_bgz=./test_input/input/CTES_clusters/scafe/cluster/Neuron_series_demo.CTES.s3_n3_c5/bed/Neuron_series_demo.CTES.s3_n3_c5.tssCluster.bed.bgz \
--signal_end5_bed_bgz=./test_input/scafe/aggregate/run_full/out/aggregate/Neuron_series_demo/bed/Neuron_series_demo.aggregate.collapse.ctss.bed.bgz \
--signal_end3_bed_bgz=./test_input/input/CTES_clusters/scafe/cluster/Neuron_series_demo.CTES.s3_n3_c5/bed/Neuron_series_demo.CTES.s3_n3_c5.tssCluster.bed.bgz \
--out_prefix=Neuron_series_demo \
--novel_model_prefix=ONTT \
--min_output_qry_count=1 \
--conf_end3_merge_flank=150 \
--conf_end5_merge_flank=75 \
--conf_end3_add_ref=yes \
--conf_end5_add_ref=yes \
--min_exon_length=1 \
--min_transcript_length=15 \
--doubtful_end_avoid_summit=yes \
--print_trnscrptID=no \
--trnscpt_set_end_priority=summit:commonest:longest \
--doubtful_end_merge_dist=150 \
--doubtful_end_avoid_summit=yes \
--min_summit_dist_split=50 \
--retain_no_qry_ref_bound_set=no \
--min_size_split=100 \
--min_frac_split=0.2 \
--max_thread=1 \
--conf_junction_bed=./resources/GENCODE_V39/junction.bed,./test_input/input/junction_extractor/pool/Neuron_series_demo.hi_qual.junct.bed \
--min_qry_score=0 \
--bedtools_bin=./resources/bin/bedtools/bedtools \
--tabix_bin=./resources/bin/tabix/tabix \
--bgzip_bin=./resources/bin/bgzip/bgzip

```

## <a name="gene_group"></a>Grouping into initial gene models
This tool annotates the 5'end-guided assembler models as genes. The following script uses the parameter --disable_ref_chain_bound_gene_anno=yes, which prevents recursive gene region extension when unfiltered transcriptional noise is present in the permissive transcript model output. At this stage, the tool assigns preliminary Gene IDs to all transcript models, which will be used in the filtering step.

```
Usage: assemble_gene_annotator_v0.1.pl [options] --qry_bed_bgz --ref_bed_bgz --out_dir
   
    --model_bed_bgz             <required> [path]    generated by end5_guided_assembler, contain info for all reference transcripts
    --model_info_gz             <required> [path]    generated by end5_guided_assembler, contain info for all reference transcripts
    --ref_model_gene_link       <required> [path]    transcript ID to gene ID link of reference transcriptome, columns 1-5: transcript_id, gene_id, transcript_type, gene_type, gene_name
    --out_dir                   <required> [path]    output directory
    --revert_ref_model_bed_bgz  (optional) [path]    if input, will revert the model without qry support back to the bounds in the revert_ref_model_bed_bgz bed input
    --out_prefix                (optional) [string]  output files prefix, if not defined, qry_bed_bgz filename will be used
    --novel_gene_prefix         (optional) [string]  prefix of the novel gene models [default=ONTC]
    --max_thread                (optional) [integer] number of threads to be used [default=5]
    --bedtools_bin              (optional) [path]    path to the binary of bedtools, if not provided, "bedtools" will be called
    --tabix_bin                 (optional) [path]    path to the binary of tabix, if not provided, "tabix" will be called
    --bgzip_bin                 (optional) [path]    path to the binary of bgzip, if not provided, "bgzip" will be called
```
Test run
```
cd /your/SALA/directory/

perl ./code/SALA/assemble_gene_annotator.pl \
--chrom_size_path=./resources/chrom.sizes.tsv \
--model_bed_bgz=./test_input/sala/transcript/Neuron_series_demo/bed/Neuron_series_demo.model.bed.bgz \
--model_info_gz=./test_input/sala/transcript/Neuron_series_demo/log/Neuron_series_demo.model.info.tsv.gz \
--revert_ref_model_bed_bgz=./resources/GENCODE_V39/transcript.bed.bgz \
--ref_model_gene_link=./resources/GENCODE_V39/transcript_to_gene.tsv \
--novel_gene_prefix=IN1G \
--disable_ref_chain_bound_gene_anno=yes \
--min_ref_exon_overlap_pct=10 \
--exon_overlap_dist=-1 \
--locus_merge_dist=100000 \
--exclude_t_type=retained_intron \
--out_prefix=Neuron_series_demo \
--out_dir=./test_input/sala/gene/table0_gene \
--max_thread=1 \
--bedtools_bin=./resources/bin/bedtools/bedtools \
--tabix_bin=./resources/bin/tabix/tabix \
--bgzip_bin=./resources/bin/bgzip/bgzip
```
## <a name="transcript_count"></a>Count matrix per transcript model
This tool parse the information from end5_guided_assembler to generate a count matrix per library across all transcript models.

```
Usage: Rscript SALA.count_matrix.R <SALA_directory> <output_directory> <ref.transcriptome_path>
```

## <a name="transcript_filter"></a>Annotating the transcript model
This step intersects reference transcript models’ 3' ends with SALA 3' end clusters, reference 5' ends with SALA 5' end clusters, and SCAFE-defined confident 5' end clusters with SALA 5' end clusters. It also annotates the promoter types of SCAFE tCREs by intersecting them with GENCODE SCREEN cCREs (The ENCODE Project Consortium et al. 2020) and converts SALA transcripts into an exon BED file. Additionally, an optional coding potential analysis using CPAT (Wang et al. 2013) will be performed if a directory for storing CPAT results is provided. CPAT must be installed separately before running this script and can be obtained from: https://github.com/liguowang/cpat. The script gathers all collected information into both a raw log table and a filtered log table. This step requires R and several R packages, including "Biostrings", "GenomicRanges", "Rsamtools", "data.table", "dplyr", "magrittr", and "tidyr". If these packages are not available, the script will attempt to install them automatically. If installation fails, manual installation will be required.

Overall, this step will do the followings:
1. Intersect 5’ ends and 3’ ends from reference transcript model to SALA 5’ end and 3’ end clusters
2. Intersect ENCODE defined promoter-type to SCAFE confident tCREs
3. Internal priming prediction, and filtering of internal primed novel transcript models
4. Perform PAS motif serach from the transcription end site of all the transcript models.
5. Add initial gene ID  
6. Incorporate full read count per transcript model
7. Add transcript length and exon number
8. Annotate 3’ end clusters, including internal priming prediction, PAS motif search and ref 3' end intersect
9. Annotate 5’ end clusters, including confident tCREs, promoter-typing from ENCODE SCREEN cCRE
10. Define promoter-type per transcript model
11. If CPAT path is provided, CPAT coding potential prediction is performed and incorporate into the output table
12. Export raw table with above details
13. Filter to final table according to the parameters input, and export

```
Usage: Rscript SALA.filter.R <SALA_directory> <out_prefix> <resource_directory> <ref_directory> <fasta_file> <read.per.rep_ref.novel.Tx> <read.per.rep_non-ref.novel.Tx> <isoform_ratio> <require.5'.confidence> <SALA_gene_path> <sample_file> <SCAFE_directory> <genonme> <CPAT_path(optional)>

	SALA_directory                <required>	path of the folder of SALA transcript annotation output
	out_prefix                    <required>	output files prefix
	resource_directory            <required>	path of the resources folder of SALA
	ref_directory                 <required>	path of the folder containing the infomation of reference transcriptome
	fasta_file                    <required>	genome fasta file
	read.per.rep_ref.novel.Tx     <required>	number of reads per replicate from a sample for novel isoform of known gene
	read.per.rep_non-ref.novel.Tx <required>	number of reads per replicate from a sample for novel transcript of novel gene
	isoform_ratio                 <required>	the ratio of a novel isoform across all the transcript in a gene, according to the full length read count
	require.5'.confidence         <required>	if Yes, all the novel transcripts are required to have their 5’ ends located inside confident TSS clusters
	SALA_gene_path                <required>	path of the folder of SALA gene annotation output
	sample_file                   <required>	txt file for the input library: column1, library prefix; column2, sample ID; column3, sampleID with replicate ID
	SCAFE_directory               <required>	path of the folder of SCAFE output (NA if SCAFE is skipped)
	genome                        <required>    genome name: hg38, mm10, mm39 or NA, control promoter typing
	keep_internal_prime           <required>    Yes or No (designed by poly-dT internal priming potential, set Yes if poly-dT priming is not used)
	CPAT_path                     <optional>	path of the folder expected for CPAT result
```

## <a name="final_gene_model"></a>Getting final gene model
In the final gene annotation step, after applying custom filtering, gene annotation should be performed with “--disable_ref_chain_bound_gene_anno=no” to allow proper gene boundary refinement.


## <a name="final_table_gtf"></a>Final log table and gtf output
This script updates the filtered transcript models with the finalized gene annotation from the previous section. It extracts the transcript class and gene class from the reference transcriptome and annotates novel transcript and gene classes if CPAT coding potential results were provided previously. Finally, it generates two GTF files from the filtered transcript models: (1) table4.All_Ref.gtf.gz, which includes filtered novel transcripts and genes along with all transcripts and genes from the reference transcriptome, and (2) table4.Detected_Ref.gtf.gz, which includes filtered novel transcripts and genes along with only the detected transcripts and genes from the reference transcriptome.

```
Usage: Rscript SALA.gene_gtf_annotation.R <SALA_directory> <out_prefix> <resource_directory> <ref_directory> <SALA_gene_path>

	SALA_directory         <required>	path of the folder of SALA transcript annotation output
	out_prefix             <required>	output files prefix
	resource_directory     <required>	path of the resources folder of SALA
	ref_directory          <required>	path of the folder containing the infomation of reference transcriptome
	SALA_gene_path         <required>	path of the folder of SALA final gene annotation output
```

# <a name="SALA_result"></a>Working with the SALA results
Major output from SALA includes the log tables for each transcript model and their annotation, and GTF files including the final novel transcripts. 

Column description of the log table:
```
column          final   if      description
                table   CPAT
                only    run
model_ID                        transcript model ID, use ENST if match with ENST
loc                             coordinate of the transcript region
completeness                    Y if both end5_conf & 3end_conf are confident
full_set_ID                     set ID of this model (intermediate info during assembling)
end5_conf                       end5 clusters derived from SCAFE or GENCODE are considered confident
end3_conf                       end3 clusters derived from paraclu (min. read per cluster > 10) or GENCODE are considered confident
junct_conf                      model contain any junctions without support is considered non-confident (junction support: any GENCODE junction, any defined short-read junction, any canonical junction, junction having all +/- 3 basepair basecalling score >10 (use max if more than one read))
end5_endtype                    position inside a 5end cluster is used as the model 5' end: S: summit; C: commonest; L: most upstream
end3_endtype                    position inside a 3end cluster is used as the model 3' end: S: summit; C: commonest; L: most downstream
novelty                         transcript model is new or from reference transcriptome
strand                          strand of the transcript model
partial_set_count               number of partial set that are collapsed into this model
partial_set_ID_str              all the partial sets that are assigned to this model
full_ref_count                  number of transcript model from GENCODE that are completely match with this model
full_qry_count                  number of ONT-CAGE reads that are completely match with this model
partial_ref_count               number of transcript model from GENCODE that are partially match with this model
partial_qry_count               number of ONT-CAGE reads that are partially match with this model
full_set_bound_str              set string of this model (intermediate info during assembling)
IN1_gene_ID                     intermediate Gene ID for grouping of transcript model; same gene if match for least one of feature: 5end cluster, 3end cluster, any junction
IN1_gene_name                   intermediate gene name, use gene name of reference transcriptome
gene_novelty                    gene model is new or from reference transcriptome: Novel / Ref
transcript_novelty              transcript model is new or from reference transcriptome: Novel / Ref
iPSC_rep1                       number of reads with full-length match with this model
iPSC_rep2                       number of reads with full-length match with this model
NSC_rep1                        number of reads with full-length match with this model
NSC_rep2                        number of reads with full-length match with this model
Neuron_rep1                     number of reads with full-length match with this model
Neuron_rep2                     number of reads with full-length match with this model
iPSC                            number of reads with full-length match with this model; sum of the provided replicates
NSC                             number of reads with full-length match with this model; sum of the provided replicates
Neuron                          number of reads with full-length match with this model; sum of the provided replicates
isoform_filter                  permissive or standard, according to the read count filter per replicate for novel isoform of reference genes
novel_gene_Tx_filter            permissive or standard, according to the read count filter per replicate for novel transcript of novel genes
ref_source                      "non_detectable_ref", "fulllength_ref", "partial_ref" or "novel_transcript"
Tx_ratio_iPSC                   transcript ratio of gene (IN1_gene), only consider full-length count
Tx_ratio_NSC                    transcript ratio of gene (IN1_gene), only consider full-length count
Tx_ratio_Neuron                 transcript ratio of gene (IN1_gene), only consider full-length count
max_T_ratio                     maximum transcript ratio of all the samples
n_exon                          number of exon 
transcript_length               transcript length (exclude intron)
n3_string                       ID of end3 cluster
internal_priming                Yes if transcript model 3' end hit internal priming site 
n3_Reference                    Yes if end3 cluster intersect with any Reference 3' end
n3_support                      a collection of support by Reference 3' end, end3 cluster is confident; labeled as internal priming if the 3'end of transcript model hit internal priming site but not a Reference annotated 3' end
TES                             ID assigned to each transcription end site, chr_start(0-base)_end_strand according to hg38
PAS                             PAS motif sequence. Serach from the order of "AATAAA", "ATTAAA", "TATAAA", "AGTAAA", "AATATA", "CATAAA", "GATAAA", "AAAAAA", "TTTAAA", "ACTAAA", "AATACA", "AATAGA", "AAGAAA", "AATAAG", "AATAAT", "AATGAA", "AATTAA", "ATTATA"
PAS_score                       the motif score of the PAS motif
PAS_3pos                        the position of the last nt of PAS motif refer to the 3' end of the transcript model
n5_string                       ID of end5 cluster
TSScluster                      the SCAFE defined TSS cluster linked to the end5 cluster
CREID                           the ID of the tCRE linked to the TSS cluster
promoter_type                   promoter-like, enhancer-like, CFCT-alone or unclassed defined from SCREEN (only if genome = hg38, mm10 or mm39)
n5_Reference                    Yes if end5 cluster intersect with any Reference 5' end
n5_support                      a collection of support by Reference 5' end and SCAFE cluster
ORF                         1   length of the best ORF by CPAT
Coding_prob                 1   coding probability of the best ORF by CPAT, 0 if no ORF was found
CPAT_class                  1   coding or non-coding, cutoff: Coding_prob < 0.364
T4_gene_ID              1       final gene ID according to table 4 transcript model
T4_gene_name            1       final gene name according to table 4 transcript model
T4_gene_novelty         1       updated gene novelty
Ref_transcriptClass     1       transcript class inherited from Reference (NA for novel transcript)
Ref_transcriptClass2    1       simplified transcript gene class showing only "protein_coding", "lncRNA" and "others" (NA for novel transcript)
Ref_geneClass           1       gene class inherited from Reference (NA for novel gene)
Ref_geneClass2          1       simplified gene class showing only "protein_coding", "lncRNA" and "others" (NA for novel gene)  
Novel_transcriptClass   1   1   class for novel transcript: "lncRNA"(>200bp), "ncRNA"(<=200bp) and "others" (potentially coding by CPAT)
overall_T_ratio         1   1   ratio of transcript model per gene, according to full-length count (pan-cell-types)
Novel_geneClass         1   1   class for novel gene: "lncRNA"(>50% weighed number of transcript is grouped as non-coding & weighed average transcript length >200bp), "ncRNA"(>50% weighed number of transcript is grouped as non-coding & weighed average transcript length <=200bp) and "others"
Ref_gene_adjust         1       if the 5' and 3' ends of Reference gene are adjusted 
Ref_transcript_adjust   1       if the 5' and 3' ends of Reference transcript are adjusted 
```

# <a name="SALA_cite"></a>Citing SALA
Please cite our preprint when using SALA:  

*CFC-seq: identification of full-length capped RNAs unveil enhancer-derived transcription.*
Chi Wai Yip, Callum Parr, Hazuki Takahashi, Kayoko Yasuzawa, Matthew Valentine, Hiromi Nishiyori-Sueki, Camilla Ugolini, Valeria Ranzani, Mitsuyoshi Murata, Masaki Kato, Wenjing Kang, Wing Hin Yip, Youtaro Shibayama, Andre Darah Sim, Ying Chen, Xufeng Shu, Jonathan Moody, Ramzan Umarov, Manli Yang, Jen-Chien Chang, Luca Pandolfini, Tsugumi Kawashima, Michihira Tagami, Tomoe Nobusada, Tsukasa Kouno, Carlos Alfonso Gonzale, Roberto Albanese, Francesco Dossena, Nejc Haberman, Kokoro Ozaki, Takeya Kasukawa, Boris Lenhard, Martin Frith, Beatrice Bodega, Francesco Nicassio, Lorenzo Calviello, Magda Bienko, Ivano Legnini, Valérie Hilgers, Stefano Gustincich, Jonathan Göke, Charles-Henri Lecellier, Jay W. Shin, Chung-Chau Hon, Piero Carninci
bioRxiv; doi: https://doi.org/10.1101/2024.10.31.620483


# <a name="contribution"></a>Contribution
Major contribution for SALA development: 

1. Main: Chung-Chau Hon, Chi Wai Yip
2. Modification of bambu: Andre Darah Sim, Ying Chen, Jonathan Göke
3. Testing SALA: Hiromi Nishiyori-Sueki, Callum Parr, Camilla Ugolini, Valeria Ranzani 


# <a name="ref"></a>References
* Chen Y, Sim A, Wan YK, Yeo K, Lee JJX, Ling MH, Love MI, Göke J. 2023. Context-aware transcript quantification from long-read RNA-seq data with Bambu. Nat Methods 20: 1187–1195.
* Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, et al. 2021. Twelve years of SAMtools and BCFtools. GigaScience 10: giab008.
* Frith MC, Valen E, Krogh A, Hayashizaki Y, Carninci P, Sandelin A. 2008. A code for transcription initiation in mammalian genomes. Genome Res 18: 1–12.
* Hon C-C, Ramilowski JA, Harshbarger J, Bertin N, Rackham OJL, Gough J, Denisenko E, Schmeier S, Poulsen TM, Severin J, et al. 2017. An atlas of human long non-coding RNAs with accurate 5’ ends. Nature 543: 199–204.
* Leonardi T. 2019. Bedparse: feature extraction from BED files. J Open Source Softw 4: 1228.
* Li H. 2011. Tabix: fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics 27: 718–719.
* Moody J, Kouno T, Chang J-C, Ando Y, Carninci P, Shin JW, Hon C-C. 2022. SCAFE: a software suite for analysis of transcribed cis-regulatory elements in single cells. Bioinforma Oxf Engl 38: 5126–5128.
* Quinlan AR, Hall IM. 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinforma Oxf Engl 26: 841–842.
* The ENCODE Project Consortium, Abascal F, Acosta R, Addleman NJ, Adrian J, Afzal V, Ai R, Aken B, Akiyama JA, Jammal OA, et al. 2020. Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature 583: 699–710.
* Wang L, Park HJ, Dasari S, Wang S, Kocher J-P, Li W. 2013. CPAT: Coding-Potential Assessment Tool using an alignment-free logistic regression model. Nucleic Acids Res 41: e74–e74.
* Wyman D, Mortazavi A. 2019. TranscriptClean: variant-aware correction of indels, mismatches and splice junctions in long-read transcripts. Bioinforma Oxf Engl 35: 340–342.
* Yao L, Wang H, Song Y, Sui G. 2017. BioQueue: a novel pipeline framework to accelerate bioinformatics analysis. Bioinforma Oxf Engl 33: 3286–3288.
* Yip CW, Parr C, Takahashi H, Yasuzawa K, Valentine M, Nishiyori-Sueki H, Ugolini C, Ranzani V, Murata M, Kato M, et al. 2024. CFC-seq: identification of full-length capped RNAs unveil enhancer-derived transcription. BioRvix DOI:doi/10.1101/2024.10.31.620483.










