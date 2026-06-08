# List of required packages
required_packages <- c(
  "Biostrings", "GenomicRanges", "Rsamtools",
  "data.table", "dplyr", "magrittr", "tidyr"
)

# Function to check and install missing packages
install_if_missing <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse=", "))
    BiocManager::install(missing_packages, ask=FALSE)
  }
}

# Ensure BiocManager is installed for Bioconductor packages
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}

# Install missing packages
install_if_missing(required_packages)

# load packages
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyr))

#SALA_directory="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/demo_output_local/sala/transcript/Neuron_series_demo"
#out_prefix="Neuron_series_demo"
#resource_directory="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/resources"
#ref_directory="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/resources/GENCODE_V39"
#fasta_file="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
#read_per_rep_isoform=5
#read_per_rep_novel=1
#n5_confid="Yes"
#SALA_gene_path="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/demo_output_local/sala/gene/table0_gene/Neuron_series_demo "
#sample_info="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/demo_output_local/input/sample.txt"
#SCAFE_directory="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/demo_output_local/scafe/aggregate/run_full/out"
#isoform_ratio=0.1
#keep_internal_prime="No"
#n3_confid_isoform="Yes"
#n3_confid_novel="No"
#cpat_path="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/demo_output_local/sala/transcript/Neuron_series_demo/cpat"
#===========

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 14) {
  stop("Usage: Rscript SALA.filter.R <SALA_directory> <out_prefix> <resource_directory> <ref_directory> <fasta_file> <read.per.rep_ref.novel.Tx> <read.per.rep_non-ref.novel.Tx> <require.5'.confidence> <SALA_gene_path> <sample_file> <SCAFE_directory> <isoform_ratio> <genome> <keep_internal_prime> <CPAT_path(optional)>\n\n",
    "SALA_directory                <required>	path of the folder of SALA transcript annotation output\n",
	"out_prefix                    <required>	output files prefix\n",
	"resource_directory            <required>	path of the resources folder of SALA\n",
	"ref_directory                 <required>	path of the folder containing the infomation of reference transcriptome\n",
	"fasta_file                    <required>	genome fasta file\n",
	"read.per.rep_ref.novel.Tx     <required>	number of reads per replicate from a sample for novel isoform of known gene\n",
	"read.per.rep_non-ref.novel.Tx <required>	number of reads per replicate from a sample for novel transcript of novel gene\n",
	"isoform_ratio                 <required>	the ratio of a novel isoform across all the transcript in a gene, according to the full length read count\n",
	"require.5'.confidence         <required>	if Yes, all the novel transcripts are required to have their 5’ ends located inside confident TSS clusters\n",
	"SALA_gene_path                <required>	path of the folder of SALA gene annotation output\n",
	"sample_file                   <required>	txt file for the input library: column1, library prefix; column2, sample ID; column3, sampleID with replicate ID\n",
	"SCAFE_directory               <required>	path of the folder of SCAFE output (NA if SCAFE is skipped)\n",
	"genome                        <required>   genome name: hg38, mm10, mm39 or NA",
	"keep_internal_prime           <required>   Yes or No (designed by poly-dT internal priming potential, set Yes if poly-dT priming is not used)",
	"n3_confid_ref.novel.Tx        <required>	Yes or No (filter for transcripts that are supported by ex3_cluster for novel isoform of known gene)\n",
	"n3_confid_non-ref.novel.Tx    <required>	Yes or No (filter for transcripts that are supported by ex3_cluster for novel transcript of novel gene)\n",
	"CPAT_path                     <optional>	path of the folder expected for CPAT result")}

# Assign arguments to variables
SALA_directory <- args[1]
out_prefix <- args[2]
resource_directory <- args[3]
ref_directory <- args[4]
fasta_file <- args[5]
read_per_rep_isoform <- args[6]
read_per_rep_novel <- args[7]
isoform_ratio <- args[8]
n5_confid <- args[9]
SALA_gene_path <- args[10]
sample_info <- args[11]
SCAFE_directory <- args[12]
genome <- args[13]
keep_internal_prime <- args[14]
keep_internal_prime <- args[14]
n3_confid_isoform <- args[15]
n3_confid_novel <- args[16]

#cpat_path <- args[17]

#path
path1 <- paste0(SALA_directory, "/bed/")
path2 <- paste0(SALA_directory, "/log/")
path3 <- paste0(SALA_directory, "/tmp/")
SALA_gene_info <- paste0(SALA_gene_path,"/log/",out_prefix,".model.info.tsv.gz")
SALA_gene_bed <- paste0(SALA_gene_path,"/bed/",out_prefix,".gene.bed.bgz")

bedtools_bin <- "bedtools"
samtools_bin <- "samtools"
tabix_bin <- "tabix"
bgzip_bin <- "bgzip"

#===========
#prepare 3'end and 5' end bed6 bed file from transcript models

table0_model <- read.delim(paste0(path1,out_prefix,".model.bed.bgz"), header=F, stringsAsFactors = F, check.names = F)
options(scipen=999)
table0_model3n <- table0_model
table0_model5n <- table0_model
table0_model5n$V3[which(table0_model5n$V6=="+")] <- table0_model5n$V2[which(table0_model5n$V6=="+")]+1
table0_model5n$V2[which(table0_model5n$V6=="-")] <- table0_model5n$V3[which(table0_model5n$V6=="-")]-1
write.table(table0_model5n[order(table0_model5n$V1,table0_model5n$V2),c(1:6)],gzfile(paste0(path1,out_prefix,".model.5n.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
table0_model3n$V3[which(table0_model3n$V6=="-")] <- table0_model3n$V2[which(table0_model3n$V6=="-")]+1
table0_model3n$V2[which(table0_model3n$V6=="+")] <- table0_model3n$V3[which(table0_model3n$V6=="+")]-1
write.table(table0_model3n[order(table0_model3n$V1,table0_model3n$V2),c(1:6)],gzfile(paste0(path1,out_prefix,".model.3n.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
rm(table0_model, table0_model3n, table0_model5n)

#===========
files=list.files(path=SCAFE_directory, pattern="CRE.coord.bed.gz", recursive=T)
out_prefix2=sapply(strsplit(files[1],"\\/"),"[",2)

system(paste0(bedtools_bin," bed12tobed6 -i ",path1,out_prefix,".model.bed.bgz | gzip > ", path1,out_prefix,".model.bed6.bed.gz"))
system(paste0(bedtools_bin," intersect -s -wa -wb -a ",path1,out_prefix,".end3.bed.bgz -b ",ref_directory,"/n3.bed.gz | gzip > ", path1,out_prefix,".end3.cluster.region.ref3n.bed.gz"))
system(paste0(bedtools_bin," intersect -s -wa -wb -a ",path1,out_prefix,".end5.bed.bgz -b ",ref_directory,"/n5.bed.gz | gzip > ", path1,out_prefix,".end5.cluster.region.ref5n.bed.gz"))
if (SCAFE_directory != "NA"){
system(paste0(bedtools_bin," intersect -s -wa -wb -a ",path1,out_prefix,".end5.bed.bgz -b ",SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".cluster.coord.bed.gz | gzip > ", path1,out_prefix,".end5.cluster.region.cluster.bed.gz"))

if (genome == "hg38"){
system(paste0(bedtools_bin," closest -a ",SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.bed.gz -b  ",resource_directory,"/CRE.from.SCREEN.encodev3/GRCh38-ELS.all.enhancer.sort.bed -D a> ", SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.e.bed"))
system(paste0(bedtools_bin," closest -a ",SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.bed.gz -b  ",resource_directory,"/CRE.from.SCREEN.encodev3/GRCh38-PLS.all.promoter.sort.bed -D a> ", SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.p.bed"))
system(paste0(bedtools_bin," closest -a ",SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.bed.gz -b  ",resource_directory,"/CRE.from.SCREEN.encodev3/GRCh38-CTCF.sort.bed -D a> ", SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.c.bed"))
system(paste0(bedtools_bin," closest -a ",SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.bed.gz -b  ",resource_directory,"/F5_enhancer/hg38_robust_enhancers.sort.bed -D a> ", SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.andersson.robust.e.bed"))
system(paste0(bedtools_bin," closest -a ",SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.bed.gz -b  ",resource_directory,"/F5_enhancer/hg38_permissive_enhancers.sort.bed -D a> ", SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.andersson.permissive.e.bed"))
}else if (genome == "mm10"){
system(paste0(bedtools_bin," closest -a ",SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.bed.gz -b  ",resource_directory,"/CRE.from.SCREEN.encodev3/mm10-ELS.all.sort.bed -D a> ", SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.e.bed"))
system(paste0(bedtools_bin," closest -a ",SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.bed.gz -b  ",resource_directory,"/CRE.from.SCREEN.encodev3/mm10-PLS.all.sort.bed -D a> ", SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.p.bed"))
system(paste0(bedtools_bin," closest -a ",SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.bed.gz -b  ",resource_directory,"/CRE.from.SCREEN.encodev3/mm10-CTCF.sort.bed -D a> ", SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.c.bed"))
}else if (genome == "mm39"){
system(paste0(bedtools_bin," closest -a ",SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.bed.gz -b  ",resource_directory,"/CRE.from.SCREEN.encodev3/mm39-ELS.all.sort.bed -D a> ", SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.e.bed"))
system(paste0(bedtools_bin," closest -a ",SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.bed.gz -b  ",resource_directory,"/CRE.from.SCREEN.encodev3/mm39-PLS.all.sort.bed -D a> ", SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.p.bed"))
system(paste0(bedtools_bin," closest -a ",SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.bed.gz -b  ",resource_directory,"/CRE.from.SCREEN.encodev3/mm39-CTCF.sort.bed -D a> ", SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.c.bed"))
}else {print("promoter type intersection skipping due to unmatched genome")}}

if (SCAFE_directory == "NA"){
if (genome == "hg38"){system(paste0(bedtools_bin," intersect -s -wa -wb -a ",path1,out_prefix,".end5.bed.bgz -b ",resource_directory,"/CTSS_cluster/hg38_fair+new_CAGE_peaks_phase1and2.bed.bgz | gzip > ", path1,out_prefix,".end5.cluster.region.cluster.bed.gz"))}
if (genome == "mm10"){system(paste0(bedtools_bin," intersect -s -wa -wb -a ",path1,out_prefix,".end5.bed.bgz -b ",resource_directory,"/CTSS_cluster/mm10_fair+new_CAGE_peaks_phase1and2.bed.bgz | gzip > ", path1,out_prefix,".end5.cluster.region.cluster.bed.gz"))}
}

system(paste0(samtools_bin, " faidx ", fasta_file))
system(paste0(bedtools_bin, " getfasta -s -nameOnly -split -fi ", fasta_file, " -bed ", path1, out_prefix, ".model.bed.bgz > ", path1, out_prefix, ".model.fasta"))
#===========
# Load reference genome (FASTA file)
#indexFa(fasta_file)
genome_seq <- FaFile(fasta_file)
open(genome_seq)
chrom_sizes <- seqlengths(genome_seq)
chrom_table <- data.frame(Chromosome = names(chrom_sizes), Length = as.numeric(chrom_sizes))

print("Internal priming prediction")
t0bed12 <- read.delim(paste0(path1,out_prefix,".model.bed.bgz"), header=F, stringsAsFactors=F)
t0bed12$V2[which(t0bed12$V6 == "+")] <- t0bed12$V3[which(t0bed12$V6 == "+")]-1
t0bed12$V3[which(t0bed12$V6 == "-")] <- t0bed12$V2[which(t0bed12$V6 == "-")]+1
t0bed12$V4 <- paste0(t0bed12$V1,"_",t0bed12$V2,"_",t0bed12$V3,"_",t0bed12$V6)
t0bed12 <- t0bed12%>%group_by(V4)%>%dplyr::mutate(V5=n())
t0bed12a3 <- unique(t0bed12[,c(1:6)])
write.table(t0bed12a3[order(t0bed12a3$V1,t0bed12a3$V2),],gzfile(paste0(path1,out_prefix,".model.n3.bed6.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

t0bed12a3$V2 <- t0bed12a3$V2-50
t0bed12a3$V3 <- t0bed12a3$V3+50
t0bed12a3 <- t0bed12a3[which(t0bed12a3$V2>=0),]
t0bed12a3 <- left_join(t0bed12a3,chrom_table,by=c("V1"="Chromosome"),copy=F)
t0bed12a3 <- t0bed12a3[which(t0bed12a3$V3 < t0bed12a3$Length),]
t0bed12a3 <- t0bed12a3[order(t0bed12a3$V1,t0bed12a3$V2),]
write.table(t0bed12a3, gzfile(paste0(path1,out_prefix,".model.n3.n101.bed6.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)


# Convert BED data to GRanges
gr <- GRanges(seqnames = t0bed12a3$V1, ranges = IRanges(start = t0bed12a3$V2+1, end = t0bed12a3$V3), strand = t0bed12a3$V6)
names(gr) <- t0bed12a3$V4
# Extract sequences
extracted_seqs <- getSeq(genome_seq, gr)

extracted_df <- data.frame(TES_ID = names(gr), sequence = as.character(extracted_seqs))
extracted_df <- left_join(extracted_df, t0bed12a3[,c(4,5)], by=c("TES_ID"="V4"), copy=F)
colnames(extracted_df)[3] <- "TES_count"

for (i in 1:101){extracted_df[,i+3]=substr(extracted_df$sequence, start = i, stop = i)}
extracted_df$fracA08=(rowSums(extracted_df[,c(55:62)] == "A"))/8
extracted_df$fracA16=(rowSums(extracted_df[,c(55:70)] == "A"))/16
extracted_df$fracA20=(rowSums(extracted_df[,c(55:74)] == "A"))/20
extracted_df$internal_prime="no"
extracted_df$internal_prime[which(extracted_df$fracA08 > 0.75)]="fracA08"
extracted_df$internal_prime[which(extracted_df$fracA16 > 0.5)]="fracA16"
with_IP=100-(sum(extracted_df$TES_count[which(extracted_df$internal_prime == "no")])/sum(extracted_df$TES_count))*100
print(paste0("% of transcript model with potential internal priming: ", signif(with_IP,3),"%"))
extracted_df=extracted_df[,c(1:3,105:108)]
write.table(extracted_df,gzfile(paste0(path2,"potential_internal_prime_TES.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
rm(t0bed12)

#====================
# PAS search
print("PAS search")
gr <- GRanges(seqnames = t0bed12a3$V1, ranges = IRanges(start = t0bed12a3$V2+1, end = t0bed12a3$V3), strand = t0bed12a3$V6)
names(gr) <- t0bed12a3$V4
# Extract sequences
extracted_seqs <- getSeq(genome_seq, gr)

extracted_df <- data.frame(TES_ID = names(gr), sequence = as.character(extracted_seqs))
extracted_df <- left_join(extracted_df, t0bed12a3[,c(4,5)], by=c("TES_ID"="V4"), copy=F)
colnames(extracted_df)[3] <- "TES_count"
# take window from -41 to -5 (-35 to -5 refer to the last nt of PAS motif)
extracted_df$sequence <- substr(extracted_df$sequence, 10, 46)

motifs <- c("AATAAA", "ATTAAA", "TATAAA", "AGTAAA", "AATATA", "CATAAA", "GATAAA", "AAAAAA", "TTTAAA", "ACTAAA", "AATACA", "AATAGA", "AAGAAA", "AATAAG", "AATAAT", "AATGAA", "AATTAA", "ATTATA")
motifs_score <- c(7.792203, 6.470447, 5.084153, 4.861009, 4.336938, 3.985540, 3.985540, 3.921002, 3.762397, 3.474715, 3.238326, 3.238326, 3.227855, 3.207235, 3.207235, 3.207235, 3.207235, 3.015182)
motif_df <- data.frame(motifs = motifs, score = motifs_score, stringsAsFactors = FALSE)

df_long <- do.call(rbind, lapply(motif_df$motifs, function(motif) {
  hits <- gregexpr(motif, extracted_df$sequence, fixed = TRUE)
  
  data.frame(
    TES_ID = extracted_df$TES_ID,
    sequence = extracted_df$sequence,
    motifs = motif,
    start = sapply(hits, function(x) if (x[1] == -1) NA else x[1])
  )
}))
df_long <- left_join(df_long, motif_df, by="motifs",copy=F)
df_long <- df_long[which(!is.na(df_long$start)),]
df_long1 <- df_long%>%group_by(TES_ID)%>%slice_max(score)
df_long1 <- df_long1%>%group_by(TES_ID)%>%slice_min(abs(15-start))
df_long1 <- df_long1%>%group_by(TES_ID)%>%slice_max(start)

df_long1$PAS_3pos <- df_long1$start-42+5
colnames(df_long1)[c(3,5)] <- c("PAS", "PAS_score")
write.table(df_long1[,c(1,3,5,6)],gzfile(paste0(path2,"PAS_TES.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

rm(extracted_df, t0bed12a3)

#====================
print("Add gene group and full-length read count")
ref_info <- fread(paste0(ref_directory,"/transcript_to_gene.tsv"), header = FALSE, stringsAsFactors = FALSE)
transcript_info <- fread(paste0(path2,out_prefix, ".model.info.tsv.gz"), header = TRUE, stringsAsFactors = FALSE)
transcript_info <- data.frame(transcript_info)

#keep only major chromosome
transcript_info$chr=sapply(strsplit(transcript_info$loc,":"),"[",1)
transcript_info=transcript_info[which(nchar(transcript_info$chr)<=5),c(1:18)]

gene_info <- fread(SALA_gene_info, header=T, stringsAsFactors = F, check.names = F)

colnames(gene_info)[c(7,8)]=c("IN1_gene_ID","IN1_gene_name")
transcript_info=left_join(transcript_info,gene_info[,c(1,7,8)],by="model_ID",copy=F)

transcript_info$gene_novelty="Novel"
transcript_info$gene_novelty[which(transcript_info$IN1_gene_ID %in% unique(ref_info$V2))]="Ref"
transcript_info$transcript_novelty="Novel"
transcript_info$transcript_novelty[which(transcript_info$model_ID %in% ref_info$V1)]="Ref"

#====================
#add read per replicate and convert library prefix to sample_rep
sample_info1 <- fread(sample_info, header = TRUE, stringsAsFactors = FALSE)
sample_info1 <- sample_info1 %>% group_by(sample_ID, sample_rep) %>% dplyr::mutate(count=n())
sample_info1a <- sample_info1[which(sample_info1$count == 1),]
sample_info1b <- sample_info1[which(sample_info1$count > 1),]

read_info3 <- read.delim(paste0(path2,"full_length_support_matrix.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
matched_indices <- match(sample_info1a$library_prefix, colnames(read_info3))
colnames(read_info3)[matched_indices] <- sample_info1a$sample_rep
num <- ncol(read_info3)
if (nrow(sample_info1b>0)){
  need <- unique(sample_info1b$sample_rep)
  for (i in 1:length(need)){
    label <- sample_info1b$library_prefix[which(sample_info1b$sample_rep == need[i])]
    read_info3$new <- rowSums(read_info3[,which(colnames(read_info3) %in% label)])
    colnames(read_info3)[num+i] <- need[i]}}

label1 <- c(colnames(read_info3)[1],unique(sample_info1$sample_rep))
read_info3 <- read_info3[,label1]

#add read count per sample
sample_info2 <- unique(sample_info1[,c(2,3)])%>%group_by(sample_ID)%>%dplyr::mutate(count=n())
need2 <- unique(sample_info2$sample_ID)
num2 <- ncol(read_info3)
for (i in 1:length(need2)){
  label <- sample_info2$sample_rep[which(sample_info2$sample_ID == need2[i])]
  read_info3$new <- rowSums(read_info3[,which(colnames(read_info3) %in% label)])
  colnames(read_info3)[num2+i] <- need2[i]}

transcript_info <- left_join(transcript_info,read_info3, by="model_ID",copy=F)
transcript_info[is.na(transcript_info)]=0

transcript_info$isoform_filter <- "permissive"
for (i in 1:length(need2)){
  label2 <- sample_info2$sample_rep[which(sample_info2$sample_ID == need2[i])]
  transcript_info$isoform_filter[which(rowSums(transcript_info[,which(colnames(transcript_info) %in% label2)] >= as.numeric(read_per_rep_isoform)) >=2)] <- "standard"}
stat1 <- transcript_info%>%group_by(isoform_filter)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percentage=count/sum(count))
stat1$label=paste0(stat1$isoform_filter,": ",stat1$count," (",signif(stat1$percentage*100,3),"%)")
#print(paste0("[Isoform_filter]: ",stat1$label))

transcript_info$novel_gene_Tx_filter <- "permissive"
for (i in 1:length(need2)){
  label2 <- sample_info2$sample_rep[which(sample_info2$sample_ID == need2[i])]
  transcript_info$novel_gene_Tx_filter[which(rowSums(transcript_info[,which(colnames(transcript_info) %in% label2)] >= as.numeric(read_per_rep_novel)) >=2)] <- "standard"}
stat2 <- transcript_info%>%group_by(novel_gene_Tx_filter)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percentage=count/sum(count))
stat2$label=paste0(stat2$novel_gene_Tx_filter,": ",stat2$count," (",signif(stat2$percentage*100,3),"%)")
#print(paste0("[Raw table transcript model content]: ",stat2$label))

transcript_info$ref_source <- "novel_transcript"
transcript_info$ref_source[which(transcript_info$model_ID %in% ref_info$V1)] <- "fulllength_ref"
transcript_info$ref_source[which(transcript_info$ref_source == "fulllength_ref" & transcript_info$full_qry_count==0 & transcript_info$partial_qry_count>0)] <- "partial_ref"
transcript_info$ref_source[which(transcript_info$ref_source == "fulllength_ref" & transcript_info$full_qry_count==0 & transcript_info$partial_qry_count==0)] <- "non_detectable_ref"
stat3 <- transcript_info%>%group_by(ref_source)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percentage=count/sum(count))
stat3$label=paste0(stat3$ref_source,": ",stat3$count," (",signif(stat3$percentage*100,3),"%)")
print(paste0("[Raw table transcript model content]: ",stat3$label))

#=====
##add transcript ratio##
print("add transcript ratio, exon number & transcript length")
need3=paste0("Tx_ratio_",need2)
for (i in 1: length(need3)){
  transcript_info <- transcript_info%>%group_by(IN1_gene_ID)%>%dplyr::mutate(!!sym(need3[i]) := !!sym(need2[i]) / sum(!!sym(need2[i])))
  transcript_info[which(transcript_info[,ncol(transcript_info)] == "NaN"),ncol(transcript_info)] <- 0
}
transcript_info$max_T_ratio <- do.call(pmax, transcript_info[,c(need3)])
#=====
#add transcript length and exon number
bed6 <- read.delim(paste0(path1,out_prefix,".model.bed6.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
bed6$length <- bed6$V3-bed6$V2
bed6a <- bed6%>%group_by(V4)%>%dplyr::summarise(n_exon=n(),transcript_length=sum(length))
transcript_info <- left_join(transcript_info, bed6a, by=c("model_ID"="V4"), copy=F)

#=====
#add n3 string
data <- separate_rows(transcript_info[,c(1,18)],full_set_bound_str,sep="_")
data <- data[grep("T",data$full_set_bound_str),]
colnames(data)[2]="n3_string"
transcript_info <- left_join(transcript_info,data,by="model_ID",copy=F)

#=====
print("add 3' end feature & internal priming")
internal_prime_sample <- read.delim(paste0(path2,"potential_internal_prime_TES.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
n3_ref <- read.delim(paste0(ref_directory,"/n3.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
internal_prime_sample$internal_prime[which(internal_prime_sample$internal_prime=="no")] <- "No"
internal_prime_sample$internal_prime[which(internal_prime_sample$internal_prime!="No")] <- "Yes"
internal_prime_sample$Reference <- "No"
internal_prime_sample$Reference[which(internal_prime_sample$TES_ID%in% n3_ref$V4)] <- "Yes"
internal_prime_sample <- internal_prime_sample[,c(1,8,7)]
colnames(internal_prime_sample) <- c("label","Reference","internal_priming")

table0_model3n <- read.delim(paste0(path1,out_prefix,".model.3n.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
table0_model3n$TES <- paste0(table0_model3n$V1,"_",table0_model3n$V2,"_",table0_model3n$V3,"_",table0_model3n$V6)
table0_model3n <- left_join(table0_model3n,internal_prime_sample, by=c("TES"="label"),copy=F)

table0_model3n$internal_priming[which(is.na(table0_model3n$internal_priming))] <- "No"
table0_model3n$Reference[which(is.na(table0_model3n$Reference))] <- "No"
table0_model3n$internal_priming[which(table0_model3n$Reference == "Yes" & table0_model3n$internal_priming == "Yes")] <- "Yes;Reference"
transcript_info <- left_join(transcript_info, table0_model3n[c(4,9)], by=c("model_ID" = "V4"), copy=F)
rm(internal_prime_sample,n3_ref,table0_model3n)

#================================
table0_model2 <- read.delim(paste0(path1,out_prefix,".end3.cluster.region.ref3n.bed.gz"), header=F, stringsAsFactors = F, check.names = F) 
transcript_info$n3_Reference <- "No"
transcript_info$n3_Reference[which(transcript_info$n3_string %in% unique(table0_model2$V4))] <- "Yes"
#transcript_info%>%group_by(n3_Reference)%>%summarise(count=n())

transcript_info$n3_support <- "cluster"
transcript_info$n3_support[grep("XT",transcript_info$n3_string)] <- "non_cluster"
#transcript_info%>%group_by(internal_priming,n3_support)%>%dplyr::summarise(count=n())

data0 <- transcript_info[,c("model_ID","n3_Reference","n3_support")]
data0$n3_support[which(data0$n3_support=="cluster")] <- "n3_cluster"
data0$n3_support[which(data0$n3_support=="non_cluster")] <- NA
data0$n3_Reference[which(data0$n3_Reference == "Yes")] <- "Reference"
data0$n3_Reference[which(data0$n3_Reference == "No")] <- NA
data0 <- reshape2::melt(data0, id=1)
data0 <- data0[which(!is.na(data0$value)),]
data01 <- data0%>%group_by(model_ID)%>%dplyr::summarise(n3_support=paste(value,collapse=" & "))
transcript_info <- transcript_info[,-which(colnames(transcript_info) == "n3_support")]

transcript_info <- left_join(transcript_info, data01, by="model_ID", copy=F)
transcript_info$n3_support[which(is.na(transcript_info$n3_support))] <- "no_support"
#transcript_info$n3_support[which(transcript_info$internal_priming == "Yes")] <- "internal_priming"
transcript_info%>%group_by(n3_support)%>%dplyr::summarise(count=n(), .groups="drop_last")

table0_model3n <- read.delim(paste0(path1,out_prefix,".model.3n.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
table0_model3n$TES <- paste0(table0_model3n$V1,"_",table0_model3n$V2,"_",table0_model3n$V3,"_",table0_model3n$V6)
table0_model3n <- left_join(table0_model3n,df_long1[,c("TES_ID","PAS","PAS_score","PAS_3pos")], by=c("TES"="TES_ID"),copy=F)
transcript_info <- left_join(transcript_info, table0_model3n[,c("V4","TES","PAS","PAS_score","PAS_3pos")], by=c("model_ID" = "V4"), copy=F)

rm(table0_model3n, df_long1)

#===================================================================================================================
if (SCAFE_directory != "NA"){
CRE <- read.delim(paste0(SCAFE_directory,"/annotate/",out_prefix2,"/log/",out_prefix2,".CRE.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
if (genome %in% c("hg38","mm10","mm39")){
print("Perform SCAFE - SCREEN connection")
CRE.p <- read.delim(paste0(SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.p.bed"), header=F, stringsAsFactors = F, check.names = F)
CRE.e <- read.delim(paste0(SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.e.bed"), header=F, stringsAsFactors = F, check.names = F)
CRE.ctcf <- read.delim(paste0(SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.c.bed"), header=F, stringsAsFactors = F, check.names = F)

CREp=unique(CRE.p$V4[which(CRE.p$V19 ==0)])
CREe=unique(CRE.e$V4[which(CRE.e$V19 ==0)])
CREc=unique(CRE.ctcf$V4[which(CRE.ctcf$V19 ==0)])

CRE$promoter=0
CRE$promoter[which(CRE$CREID %in% CREp)]=1
CRE$enhancer=0
CRE$enhancer[which(CRE$CREID %in% CREe)]=1
CRE$CTCF=0
CRE$CTCF[which(CRE$CREID %in% CREc)]=1

CRE$promoter_type="unclassed"
CRE$promoter_type[which(CRE$CTCF ==1)]="CTCF-alone"
CRE$promoter_type[which(CRE$enhancer ==1)]="enhancer-like"
CRE$promoter_type[which(CRE$promoter ==1)]="promoter-like"

###sample-specific bi-directional enhancer###
CREdir=read.delim(paste0(SCAFE_directory,"/directionality/",out_prefix2,"/log/",out_prefix2,".directionality.log.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CREdir$fwd_rev_count=CREdir$fwd_count+CREdir$rev_count
CREdir$orientation[which(abs(CREdir$directionality)<=0.2)]="unidirectional"
CREdir$orientation[which(CREdir$orientation == "divergent")]="bidirectional"
CRE=left_join(CRE,CREdir[,c(1,4,5,7,8,13,10)],by="CREID",copy=F)
write.table(CRE, paste0(SCAFE_directory,"/annotate/",out_prefix2,"/log/",out_prefix2,".CRE.info.p.e.se.tsv"), col.names=T, row.names=F, sep="\t", quote=F)}

if (genome == "hg38"){
re=read.delim(paste0(SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.andersson.robust.e.bed"), header=F, stringsAsFactors = F, check.names = F)
pe=read.delim(paste0(SCAFE_directory,"/annotate/",out_prefix2,"/bed/",out_prefix2,".CRE.coord.andersson.permissive.e.bed"), header=F, stringsAsFactors = F, check.names = F)
re_e=unique(re$V4[which(re$V25 ==0)])
pe_e=unique(pe$V4[which(pe$V25 ==0)])
CRE$Andersson_robust=0
CRE$Andersson_robust[which(CRE$CREID %in% re_e)]=1
CRE$Andersson_permissive=0
CRE$Andersson_permissive[which(CRE$CREID %in% pe_e)]=1
write.table(CRE, paste0(SCAFE_directory,"/annotate/",out_prefix2,"/log/",out_prefix2,".CRE.info.p.e.se.tsv"), col.names=T, row.names=F, sep="\t", quote=F)}}

#==================================================
print("Add n5 cluster annotation & promoter-type")
#Use cluster to cluster intersect and link to tCRE
transcript_info$n5_string <- sapply(strsplit(transcript_info$full_set_bound_str,"_"),"[",1)
table0_model4 <- read.delim(paste0(path1,out_prefix,".end5.cluster.region.cluster.bed.gz"), header=F, stringsAsFactors = F, check.names = F) 



if (SCAFE_directory != "NA"){
cluster <- read.delim(paste0(SCAFE_directory,"/annotate/",out_prefix2,"/log/",out_prefix2,".cluster.info.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
table0_model4 <- left_join(table0_model4,cluster[,c(1,16)],by=c("V16"="clusterID"),copy=F)
table0_model4a <- unique(table0_model4[,c(4,16,25)])%>%group_by(V4)%>%dplyr::summarise(TSScluster=paste(unique(V16), collapse=";"),
                                                                                    CREID=paste(unique(CREID), collapse=";"))

if (genome %in% c("hg38","mm10","mm39")){
table0_model4a <- left_join(table0_model4a,CRE[,c("CREID","promoter_type")], by="CREID",copy=F)}}

if (SCAFE_directory == "NA"){
table0_model4a <- unique(table0_model4[,c(4,16,22)])%>%group_by(V4)%>%dplyr::summarise(TSScluster=paste(unique(V16), collapse=";"),
                                                                                    CREID=paste(unique(V22), collapse=";"))
if (genome == "hg38"){CRE=read.delim(paste0(resource_directory,"/CTSS_cluster/hg38.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table0_model4a <- left_join(table0_model4a,CRE[,c("CREID","promoter_type")], by="CREID",copy=F)}
if (genome == "mm10"){CRE=read.delim(paste0(resource_directory,"/CTSS_cluster/mm10.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table0_model4a <- left_join(table0_model4a,CRE[,c("CREID","promoter_type")], by="CREID",copy=F)}
}


transcript_info <- left_join(transcript_info, table0_model4a, by=c("n5_string"="V4"),copy=F)

table0_model5 <- read.delim(paste0(path1,out_prefix,".end5.cluster.region.ref5n.bed.gz"), header=F, stringsAsFactors = F, check.names = F) 
transcript_info$n5_Reference <- "no"
transcript_info$n5_Reference[which(transcript_info$n5_string %in% unique(table0_model5$V4))] <- "yes"
#transcript_info%>%group_by(n5_Reference,promoter_type)%>%summarise(count=n())

data1 <- transcript_info[,c("model_ID","CREID","n5_Reference")]
data1$CREID[!is.na(data1$CREID)] <- "SCAFE"
data1$n5_Reference[which(data1$n5_Reference == "yes")] <- "Reference"
data1$n5_Reference[which(data1$n5_Reference == "no")] <- NA
data1 <- reshape2::melt(data1, id=1)
data1 <- data1[which(!is.na(data1$value)),]
data2 <- data1%>%group_by(model_ID)%>%dplyr::summarise(n5_support=paste(value,collapse=" & "), .groups="drop_last")
transcript_info <- left_join(transcript_info, data2, by="model_ID", copy=F)
transcript_info$n5_support[which(is.na(transcript_info$n5_support))] <- "no_support"
transcript_info%>%group_by(n5_support)%>%dplyr::summarise(count=n(), .groups="drop_last")
#write.table(transcript_info, gzfile("table0.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#==================
#coding potential by CPAT
if (length(args) >= 17 && nchar(args[17]) > 0) {
  cpat_path <- args[17]
  
  if (!file.exists(paste0(cpat_path,"/output.ORF_prob.best.tsv"))) {
    print("Running CPAT...")
    if (!dir.exists(cpat_path)) {dir.create(cpat_path)}
    system(paste0("cpat -x ",resource_directory,"/CPAT/Human_Hexamer.tsv -d ",resource_directory,"/CPAT/Human_logitModel.RData  --top-orf=5 -g ",path1,out_prefix,".model.fasta -o ", cpat_path, "/output --log-file ", cpat_path, "/CPAT_run_info.log >/dev/null 2>&1"))
    print("CPAT finished")
  } else {
    message("CPAT result already exists. Skip running CPAT.")
  }
  
  CPAT <- read.delim(paste0(cpat_path,"/output.ORF_prob.best.tsv"), header=T, stringsAsFactors=F, check.names=F)
  CPAT$seq_ID <- sapply(strsplit(CPAT$seq_ID, "\\("), "[", 1)
  transcript_info <- left_join(transcript_info, CPAT[, c(1, 8, 11)], by=c("model_ID"="seq_ID"), copy=F)
  transcript_info$Coding_prob[is.na(transcript_info$Coding_prob)] <- 0
  transcript_info$CPAT_class <- "coding"
  transcript_info$CPAT_class[which(transcript_info$Coding_prob < 0.364)] <- "non-coding"
  transcript_info %>% group_by(CPAT_class) %>% dplyr::summarise(count=n())
} else {
  message("Skipping CPAT integration (<cpat_path> not provided)")
}

#==================
#generate raw & filtered table
write.table(transcript_info, gzfile(paste0(path2,out_prefix,".table0_raw.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
print(paste0(out_prefix,".table0_raw.tsv.gz exported to ",path2))
#remove internal priming
if (keep_internal_prime == "No") {
  transcript_info1 <- transcript_info[-which(transcript_info$internal_priming == "Yes" & transcript_info$ref_source == "novel_transcript"),]
  print("potential internal primed transcrpts are removed")
} else {
  transcript_info1=transcript_info
  print("potential internal primed transcrpts are kept")
}

#n5_confid 
if (n5_confid=="Yes"){transcript_info1=transcript_info1[which(transcript_info1$n5_support != "no_support"),]}

#n3_confid
if (n3_confid_isoform=="Yes"){transcript_info1=transcript_info1[-which(transcript_info1$n3_support == "no_support" & transcript_info1$gene_novelty == "Ref" & transcript_info1$ref_source == "novel_transcript"),]}
if (n3_confid_novel=="Yes"){transcript_info1=transcript_info1[-which(transcript_info1$n3_support == "no_support" & transcript_info1$gene_novelty == "Novel" & transcript_info1$ref_source == "novel_transcript"),]}


#filter table by read count
transcript_info_finala=transcript_info1[which(transcript_info1$ref_source!="novel_transcript"),]

transcript_info_finalb=transcript_info1[which(transcript_info1$isoform_filter=="standard" & transcript_info1$max_T_ratio >= as.numeric(isoform_ratio) & transcript_info1$gene_novelty == "Ref" & transcript_info1$ref_source == "novel_transcript"),]
transcript_info_finalc=transcript_info1[which(transcript_info1$novel_gene_Tx_filter=="standard" & transcript_info1$gene_novelty == "Novel"),]
transcript_info_final=rbind(transcript_info_finala,transcript_info_finalb,transcript_info_finalc)

write.table(transcript_info_final, gzfile(paste0(path2,out_prefix,".table4_filtered.All_Ref.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
print(paste0(out_prefix,".table4_filtered.All_Ref.tsv.gz exported to ",path2))

#bed12 for filtered table
table0.bed12=read.delim(paste0(path1,out_prefix,".model.bed.bgz"), header=F, stringsAsFactors = F, check.names = F)
table4.bed12=table0.bed12[which(table0.bed12$V4 %in% transcript_info_final$model_ID),]
write.table(table4.bed12[order(table4.bed12$V1,table4.bed12$V2),], gzfile(paste0(path1,out_prefix,".table4.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
system(paste0("gunzip -c ",path1,out_prefix,".table4.bed12.bed.gz | ",bgzip_bin, " > ",path1,out_prefix,".table4.bed12.bed.bgz"))
system(paste0(tabix_bin," -p bed ",path1,out_prefix,".table4.bed12.bed.bgz"))

print(paste0(out_prefix,".table4.bed12.bed.bgz exported to ",path1,". It is used for gene annotation."))
