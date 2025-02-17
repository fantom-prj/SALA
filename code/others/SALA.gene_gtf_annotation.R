# List of required packages
required_packages <- c(
  "data.table", "dplyr", "magrittr", "tidyr" ,"stringr"
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
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))

#SALA_directory="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/demo_output_local/sala/transcript/Neuron_series_demo"
#out_prefix="Neuron_series_demo"
#resource_directory="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/resources"
#ref="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/resources/GENCODE_V39/transcript_to_gene.tsv"
#ref_gene_bed="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/resources/GENCODE_V39/gene.bed"
#ref_transcript_bed="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/resources/GENCODE_V39/transcript.bed.bgz"
#ref_gtf="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/resources/GENCODE_V39/gencode.v39.annotation.gtf.gz"
#SALA_gene_path="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/demo_output_local/sala/gene/table_filtered_gene/Neuron_series_demo"
#===========

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  stop("Usage: Rscript SALA.gene_gtf_annotation <SALA_directory> <out_prefix> <resource_directory> <ref_directory> <SALA_gene_path>\n
  	\n
  	SALA_directory         <required>	path of the folder of SALA transcript annotation output\n
	out_prefix             <required>	output files prefix\n
	resource_directory     <required>	path of the resources folder of SALA\n
	ref_directory          <required>	path of the folder containing the infomation of reference transcriptome\n
	SALA_gene_path         <required>	path of the folder of SALA final gene annotation output")
}

# Assign arguments to variables
SALA_directory <- args[1]
out_prefix <- args[2]
resource_directory <- args[3]
ref <- args[4]
SALA_gene_path <- args[5]

#path
path1 <- paste0(SALA_directory, "/bed/")
path2 <- paste0(SALA_directory, "/log/")
path3 <- paste0(SALA_directory, "/tmp/")
ref_gene_bed <- paste0(ref,"/gene.bed")
ref_transcript_bed <- paste0(ref,"/transcript.bed.bgz")
ref_gtf <- paste0(ref,"/gencode.v39.annotation.gtf.gz")
SALA_gene_info <- paste0(SALA_gene_path,"/log/",out_prefix,".model.info.tsv.gz")
SALA_gene_bed <- paste0(SALA_gene_path,"/bed/",out_prefix,".gene.bed.bgz")

print("Start running...")

#=====
#add gene model
ref_info <- fread(ref, header = FALSE, stringsAsFactors = FALSE)
transcript_info_final <- read.delim(paste0(path2,out_prefix,".table4_filtered.noIP.All_Ref.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
t4gene <- read.delim(SALA_gene_info, header=T, stringsAsFactors = F, check.names = F)
colnames(t4gene)[c(7,8)] <- c("T4_gene_ID","T4_gene_name")
transcript_info_final <- left_join(transcript_info_final,t4gene[,c(1,7,8)], by="model_ID", copy=F)
transcript_info_final$T4_gene_novelty <- "Novel"
transcript_info_final$T4_gene_novelty[which(transcript_info_final$T4_gene_ID %in% ref_info$V2)] <- "Ref"

#=====
#add Ref Transcript Class & Gene Class
colnames(ref_info)[c(3,4)] <- c("Ref_transcriptClass","Ref_geneClass")
transcript_info_final <- left_join(transcript_info_final, unique(ref_info[,c(1,3)]), by=c("model_ID"="V1"),copy=F)
transcript_info_final$Ref_transcriptClass2 <- "NA"
transcript_info_final$Ref_transcriptClass2[which(!is.na(transcript_info_final$Ref_transcriptClass))] <- "others"
transcript_info_final$Ref_transcriptClass2[which(transcript_info_final$Ref_transcriptClass == "lncRNA")] <- "lncRNA"
transcript_info_final$Ref_transcriptClass2[which(transcript_info_final$Ref_transcriptClass == "protein_coding")] <- "protein_coding"
transcript_info_final <- left_join(transcript_info_final, unique(ref_info[,c(2,4)]), by=c("T4_gene_ID"="V2"),copy=F)
transcript_info_final$Ref_geneClass2 <- "NA"
transcript_info_final$Ref_geneClass2[which(!is.na(transcript_info_final$Ref_geneClass))] <- "others"
transcript_info_final$Ref_geneClass2[which(transcript_info_final$Ref_geneClass == "lncRNA")] <- "lncRNA"
transcript_info_final$Ref_geneClass2[which(transcript_info_final$Ref_geneClass == "protein_coding")] <- "protein_coding"

#=====
#add Novel Transcript Class & Gene Class 
if ("CPAT_class" %in% colnames(transcript_info_final)) {
  transcript_info_final$Novel_transcriptClass <- "NA"
  transcript_info_final$Novel_transcriptClass[which(transcript_info_final$transcript_novelty == "Novel" )] <- "others"
  transcript_info_final$Novel_transcriptClass[which(transcript_info_final$transcript_novelty == "Novel" & transcript_info_final$CPAT_class == "non-coding")] <- "ncRNA"
  transcript_info_final$Novel_transcriptClass[which(transcript_info_final$transcript_novelty == "Novel" & transcript_info_final$CPAT_class == "non-coding"  & transcript_info_final$transcript_length > 200)] <- "lncRNA"
  transcript_info_final%>%group_by(Novel_transcriptClass)%>%dplyr::summarise(count=n())
  
  transcript_info_final <- transcript_info_final%>%group_by(T4_gene_ID)%>%dplyr::mutate(overall_T_ratio=full_qry_count/sum(full_qry_count))
  transcript_infoa <- transcript_info_final[which(transcript_info_final$Novel_transcriptClass %in% c("ncRNA","lncRNA") & transcript_info_final$T4_gene_novelty == "Novel"),]
  transcript_infoa <- transcript_infoa%>%group_by(T4_gene_ID)%>%dplyr::mutate(overall_T_ratio2=full_qry_count/sum(full_qry_count))
  transcript_infob <- transcript_infoa%>%group_by(T4_gene_ID)%>%dplyr::summarise(ncRNA_rate=sum(overall_T_ratio), weighted_average_length=sum(overall_T_ratio2*transcript_length))
  transcript_infob$Novel_geneClass <- "others"
  transcript_infob$Novel_geneClass[which(transcript_infob$ncRNA_rate > 0.5)] <- "ncRNA"
  transcript_infob$Novel_geneClass[which(transcript_infob$ncRNA_rate > 0.5 & transcript_infob$weighted_average_length > 200)] <- "lncRNA"
  transcript_info_final <- left_join(transcript_info_final,transcript_infob[,c(1,4)], by="T4_gene_ID", copy=F)
  transcript_info_final$Novel_geneClass[which(is.na(transcript_info_final$Novel_geneClass) & transcript_info_final$T4_gene_novelty == "Novel")] <- "others"
} else {
  message("Skipping step: Novel Transcript Class & Gene Class as CPAT was not performed")
}
#=====
##Ref gene with re-adjusted 3/5'end
gene_model <- read.delim(SALA_gene_bed, header=F, stringsAsFactors = F, check.names = F)
gene_model$V4 <- sapply(strsplit(gene_model$V4,"\\|"),"[",1)
gencode_gene <- read.delim(ref_gene_bed, header=F, stringsAsFactors = F, check.names = F)
gencode_gene <- left_join(gencode_gene[,c(1:6)], gene_model[,c(1:6)], by="V4", copy=F, suffix=c("_ori","_new")) #61533
gencode_gene <- gencode_gene[which(!is.na(gencode_gene$V1_new)),]
gencode_gene$n5_adjust <- "5n_adjust"
gencode_gene$n3_adjust <- "3n_adjust"
gencode_gene$n5_adjust[union(which(gencode_gene$V2_ori==gencode_gene$V2_new & gencode_gene$V6_ori=="+"), which(gencode_gene$V3_ori==gencode_gene$V3_new & gencode_gene$V6_ori=="-"))] <- NA
gencode_gene$n3_adjust[union(which(gencode_gene$V2_ori==gencode_gene$V2_new & gencode_gene$V6_ori=="-"), which(gencode_gene$V3_ori==gencode_gene$V3_new & gencode_gene$V6_ori=="+"))] <- NA
gencode_gene2 <- reshape2::melt(gencode_gene[,c(4,12,13)], id=1)
gencode_gene2 <- gencode_gene2[which(!is.na(gencode_gene2$value)),]
gencode_gene2 <- gencode_gene2%>%group_by(V4)%>%dplyr::summarise(Ref_gene_adjust=paste(value,collapse=" & "))
gencode_gene <- left_join(gencode_gene,gencode_gene2, by="V4", copy=F)
gencode_gene$Ref_gene_adjust[is.na(gencode_gene$Ref_gene_adjust)] <- "No"
transcript_info_final <- left_join(transcript_info_final,gencode_gene[,c(4,14)], by=c("T4_gene_ID" = "V4"), copy=F)

#=====
##Ref transcript with re-adjusted 3/5'end
transcript_info1aa <- transcript_info_final[which(transcript_info_final$ref_source == "fulllength_ref"),]
transcript_model <- read.delim(paste0(path1,out_prefix,".model.bed.bgz"), header=F, stringsAsFactors = F, check.names = F)
transcript_model_genecode <- transcript_model[which(transcript_model$V4 %in% transcript_info1aa$model_ID),]
gencode_tran <- read.delim(ref_transcript_bed, header=F, stringsAsFactors = F, check.names = F)
gencode_tran1 <- gencode_tran[which(gencode_tran$V4 %in% transcript_info1aa$model_ID),]
gencode_tran1 <- left_join(gencode_tran1[,c(1:6)], transcript_model_genecode[,c(1:6)], by="V4", copy=F, suffix=c("_ori","_new"))
gencode_tran1$n5_adjust <- "5n_adjust"
gencode_tran1$n3_adjust <- "3n_adjust"
gencode_tran1$n5_adjust[union(which(gencode_tran1$V2_ori==gencode_tran1$V2_new & gencode_tran1$V6_ori=="+"), which(gencode_tran1$V3_ori==gencode_tran1$V3_new & gencode_tran1$V6_ori=="-"))] <- NA
gencode_tran1$n3_adjust[union(which(gencode_tran1$V2_ori==gencode_tran1$V2_new & gencode_tran1$V6_ori=="-"), which(gencode_tran1$V3_ori==gencode_tran1$V3_new & gencode_tran1$V6_ori=="+"))] <- NA
gencode_tran2 <- reshape2::melt(gencode_tran1[,c(4,12,13)], id=1)
gencode_tran2 <- gencode_tran2[which(!is.na(gencode_tran2$value)),]
gencode_tran2 <- gencode_tran2%>%group_by(V4)%>%dplyr::summarise(Ref_transcript_adjust=paste(value,collapse=" & "))
gencode_tran1 <- left_join(gencode_tran1,gencode_tran2, by="V4", copy=F)
gencode_tran1$Ref_transcript_adjust[is.na(gencode_tran1$Ref_transcript_adjust)] <- "No"
transcript_info_final <- left_join(transcript_info_final,gencode_tran1[,c(4,14)], by=c("model_ID" = "V4"), copy=F)
write.table(transcript_info_final, gzfile(paste0(path2,out_prefix,".table4_filtered.noIP.All_Ref_updated.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
rm(gene_model,gencode_gene,transcript_model,gencode_tran)
print(paste0(out_prefix,".table4_filtered.noIP.All_Ref_updated.tsv.gz exported to ",path2))
print("start building gtf")

###make gtf from bed12, bed6 exon, gene bed, gene info and GENCODE gtf

gencode.gtf <- fread(ref_gtf, header=F, stringsAsFactors = F)
table0.bed12 <- read.delim(paste0(path1,out_prefix,".model.bed.bgz"),header=F, stringsAsFactors = F, check.names = F)
table0.bed6 <- read.delim(paste0(path1,out_prefix,".model.bed6.bed.gz"),header=F, stringsAsFactors = F, check.names = F)
table4.genebed <- read.delim(SALA_gene_bed, header=F, stringsAsFactors = F, check.names = F)

transcript_info1aa <- transcript_info_final[which(transcript_info_final$Ref_transcript_adjust != "No" & !is.na(transcript_info_final$Ref_transcript_adjust)),]
transcript_info1bb <- transcript_info_final[which(transcript_info_final$Ref_gene_adjust != "No" & !is.na(transcript_info_final$Ref_gene_adjust)),]
gencode.gtf.g <- gencode.gtf[which(gencode.gtf$V3 == "gene"),]
gencode.gtf.t <- gencode.gtf[which(gencode.gtf$V3 == "transcript"),]
gencode.gtf.e <- gencode.gtf[which(gencode.gtf$V3 == "exon"),]

gencode.gtf.g$gene_ID <- str_match(gencode.gtf.g$V9, 'gene_id "([^"]+)"')[,2]
gencode.gtf.g$gene_type <- str_match(gencode.gtf.g$V9, 'gene_type "([^"]+)"')[,2]
gencode.gtf.g$gene_name <- str_match(gencode.gtf.g$V9, 'gene_name "([^"]+)"')[,2]
gencode.gtf.ga <- gencode.gtf.g[-which(gencode.gtf.g$gene_ID %in% unique(transcript_info1bb$T4_gene_ID)),] 
gencode.gtf.ga <- gencode.gtf.ga[,c(1:8,10:12)] #component1
gencode.gtf.ga$V2 <- "Reference"
gencode.gtf.ga$gene_novelty <- "Reference"

table4.genebed$gene_ID <- sapply(strsplit(table4.genebed$V4,"\\|"),"[",1)
table4.genebed$label <- "SALA"
table4.genebed$type <- "gene"
table4.genebed <- table4.genebed[,c(1,14,15,2,3,5,6,5,13)]
table4.genebed$V2 <- table4.genebed$V2+1
table4.genebed$V5 <- "."
table4.genebed$V5.1 <- "."
table4a.genebed <- table4.genebed[-which(table4.genebed$gene_ID %in% unique(ref_info$V2)),]
if ("Novel_geneClass" %in% colnames(transcript_info_final)){
  table4a.genebed <- left_join(table4a.genebed, unique(transcript_info_final[,c("T4_gene_ID","Novel_geneClass")]), by=c("gene_ID"="T4_gene_ID"),copy=F)
} else {
  table4a.genebed$Novel_geneClass <- "novel_gene"}
table4a.genebed$gene_name <- sapply(strsplit(table4a.genebed$gene_ID,"\\."),"[",1)
table4a.genebed$gene_novelty <- "novel"
table4b.genebed <- table4.genebed[which(table4.genebed$gene_ID %in% unique(transcript_info1bb$T4_gene_ID)),]
table4b.genebed <- left_join(table4b.genebed,gencode.gtf.g[,c(10:12)],by="gene_ID",copy=F)
table4b.genebed$gene_novelty <- "Reference_updated"
table4b.genebed$label <- "Reference"
colnames(table4a.genebed) <- colnames(gencode.gtf.ga)
colnames(table4b.genebed) <- colnames(gencode.gtf.ga)
table4.genebed <- rbind(table4a.genebed,table4b.genebed,gencode.gtf.ga)
table4.genebed$V9 <- paste0("gene_id \"",table4.genebed$gene_ID,"\"; gene_type \"",table4.genebed$gene_type,"\"; gene_name \"",table4.genebed$gene_name,"\"; gene_novelty \"",table4.genebed$gene_novelty,"\";")
#

gencode.gtf.t$transcript_ID <- str_match(gencode.gtf.t$V9, 'transcript_id "([^"]+)"')[,2]
gencode.gtf.t$gene_ID <- str_match(gencode.gtf.t$V9, 'gene_id "([^"]+)"')[,2]
gencode.gtf.t$transcript_type <- str_match(gencode.gtf.t$V9, 'transcript_type "([^"]+)"')[,2]
gencode.gtf.t$transcript_name <- str_match(gencode.gtf.t$V9, 'transcript_name "([^"]+)"')[,2]
gencode.gtf.ta <- gencode.gtf.t[-which(gencode.gtf.t$transcript_ID %in% transcript_info1aa$model_ID),] 
gencode.gtf.ta <- gencode.gtf.ta[,c(1:8,10:13)] #component1
gencode.gtf.ta$V2 <- "Reference"
gencode.gtf.ta$transcript_novelty <- "Reference"

table4.bed12 <- table0.bed12[which(table0.bed12$V4 %in% transcript_info_final$model_ID),]
table4.bed12$label <- "SALA"
table4.bed12$type <- "transcript"
table4.bed12 <- table4.bed12[,c(1,13,14,2,3,5,6,5,4)]
table4.bed12$V2 <- table4.bed12$V2+1
table4.bed12$V5 <- "."
table4.bed12$V5.1 <- "."
table4a.bed12 <- table4.bed12[-which(table4.bed12$V4 %in% ref_info$V1),]
table4a.bed12 <- left_join(table4a.bed12, transcript_info_final[,c("model_ID","T4_gene_ID")], by=c("V4"="model_ID"),copy=F)
if ("Novel_transcriptClass" %in% colnames(transcript_info_final)){
  table4a.bed12 <- left_join(table4a.bed12, transcript_info_final[,c("model_ID","Novel_transcriptClass")], by=c("V4"="model_ID"),copy=F)
} else {
  table4a.bed12$Novel_transcript_class <- "novel_transcript"}
table4a.bed12$transcript_name <- sapply(strsplit(table4a.bed12$V4,"\\."),"[",1)
table4a.bed12$transcript_novelty <- "novel"
table4b.bed12 <- table4.bed12[which(table4.bed12$V4 %in% transcript_info1aa$model_ID),]
table4b.bed12 <- left_join(table4b.bed12,gencode.gtf.t[,c(10:13)],by=c("V4"="transcript_ID"),copy=F) #take geneID, transcript_type, transcript_name
table4b.bed12$transcript_novelty <- "Reference_updated"
table4b.bed12$label <- "Reference"
colnames(table4a.bed12) <- colnames(gencode.gtf.ta)
colnames(table4b.bed12) <- colnames(gencode.gtf.ta)
table4.bed12 <- rbind(table4a.bed12,table4b.bed12,gencode.gtf.ta)
table4.bed12 <- left_join(table4.bed12,table4.genebed[,c(9:12)], by="gene_ID", copy=F) #get back gene info from last table
table4.bed12$V9 <- paste0("gene_id \"",table4.bed12$gene_ID,"\"; transcript_id \"",table4.bed12$transcript_ID,"\"; gene_type \"",table4.bed12$gene_type,"\"; gene_name \"",table4.bed12$gene_name,"\"; transcript_type \"",table4.bed12$transcript_type,"\"; transcript_name \"",table4.bed12$transcript_name,"\"; gene_novelty \"",table4.bed12$gene_novelty,"\"; transcript_novelty \"",table4.bed12$transcript_novelty,"\";")
#

gencode.gtf.e$transcript_ID <- str_match(gencode.gtf.e$V9, 'transcript_id "([^"]+)"')[,2]
gencode.gtf.e$exon_number <- str_match(gencode.gtf.e$V9, 'exon_number ([^"]+);')[,2]
gencode.gtf.e$exon_id <- str_match(gencode.gtf.e$V9, 'exon_id "([^"]+)"')[,2]
gencode.gtf.ea <- gencode.gtf.e[-which(gencode.gtf.e$transcript_ID %in% transcript_info1aa$model_ID),]  
gencode.gtf.ea <- gencode.gtf.ea[,c(1:8,10:11)] #component1
gencode.gtf.ea$exon_number <- as.character(gencode.gtf.ea$exon_number)

table4.bed6 <- table0.bed6[which(table0.bed6$V4 %in% transcript_info_final$model_ID),]
table4.bed6$label <- "SALA"
table4.bed6$type <- "exon"
table4.bed6 <- table4.bed6[,c(1,7,8,2,3,5,6,5,4)]
table4.bed6$V2 <- table4.bed6$V2+1
table4.bed6$V5 <- "."
table4.bed6$V5.1 <- "."
table4.bed6p <- table4.bed6[which(table4.bed6$V6 == "+"),]%>%group_by(V4)%>%dplyr::arrange(V2)%>%dplyr::mutate(exon_number = 1: n())
table4.bed6n <- table4.bed6[which(table4.bed6$V6 == "-"),]%>%group_by(V4)%>%dplyr::arrange(desc(V2))%>%dplyr::mutate(exon_number = 1: n())
table4.bed6 <- rbind(table4.bed6p,table4.bed6n)
table4.bed6$exon_number <- as.character(table4.bed6$exon_number)

table4a.bed6 <- table4.bed6[-which(table4.bed6$V4 %in% ref_info$V1),]
table4b.bed6 <- table4.bed6[which(table4.bed6$V4 %in% transcript_info1aa$model_ID),]
colnames(table4a.bed6) <- colnames(gencode.gtf.ea)
colnames(table4b.bed6) <- colnames(gencode.gtf.ea)
table4.bed6 <- rbind(table4a.bed6,table4b.bed6,gencode.gtf.ea)

table4.bed6 <- left_join(table4.bed6, table4.bed12[,c(9:12,14:15)], by="transcript_ID", copy=F)
table4.bed6$V9 <- paste0("gene_id \"",table4.bed6$gene_ID,"\"; transcript_id \"",table4.bed6$transcript_ID,"\"; gene_type \"",table4.bed6$gene_type,"\"; gene_name \"",table4.bed6$gene_name,"\"; transcript_type \"",table4.bed6$transcript_type,"\"; transcript_name \"",table4.bed6$transcript_name,"\"; exon_number \"",table4.bed6$exon_number,"\";")
#

options(scipen=999)
final <- rbind(table4.genebed[,c(1:8,13)],table4.bed12[,c(1:8,17)], table4.bed6[,c(1:8,16)])
final <- final[order(final$V1,final$V4),]
write.table(final,gzfile(paste0(path2,"table4.All_Ref.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

print(paste0("table4.All_Ref.gtf.gz exported to ",path2))

#=========
#validate the gtf
#gzip -d -c table4.All_Ref.gtf.gz | bedparse gtf2bed | gzip > table4.All_Ref.bed12.bed.gz
#================================================================================

#============================================================
#detectable alone
transcript_info_final1 <- transcript_info_final[which(transcript_info_final$ref_source != "non_detectable_ref"),]
t4transcript <- unique(transcript_info_final1$model_ID)
t4gene <- unique(transcript_info_final1$T4_gene_ID)

#===split the gtf into transcript_info4 and revise gene range for all (all contain all ENSG and ENST)===
gtf <- final
gtf.gene <- final[which(final$V3=="gene"),]
gtf.nogene <- final[which(final$V3!="gene"),]
gtf.gene$gene_ID <- str_match(gtf.gene$V9, 'gene_id "([^"]+)"')[,2]
gtf.nogene$transcript_ID <- str_match(gtf.nogene$V9, 'transcript_id "([^"]+)"')[,2]
gtf.t <- gtf.nogene[which(gtf.nogene$V3=="transcript"),]
gtf.gene4 <- gtf.gene[which(gtf.gene$gene_ID %in% t4gene),]
gtf.nogene4 <- gtf.nogene[which(gtf.nogene$transcript_ID %in% t4transcript),]

#revise gene region start end
gtf.t$gene_ID <- str_match(gtf.t$V9, 'gene_id "([^"]+)"')[,2]
gtf.t <- gtf.t[,c(4,5,10,11)]
gtf.t4 <- gtf.t[which(gtf.t$transcript_ID %in% t4transcript),]
gtf.g4 <- gtf.t4%>%group_by(gene_ID)%>%dplyr::summarise(gene_start=min(V4), gene_end=max(V5))

gtf.gene4 <- left_join(gtf.gene4,gtf.g4,by="gene_ID", copy=F)
gtf.gene4 <- gtf.gene4[,c(1:3,11,12,6:9)]
colnames(gtf.gene4) <- colnames(gtf.nogene4)[c(1:9)]

gtf4 <- rbind(gtf.gene4,gtf.nogene4[,c(1:9)])
write.table(gtf4[order(gtf4$V1,gtf4$V4),],gzfile(paste0(path2,"table4.Detected_Ref.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

print(paste0("table4.Detected_Ref.gtf.gz exported to ",path2))
#==========================================================================================
#validate the gtf
#gzip -d -c table4.Detected_Ref.gtf.gz | bedparse gtf2bed | gzip > table4.Detected_Ref.bed12.bed.gz
#==========================================================================================







