required_packages <- c(
  "dplyr", "magrittr"
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
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))

#===============================================================================
tag="CAT.Neuron_series_demo"
transcript_prefix="ONTT"
gene_prefix="ONTG"
compare_path="~/tool2026/SALA/code/SALA_compare/CAT"

#===========
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript extract_identical.R <out_prefix> <transcript_prefix> <gene_prefix> <compare_path> \n\n",
       "out_prefix                  <required>	output files prefix\n",
       "transcript_prefix           <required>	prefix of novel transcript models\n",
       "gene_prefix                 <required>	prefix of gene transcript models\n",
       "compare_path                <required>	path of the comparison")}

# Assign arguments to variables
tag <- args[1]
transcript_prefix <- args[2]
gene_prefix <- args[3]
compare_path <- args[4]

#=====
compare_path_log=paste0(compare_path,"/transcript/",tag,"/log/")
compare_path_log_gene=paste0(compare_path,"/gene/",tag,"/log/")


#===============================================================================
# transcript level
read_info=read.delim(paste0(compare_path_log,tag,".trnscpt.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
read_info$group1=substr(read_info$trnscpt_ID, start = 1, stop = 4)
read_info$group2=substr(read_info$model_ID_str, start = 1, stop = 4)

read_info1=read_info[-which(read_info$group1 %in% c("ENST",transcript_prefix)),]
read_info2=read_info1[which(read_info1$group2 == transcript_prefix),]%>%group_by(model_ID_str)%>%dplyr::summarise(annot_transcriptID=paste(trnscpt_ID,collapse=";"))
write.table(read_info2,gzfile(paste0(compare_path,"/match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
sala_same=read_info2$model_ID_str
print(paste0("match_transcript.tsv.gz is exported to ",compare_path))

#===============================================================================
# gene level
gene.info=read.delim(paste0(compare_path_log_gene,tag,".model.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
gene.info$group1=substr(gene.info$model_ID, start = 1, stop = 4)
gene.info$group2=substr(gene.info$gene_ID, start = 1, stop = 4)

gene.info1=gene.info[which(gene.info$group1 == "NEWT"),]
ONTGout1=unique(gene.info$gene_ID[which(gene.info$group2 == gene_prefix & gene.info$group1 == "NEWT")])
ONTGout2=unique(gene.info$gene_ID[which(gene.info$model_ID %in% sala_same)])
ONTGout=union(ONTGout1,ONTGout2)

#show which ONTG contain transcript of other database
gene.info3=gene.info[which(gene.info$gene_ID %in% ONTGout),]
read_info3=read_info[which(read_info$model_ID_str %in% gene.info3$model_ID),c(1,7,15,16)]
read_info3=read_info3[-which(read_info3$group1 %in% c("ENST",transcript_prefix)),]
read_info_collap=read_info3%>%group_by(model_ID_str)%>%dplyr::summarise(annot_transcriptID=paste(trnscpt_ID, collapse=";"))
gene.info3=left_join(gene.info3,read_info_collap, by=c("model_ID"="model_ID_str"), copy=F)
gene_list=gene.info3[which(!is.na(gene.info3$annot_transcriptID)),]%>%group_by(gene_ID)%>%dplyr::summarise(annot_transcriptID=paste(annot_transcriptID, collapse=";"))
write.table(gene_list,gzfile(paste0(compare_path,"/Gene_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
print(paste0("Gene_match_transcript.tsv.gz is exported to ",compare_path))

