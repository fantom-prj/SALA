# List of required packages
required_packages <- c(
  "data.table", "magrittr" ,"stringr",  "Rsamtools"
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
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Rsamtools))

#ref_gtf <- "/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/resources/GENCODE_V39/gencode.v39.annotation.gtf"
#output_path <- "/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/resources/GENCODE_V39"
#fasta_file="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript SALA.gtf_extract.R <ref_gtf> <output_path> <fasta_file>")
}

# Assign arguments to variables
ref_gtf <- args[1]
output_path <- args[2]
fasta_file <- args[3]
fasta_path <- dirname(fasta_file)

#=====
gencode.gtf <- fread(ref_gtf, header=F, stringsAsFactors = F)
gencode.gtf.t <- gencode.gtf[which(gencode.gtf$V3 == "transcript"),]
gencode.gtf.t$transcript_ID <- str_match(gencode.gtf.t$V9, 'transcript_id "([^"]+)"')[,2]
gencode.gtf.t$gene_ID <- str_match(gencode.gtf.t$V9, 'gene_id "([^"]+)"')[,2]
gencode.gtf.t$transcript_type <- str_match(gencode.gtf.t$V9, 'transcript_type "([^"]+)"')[,2]
gencode.gtf.t$gene_type <- str_match(gencode.gtf.t$V9, 'gene_type "([^"]+)"')[,2]
gencode.gtf.t$gene_name <- str_match(gencode.gtf.t$V9, 'gene_name "([^"]+)"')[,2]
gencode.gtf.t1 <- gencode.gtf.t[,c("transcript_ID","gene_ID","transcript_type","gene_type","gene_name")]
write.table(gencode.gtf.t1, paste0(output_path,"/transcript_to_gene.tsv"), col.names=F, row.names=F, sep="\t", quote=F)
print(paste0(output_path, "/transcript_to_gene.tsv is generated"))

#===========
#prepare n3 and n5 bed6 bed file for reference transcriptome
gencode.gtf.t_3n=gencode.gtf.t
gencode.gtf.t_5n=gencode.gtf.t

gencode.gtf.t_5n$V5[which(gencode.gtf.t_5n$V7=="+")] <- gencode.gtf.t_5n$V4[which(gencode.gtf.t_5n$V7=="+")]
gencode.gtf.t_5n$V4[which(gencode.gtf.t_5n$V7=="-")] <- gencode.gtf.t_5n$V5[which(gencode.gtf.t_5n$V7=="-")]
gencode.gtf.t_5n$V4=gencode.gtf.t_5n$V4-1
gencode.gtf.t_5n$label <- paste0(gencode.gtf.t_5n$V1,"_",gencode.gtf.t_5n$V4,"_",gencode.gtf.t_5n$V5,"_",gencode.gtf.t_5n$V7)
gencode.gtf.t_5n=unique(gencode.gtf.t_5n[order(gencode.gtf.t_5n$V1,gencode.gtf.t_5n$V4),c(1,4,5,15,6,7)])
write.table(gencode.gtf.t_5n,gzfile(paste0(output_path,"/n5.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

gencode.gtf.t_3n$V5[which(gencode.gtf.t_3n$V7=="-")] <- gencode.gtf.t_3n$V4[which(gencode.gtf.t_3n$V7=="-")]
gencode.gtf.t_3n$V4[which(gencode.gtf.t_3n$V7=="+")] <- gencode.gtf.t_3n$V5[which(gencode.gtf.t_3n$V7=="+")]
gencode.gtf.t_3n$V4=gencode.gtf.t_3n$V4-1
gencode.gtf.t_3n$label <- paste0(gencode.gtf.t_3n$V1,"_",gencode.gtf.t_3n$V4,"_",gencode.gtf.t_3n$V5,"_",gencode.gtf.t_3n$V7)
gencode.gtf.t_3n=unique(gencode.gtf.t_3n[order(gencode.gtf.t_3n$V1,gencode.gtf.t_3n$V4),c(1,4,5,15,6,7)])
write.table(gencode.gtf.t_3n,gzfile(paste0(output_path,"/n3.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
rm(gencode.gtf.t_5n, gencode.gtf.t_3n)

#===========
#prepare gene bed
gencode.gtf.g <- gencode.gtf[which(gencode.gtf$V3 == "gene"),]
gencode.gtf.g$gene_ID <- str_match(gencode.gtf.g$V9, 'gene_id "([^"]+)"')[,2]
gencode.gtf.g$V4 <- gencode.gtf.g$V4-1
gencode.gtf.g <- gencode.gtf.g[order(gencode.gtf.g$V1,gencode.gtf.g$V4),c(1,4,5,10,6,7)]
write.table(gencode.gtf.g,gzfile(paste0(output_path,"/gene.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
print(paste0("n5.bed.gz, n3.bed.gz & gene.bed.gz are prepared in ",output_path))

#==========
#prepare chrom.sizes.tsv
indexFa(fasta_file)
genome <- FaFile(fasta_file)
open(genome)
chrom_sizes <- seqlengths(genome)
chrom_table <- data.frame(Chromosome = names(chrom_sizes), Length = as.numeric(chrom_sizes))
write.table(chrom_table, paste0(fasta_path,"/chrom.sizes.tsv"), col.names=F, row.names=F, sep="\t", quote=F)
print(paste0("chrom.sizes.tsv is prepared in ",fasta_path))

