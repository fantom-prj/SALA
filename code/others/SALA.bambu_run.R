# List of required packages
required_packages <- "bambu"


# Function to check and install missing packages
install_if_missing <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse=", "))
    devtools::install_github("GoekeLab/bambu", ref = "test_split_read_classes")
  }
}

# Ensure BiocManager is installed for Bioconductor packages
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org", dependencies = TRUE)
}

# Install missing packages
install_if_missing(required_packages)

# load packages
suppressPackageStartupMessages(library(bambu))

#bam_path <- "/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/data/SCAFE.step1.bam.list.txt"
#fa.file <- "/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/resources/chr17_chr18.fasta"
#SALA_gtf <- "/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/demo_output_local/sala/transcript/Neuron_series_demo/log/table4.Detected_Ref.gtf.gz"
#results_dir <- "/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/demo_output_local/sala/transcript/Neuron_series_demo/bambu"


# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript SALA.bambu_run.R <bam_path> <genome_fasta> <SALA_gtf> <bambu_results_dir>")
}

# Assign arguments to variables
bam_path <- args[1]
fa.file <- args[2]
SALA_gtf <- args[3]
results_dir <- args[4]

#=====
#run bambu
print("starting running bambu...")
bam <- read.delim(bam_path, header = F, stringsAsFactors = F, check.names = F)
test.bam <- bam$V3
annotations <- prepareAnnotations(SALA_gtf)
ID_link=unique(data.frame(cbind(mcols(annotations)$GENEID, bambu:::assignGeneIds(annotations, GRangesList())$GENEID)))
ID_link=ID_link%>%group_by(X2)%>%dplyr::summarise(geneID=paste(X1,collapse=";"))
mcols(annotations)$GENEID <- bambu:::assignGeneIds(annotations, GRangesList())$GENEID
dir.create(paste0(results_dir,"/rcOut"), recursive=TRUE)
se <- bambu(reads =  test.bam, 
                    annotations = annotations, 
                    genome = fa.file, 
                    discovery = FALSE, 
                    opt.discovery = list(min.exonDistance = 0), 
                    rcOutDir = paste0(results_dir,"/rcOut"),returnDistTable=TRUE)
saveRDS(se, paste0(results_dir,"/se.rds"))
writeBambuOutput(se, results_dir)
print(paste0("finish running bambu. results located in ",results_dir))

#=====
gene_count <- read.delim(paste0(results_dir,"/counts_gene.txt"), header=T)
transcript_count <- read.delim(paste0(results_dir,"/counts_transcript.txt"), header=T)
transcript_CMP <- read.delim(paste0(results_dir,"/CPM_transcript.txt"), header=T)
gene_count <- left_join(gene_count, ID_link, by=c("GENEID"="X2"), copy=F)
transcript_count <- left_join(transcript_count, ID_link, by=c("GENEID"="X2"), copy=F)
transcript_CMP <- left_join(transcript_CMP, ID_link, by=c("GENEID"="X2"), copy=F)
write.table(gene_count, paste0(results_dir,"/counts_gene.txt"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(transcript_count, paste0(results_dir,"/counts_transcript.txt"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(transcript_CMP, paste0(results_dir,"/CPM_transcript.txt"), col.names=T, row.names=F, sep="\t", quote=F)







