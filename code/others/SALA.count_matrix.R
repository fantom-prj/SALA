# List of required packages
required_packages <- c( 
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
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyr))

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript script_name.R <SALA_directory> <output_directory> <ref.transcriptome_path>")
}

# Assign arguments to variables
SALA_directory <- args[1]
output_directory <- args[2]
ref <- args[3]

SALA_directory="/Users/yip/Library/CloudStorage/OneDrive-Personal/Documents/my_gid/SALA/demo_output_local/sala/transcript/Neuron_series_demo"
output_directory="/Users/yip/Library/CloudStorage/OneDrive-Personal/Documents/my_gid/SALA/demo_output_local/sala/transcript/Neuron_series_demo/log"
ref="/Users/yip/Library/CloudStorage/OneDrive-Personal/Documents/my_gid/SALA/resources/GENCODE_V39/transcript_to_gene.tsv"
#SALA_directory="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/demo_output_local/sala/transcript/Neuron_series_demo"
#output_directory="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/demo_output_local/sala/transcript/Neuron_series_demo/log"
#ref="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/resources/GENCODE_V39/transcript_to_gene.tsv"

# Set paths
path2 <- paste0(SALA_directory, "/log/")
path3 <- paste0(SALA_directory, "/tmp/")

###### Read log table
info <- list.files(path = path2, pattern = "model.info.tsv.gz")
transcript_info <- fread(paste0(path2, info[1]), header = TRUE, stringsAsFactors = FALSE)
ref_info <- fread(ref, header = FALSE, stringsAsFactors = FALSE)

files <- list.files(path = path3, pattern = "transcript.info.tsv.gz", recursive = TRUE)
files.names <- sapply(strsplit(files, "\\/"), "[", 1)

for (i in 1:length(files)){
  read_info = tryCatch({
    fread(paste0(path3, files[i]), header = F, stringsAsFactors = F, select = c(1, 4, 7))
  }, error = function(e) NULL)  # If fread fails, return NULL
  
  # Check if the file is empty using nrow()
  if (is.null(read_info) || nrow(read_info) == 0) {
    print(paste0("Skipping empty or invalid file: ", files[i]))
    next  # Skip to the next file in the loop
  }
  
  colnames(read_info)=c("trnscpt_ID","set_ID","model_ID_str")
  read_info=read_info[-which(read_info$trnscpt_ID %in% ref_info$V1),]
  read_info$rep = sapply(strsplit(read_info$trnscpt_ID,"_"),"[",1)
  
  read_info2=read_info[which(read_info$set_ID %in% unique(transcript_info$full_set_ID)),c(2,4)]
  read_info2=read_info2%>%group_by(rep,set_ID)%>%dplyr::summarise(count=n(), .groups = "drop_last")
  read_info3=spread(read_info2, key=1, value=3)
  read_info3=right_join(transcript_info[,c(1,4)],read_info3, by=c("full_set_ID"="set_ID"),copy=F)
  read_info3=read_info3[,-2]
  read_info3[is.na(read_info3)]=0
  write.table(read_info3, gzfile(paste0(path3, files.names[i],"_full_length_support_matrix.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
  
  read_info1a=read_info[,c(4,3)]%>%group_by(rep,model_ID_str)%>%dplyr::summarise(count=n(), .groups = "drop_last")
  read_info1a=separate_rows(read_info1a, model_ID_str, sep = ";", convert = FALSE)
  read_info1a=read_info1a%>%group_by(rep,model_ID_str)%>%dplyr::summarise(count=sum(count), .groups = "drop_last")
  read_info2a=spread(read_info1a, key=1, value=3)
  read_info2a=right_join(transcript_info[,c(1,4)],read_info2a, by=c("model_ID"="model_ID_str"), copy=F)
  read_info2a=read_info2a[,-2]
  read_info2a[is.na(read_info2a)]=0
  write.table(read_info2a, gzfile(paste0(path3, files.names[i],"_partial_length_support_matrix.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
  print(paste0("Processed: ",files.names[i]))}

# Aggregate full-length support matrices
files <- list.files(path = path3, pattern = "full_length_support_matrix.tsv.gz")
matr <- read.delim(paste0(path3,files[1]), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
for (i in seq(2, length(files))) {
  matr1 <- read.delim(paste0(path3,files[i]), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  matr <- rbind(matr, matr1)
}
write.table(matr, gzfile(paste0(output_directory, "/full_length_support_matrix.tsv.gz")), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Similarly for partial-length support matrices
files <- list.files(path = path3, pattern = "partial_length_support_matrix.tsv.gz")
matr <- read.delim(paste0(path3,files[1]), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
for (i in seq(2, length(files))) {
  matr1 <- read.delim(paste0(path3,files[i]), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  matr <- rbind(matr, matr1)
}
write.table(matr, gzfile(paste0(output_directory, "/partial_length_support_matrix.tsv.gz")), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)



