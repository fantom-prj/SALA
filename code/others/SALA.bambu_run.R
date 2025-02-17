# List of required packages
required_packages <- c(
  "bambu"
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
suppressPackageStartupMessages(library(bambu))

bam_path <- "/Users/yip/Library/CloudStorage/OneDrive-Personal/Documents/my_gid/SALA/data/SCAFE.step1.bam.list.txt"
fa.file <- "/Users/yip/Library/CloudStorage/OneDrive-Personal/Documents/my_gid/SALA/resources/for_demo/chr20_subset.fasta"
SALA_gtf <- "/Users/yip/Library/CloudStorage/OneDrive-Personal/Documents/my_gid/SALA/demo_output_local/sala/transcript/Neuron_series_demo/log/table4.Detected_Ref.gtf.gz"
results_dir <- "/Users/yip/Library/CloudStorage/OneDrive-Personal/Documents/my_gid/SALA/demo_output_local/sala/transcript/Neuron_series_demo/bambu"

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
