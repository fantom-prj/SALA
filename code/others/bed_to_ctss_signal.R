library(GenomicRanges)
library(data.table)
library(dplyr)
library(magrittr)
library(knitr)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript bed_to_ctss_singal.R <bed_file_list> <output_file> <resource_directory>")
}

# Assign arguments to variables

bed_file_path <- args[1]
output_file_path <- args[2]
resource_directory <- args[3]
bgz_path <- gsub("gz","bgz",output_file_path)

# Read the input files
path_file <- read.delim(bed_file_path, header = FALSE, stringsAsFactors = FALSE)

# Process each file in path_file
data <- read.delim(path_file$V2[1], header = FALSE, stringsAsFactors = FALSE)
data$V3[which(data$V6 == "+")] <- data$V2[which(data$V6 == "+")] + 1
data$V2[which(data$V6 == "-")] <- data$V3[which(data$V6 == "-")] - 1
dataa=data%>%group_by(V1,V2,V3,V6)%>%dplyr::summarise(count=n(), .groups = "drop_last")
for (i in 2:nrow(path_file)) {
  data1 <- read.delim(path_file$V2[i], header = FALSE, stringsAsFactors = FALSE)
  data1$V3[which(data1$V6 == "+")] <- data1$V2[which(data1$V6 == "+")] + 1
  data1$V2[which(data1$V6 == "-")] <- data1$V3[which(data1$V6 == "-")] - 1
  data1a=data1%>%group_by(V1,V2,V3,V6)%>%dplyr::summarise(count=n(), .groups = "drop_last")
  dataa=rbind(dataa,data1a)}

data2=dataa%>%group_by(V1,V2,V3,V6)%>%dplyr::summarise(count=sum(count), .groups = "drop_last")

# Write the output file
write.table(data2[order(data2$V1,data2$V2), c("V1", "V2", "V3", "count", "count", "V6")], 
            gzfile(output_file_path), 
            col.names = FALSE, 
            row.names = FALSE, 
            sep = "\t", 
            quote = FALSE)

#index bed file
tabix_bin <- paste0(resource_directory,"/bin/tabix/tabix")
bgzip_bin <- paste0(resource_directory,"/bin/bgzip/bgzip")
system(paste0("zcat ",output_file_path, "| ",bgzip_bin, " > ",bgz_path))
system(paste0(tabix_bin," -p bed ",bgz_path))


