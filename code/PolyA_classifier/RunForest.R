#!/usr/bin/env Rscript

require(GenomicRanges)
require(GenomicFeatures)
require(dplyr) 
require(tidyr)
require(ggplot2)
require(Biostrings)
require(rtracklayer)
require(readr)
require(stringr)
require(caTools)
require(ranger)
require(patchwork)
require(memes)
require(gridExtra)
require(plotROC)
source("code/rf_helpers.R")

# ---------------------------------------------------------------------------- #


library(optparse)

option_list <- list(
  make_option(c("-g", "--genome"), type = "character", help = "Genome fasta file"),
  make_option(c("-r", "--reference_set"), type = "character", help = "Reference dataset BED file"),
  make_option(c("-t", "--test_set"), type = "character", help = "Test dataset BED file"),
  make_option(c("-e", "--output_gtf"), type = "character", help = "Output GTF with ends"),
  make_option(c("-m", "--output_model"), type = "character", help = "Output RDS file with trained model"),
  make_option(c("-f", "--output_figure"), type = "character", help = "Output directory for figure PDFs"),
  make_option(c("-k", "--threads"), type = "integer", default = 1, help = "Number of threads used by ranger()"),
  make_option(c("-s", "--seed"), type = "integer", default = 123, help = "Seed [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

genome_fasta   <- opt$genome
test_ends      <- opt$test_set
reference_ends <- opt$reference_set
output_gtf     <- opt$output_gtf
output_model   <- opt$output_model
outfig         <- opt$output_figure

# ---------------------------------------------------------------------------- #

dir.create(dirname(output_gtf), recursive = TRUE)
dir.create(dirname(output_model), recursive = TRUE)
dir.create(outfig, recursive = TRUE)

# ---------------------------------------------------------------------------- #
# Load reference genome

reference_genome        <- readDNAStringSet(genome_fasta)
names(reference_genome) <- gsub(" .+", "", names(reference_genome))

reference_chr <- names(reference_genome)

# ---------------------------------------------------------------------------- #
# Load the test dataset transcript coordinates and extract 3' ends

tx_ends_gr <- import.bed(test_ends)
tx_ends_gr <- tx_ends_gr[seqnames(tx_ends_gr) %in% reference_chr, ]
seqlevels(tx_ends_gr)  <- names(reference_genome)
seqlengths(tx_ends_gr) <- width(reference_genome[seqlevels(tx_ends_gr)])

#remove TES entries exceeding chromosome boundaries +/- 150 nt
out_of_bound_end   <- (end(tx_ends_gr) + 150) > seqlengths(tx_ends_gr)[as.character(seqnames(tx_ends_gr))]
out_of_bound_start <- (start(tx_ends_gr) - 150) < 1
tx_ends_gr <- tx_ends_gr[!out_of_bound_end & !out_of_bound_start]

# ---------------------------------------------------------------------------- #
# Clustering test dataset 3'ends

# Here we are collapsing 3' ends in the same position (or optionally within a
# window if variable `window` is set to greater than 0 nt). We also plot the 
# nucleotide frequency around the 3' ends.

# Optionally reduce n. of 3' end by their proximity
window <- 0 # increase 'window' to add proximity reduction 
# Reduce the ranges
reduced_gr <- reduce(tx_ends_gr, min.gapwidth = window)
reduced_gr$cluster_id <- paste0("clusterID_", 
                                seqnames(reduced_gr), 
                                ":", 
                                start(reduced_gr),"-", end(reduced_gr),
                                "|",
                                strand(reduced_gr))
three_prime_ends <- ifelse(strand(reduced_gr) == "+", end(reduced_gr), start(reduced_gr))

assembly_clusters <- GRanges(
  seqnames = seqnames(reduced_gr),
  ranges = IRanges(
    start = three_prime_ends,
    width = 1
  ),
  strand = strand(reduced_gr)
)
assembly_clusters$cluster_id <- reduced_gr$cluster_id 

stopifnot(all(seqnames(assembly_clusters) %in% names(reference_genome)))

# Plot nucleotide frequency
a <- plotNucleotides(assembly_clusters, reference_genome) +
  ggtitle("All test 3' ends") +
  guides(fill=guide_legend(title="Nt")) +
  ylim(0,1) + theme_bw() + theme(legend.title=element_blank())

# ---------------------------------------------------------------------------- #
# Input reference set (poly(A) databases)
combined_database <- import(reference_ends)

# Check if chr in the genome
if (any(!seqnames(combined_database) %in% names(reference_genome))) {
  warning("Removing reference 3' ends on chromosomes not present in the genome assembly provided")
  combined_database <- combined_database[seqnames(combined_database) %in% names(reference_genome),]
}

# Nucleotide probabilities in reference database
stopifnot(all(seqnames(combined_database) %in% names(reference_genome)))
c <- plotNucleotides(combined_database, reference_genome) +
  ggtitle("Reference 3' ends") +
  guides(fill=guide_legend(title="Nt")) +
  ylim(0,1) + theme_bw() + theme(legend.title=element_blank())

# ---------------------------------------------------------------------------- #
# Estimate quality of reference datasets 
# Nucleotide frequency of test 3' ends also found in the training database.

clustering_window <- 10
assembly_clusters.inReference <- subsetByOverlaps(assembly_clusters, 
                                                  combined_database, 
                                                  maxgap = clustering_window) # note optimized window

stopifnot(all(seqnames(assembly_clusters.inReference) %in% names(reference_genome)))

d <- plotNucleotides(assembly_clusters.inReference, reference_genome) +
  ggtitle("Test 3' ends in reference") +
  guides(fill=guide_legend(title="Nt")) +
  ylim(0,1) +
  theme_bw() +
  theme(legend.title=element_blank())

ntprobinput <- a | c | d

ggsave(filename = paste0(outfig, "/ntprob_input.pdf"),
       plot = ntprobinput, device = "pdf", width = 9, height = 3, dpi = 300)

# ---------------------------------------------------------------------------- #
# Add label for 3' ends that are presence in the reference database
# This step is key as the training set will be made by comparing the 3'ends
# in the reference to those not
assembly_clusters$inRef <- ifelse(
  assembly_clusters$cluster_id %in% 
    assembly_clusters.inReference$cluster_id, 
  TRUE,
  FALSE)

# ---------------------------------------------------------------------------- #
# Retrieve features to train the model or run it. 

features_in_assembly.3pends <- extract_Seqfeatures(clusters3seq = assembly_clusters,
                                                   genomeSeq = reference_genome, 
                                                   dwindow = 40,  #note optimized parameter
                                                   upwindow = 40)

features_in_assembly.3pends <- features_in_assembly.3pends %>%
  mutate(signalDistance = ifelse(is.na(signalDistance), 
                                 -99, 
                                 signalDistance))

# Polyadenylation signals (PAS) associated with 3'ends (in vs not in Reference) 

# PAS empirically ordered by frequency
pas_order <- c("AATAAA", "ATTAAA", "AATATA", "TATAAA", "TTTAAA", "AATACA",
               "AAGAAA", "CATAAA", "AATGAA", "GATAAA", "ACTAAA", "AATAGA",
               "not detected")

# Roerder PASs
features_in_assembly.3pends$polyA.Signal <- factor(features_in_assembly.3pends$polyA.Signal,
                                                   levels = pas_order)

e <- ggplot(features_in_assembly.3pends) +
  geom_bar(aes(x=inRef, 
               fill=polyA.Signal),position="fill") + 
  theme_classic() + 
  ylab("Fraction of 3' ends") + 
  xlab("TES in Reference") +
  theme(text = element_text(size=20))+ 
  scale_fill_viridis_d()+
  ggtitle("Polyadenylation signals in the test dataset")+
  guides(fill=guide_legend(title="PAS motif"))

ggsave(filename = paste0(outfig, "/pas_input.pdf"),
       plot = e, device = "pdf", width = 8, height = 8, dpi = 300)

# ---------------------------------------------------------------------------- #
## Random forest model training

### Prepare training and test sets

# Due to the significantly larger number of non-overlapping sets compared to the
# reference, we utilized a smaller sampling proportion of 3' ends not found in
# the reference for training, in order to balance the dataset, sampling a 
# fraction 0.8 of TRUE positives and 0.2 of the set of 3'ends not found in reference. 
# Depending on the size of the two sets, these paremeters should be adjusted
# accordingly to make the samples balanced.

# make training sets 
positive <- features_in_assembly.3pends %>% filter(inRef==TRUE)
negative <- features_in_assembly.3pends %>% filter(inRef==FALSE)

# sampling
set.seed(opt$seed)
split_positive <- sample.split(positive, SplitRatio = 0.8)
split_negative <- sample.split(negative, SplitRatio = 0.2)

#### train set 
# subset
positive_train <- subset(positive, split_positive == "TRUE")
negative_train <- subset(negative, split_negative == "TRUE")

# training
train <- rbind(positive_train, negative_train)
train$inRef <- as.factor(train$inRef)
rownames(train) <- make.unique( train$cluster_id)

train <- train %>%
  dplyr::select(-c(cluster_id)) %>%
  mutate(signalDistance = ifelse(is.na(signalDistance), 
                                 -99, 
                                 signalDistance))

rownames(features_in_assembly.3pends) <- make.unique(features_in_assembly.3pends$cluster_id)

# Testing data 
test_data <- features_in_assembly.3pends %>% 
  filter(!cluster_id %in% rownames(train)) %>% 
  dplyr::select(
    -c(
      cluster_id
    )
  ) %>%
  mutate(signalDistance = ifelse(is.na(signalDistance), 
                                 -99, 
                                 signalDistance))

# ---------------------------------------------------------------------------- #
### Training parameters

# In our analysis, we utilized the `ranger` function to develop a Random Forest
# classifier, named `classifier_mtry10`, trained on the `train` dataset with
# `inRef` as the classification variable and all other variables as predictors.
# The model construction involved generating 1000 trees (`num.trees = 1000`),
# each considering a random subset of 10 predictors at every split (`mtry = 10`),
# and assessing feature importance based on impurity reduction (`importance = "impurity"`).
# We also addressed potential class imbalance and refined our probability
# threshold for enhanced accuracy and sensitivity by applying weights
# (`weights <- c(FALSE = 0.03, TRUE = 1)`) via the `class.weights` parameter.
# Additionally, the `probability = TRUE` parameter was activated for obtaining 
# detailed class probabilities. This approach, especially given our clear
# a priori readouts (polyA signals and nucleotide distributions), enabled us to
# determine an optimal probability cutoff, effectively balancing accuracy and
# sensitivity in the model. The `seed` parameter was set to `NULL` to allow for
# random seeding, contributing to the model’s robustness. 

weights <- c(`FALSE` = 0.03, `TRUE` = 1)

if (file.exists(output_model)) {
  message(paste0("Model exists in ", output_model, ". Loading instead of re-running."))
  # The model can be loaded from output_model instead of being ran from scratch
  rf_classifier <- readRDS(output_model)

} else {
  rf_classifier <- ranger(inRef   ~ ., 
                          data=train, 
                          num.trees = 500, 
                          mtry=10, 
                          importance = "impurity",
                          seed = 42,
                          probability = TRUE,
                          class.weights = weights,
                          max.depth = 20,
                          num.threads = opt$threads)
  
  saveRDS(rf_classifier, output_model)
}

# ---------------------------------------------------------------------------- #
### Features of importance
feature_plot <- as.data.frame(rf_classifier$variable.importance) %>% 
  dplyr::mutate(feature = rownames(.)) %>%
  dplyr::rename(importance = 1) %>% 
  dplyr::mutate(importance = as.numeric(importance)) %>%
  filter(importance>0) %>%
  arrange(desc(importance)) %>% 
  ggplot(data=., aes(x=reorder(feature, importance) , 
                     y=importance)) +
  geom_bar(stat="identity", fill="dodgerblue") + 
  coord_flip() + 
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Features") + 
  ylab("Importance")

ggsave(filename = paste0(outfig, "/feature_importance.pdf"),
       plot = feature_plot, device = "pdf", width = 6, height = 7, dpi = 300)

# ---------------------------------------------------------------------------- #
### Accuracy and sensitivity on test data: recovered 3'ends per cut-off

# Load pre-trained model with above parameters 
model.predictions <- predict(rf_classifier, data = test_data)

# set TRUE classification probabilities for better readability
cut_offs <- c(0.001,0.01,0.05,0.1,0.25,0.5,0.75,0.95)
names(cut_offs) <- paste0( ">", cut_offs) 

# Plot recovered 3'ends per cut-off 
t1 <- lapply(cut_offs, function(x) {
  model.predictions$predictions %>% 
    as.data.frame(.) %>% 
    filter(`TRUE` > x) %>% 
    nrow(.)
}) %>% 
  unlist(.) %>% 
  as.data.frame(.) %>%
  dplyr::rename(num_of_3_ends = 1) %>%
  mutate(probability = gsub("prob", "", rownames(.)))  %>%
  ggplot(., aes(fill = probability, y = num_of_3_ends, x = probability)) +
  geom_bar(position = "dodge", stat = "identity") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_fill_brewer(palette="RdBu")+
  ylab("N. 3' ends above threshold")+
  xlab("Probability")

### Recovered poly(A) signals per cut-off

test_data$prediction_prob <- model.predictions$predictions[,2]
t2 <- lapply(cut_offs, function(x){
  test_data %>% 
    mutate(predicted_true = ifelse(prediction_prob>x, TRUE,FALSE)) %>% 
    group_by(polyAsignalClass, predicted_true) %>% 
    tally() %>% 
    group_by(predicted_true) %>%
    mutate(percentage = n / sum(n)) %>% 
    as.data.frame(.)}) %>% 
  do.call(rbind,.) %>% 
  mutate(probability = paste0("0." , 
                              str_split_fixed(rownames(.), 
                                              "\\.", 
                                              n=3)[,2])) %>%
  ggplot(., aes(fill=polyAsignalClass, 
                y=percentage, 
                x=predicted_true)) + 
  geom_bar(position="stack", 
           stat="identity") +
  theme_classic() + 
  facet_grid(cols=vars(probability)) + 
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust=1))+
  scale_fill_viridis_d()+
  ylab("Fraction of 3'ends")+
  xlab("Prediction")

# ---------------------------------------------------------------------------- #
### ROC & AUC

#with plotROC
roc.plot <- ggplot(test_data, aes(d=inRef, m=prediction_prob)) +
  geom_roc(labels=T, pointalpha=1, cutoffs.at = cut_offs, labelround = 3) + 
  theme_bw() +
  xlab("False Positive Fraction") + ylab("True positive Fraction")

auc_value <- calc_auc(roc.plot)[,3]
roc.plot <- roc.plot + ggtitle(paste("ROC curve on test set. AUC =",
                                     round(auc_value, 3)))

ggsave(filename = paste0(outfig, "/roc_testset.pdf"),
       plot = roc.plot, device = "pdf",
       width = 6, height = 6, dpi = 300)

### Sensitivity

# Check that the model is not overfitted. The predicted 3' ends should include
# both previously known 3' ends (i.e. those already in the reference) as well as new 3' ends.

t3 <- lapply(cut_offs, function(x){
  test_data %>% 
    mutate(predicted_true = ifelse(prediction_prob>x, TRUE,FALSE)) %>% 
    group_by(inRef, predicted_true) %>% 
    tally() %>% 
    group_by(predicted_true) %>%
    mutate(percentage = n / sum(n)) %>% 
    as.data.frame(.)}) %>% 
  do.call(rbind,.) %>% 
  mutate(probability = paste0("0." , str_split_fixed(rownames(.), "\\.", n=3)[,2])) %>%
  ggplot(., aes(fill=inRef, y=percentage, x=predicted_true)) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(cols=vars(probability))+
  scale_fill_manual(values=c("dodgerblue","firebrick"))+
  ylab("Fraction of 3'ends")+
  xlab("Prediction")

qctest<-(t1 | t2 | t3)

ggsave(filename = paste0(outfig, "/qc_testset.pdf"),
       plot = qctest, device = "pdf", width = 16, height = 4, dpi = 300)

### Nucleotide compositions per cut-off

# Make plot for each probability cut-off
nt.probs <- lapply(cut_offs, function(x){
  true_probabilities <- model.predictions$predictions %>% 
    as.data.frame(.) %>% dplyr::pull(2)
  test_data$true_probabilities <- true_probabilities 
  class.true.clusters <- test_data %>% 
    mutate(cluster_id = rownames(.)) %>% 
    filter(true_probabilities>x) %>% 
    dplyr::pull(cluster_id)
  class.false.clusters <- test_data %>% 
    mutate(cluster_id = rownames(.)) %>% 
    filter(!true_probabilities>x) %>% dplyr::pull(cluster_id)
  res <- list()
  message("Found in cut off:",  x  , " ","found: ", length(assembly_clusters[assembly_clusters$cluster_id %in% class.true.clusters]), " ","positive 3'-ends")
  res$true.3p_ends <- plotNucleotides(  assembly_clusters[assembly_clusters$cluster_id %in% class.true.clusters] , reference_genome) 
  res$false.3p_ends <- plotNucleotides(  assembly_clusters[assembly_clusters$cluster_id %in% class.false.clusters] , reference_genome)
  return(res)
})

# Print plot for 4 arbitrarily chosen probability cut-off
ntp1 <-  nt.probs$`>0.001`$true.3p_ends +ggtitle("Prob > 0.001")+ylim(0,1)+theme_bw()
ntp2 <-  nt.probs$`>0.01`$true.3p_ends +ggtitle("Prob > 0.01")+ylim(0,1)+theme_bw()
ntp3 <-  nt.probs$`>0.1`$true.3p_ends +ggtitle("Prob > 0.1")+ylim(0,1)+theme_bw()
ntp4 <-  nt.probs$`>0.25`$true.3p_ends +ggtitle("Prob > 0.25")+ylim(0,1)+theme_bw()

nttest <- ((ntp1 | ntp2) / (ntp3 | ntp4))
nttest
ggsave(filename = paste0(outfig, "/ntprob_testset.pdf"),
       plot = nttest, device = "pdf", width = 8, height = 8, dpi = 300)

# Define thresholds: an upper threshold such that 95% of the 3' ends in
# reference and contain a PAS (for increased stringency) are predicted to be true
# and a lower threshold such that 95% of the 3' ends that are not in reference
# and do not contain a PAS are predicted to be false. The upper threshold can be
# used to call a "bona fide" 3' end.

#find threshold of probability so that 95% of test inRef and PAS detected are kept
inRef.test_data <- test_data[test_data$inRef == TRUE & test_data$polyA.Signal != "not detected", ]
inRef.test_data.sorted <- inRef.test_data[order(inRef.test_data$prediction_prob, decreasing = TRUE), ]
true_index <- ceiling(0.95 * nrow(inRef.test_data))
true_threshold <- inRef.test_data.sorted$prediction_prob[true_index]

#find threshold of probability so that 95% of non inRef no PAS are kept
noRef.test_data <- test_data[test_data$inRef == FALSE, ]
noRef.test_data.sorted <- noRef.test_data[order(noRef.test_data$prediction_prob, decreasing = TRUE), ]
false_index <- ceiling(0.95 * nrow(noRef.test_data))
false_threshold <- noRef.test_data.sorted$prediction_prob[false_index]

my_message <- paste("The upper threshold is", round(true_threshold, 2),
                    "Over this threshold",
                    round(100*nrow(test_data[test_data$prediction_prob>true_threshold,])/nrow(test_data), 1),
                    "% of 3' ends are predicted to be true. The lower threshold is",
                    round(false_threshold, 6))

print(my_message)

# ---------------------------------------------------------------------------- #

# Output results 

# Here we run the trained classifier on the whole test dataset, produce some
#  QC plots and export gff files with the annotated 3' ends. 

## Run classifier and annotate 3' ends

# Run classification for the entire test dataset
rownames(features_in_assembly.3pends) <- make.unique(features_in_assembly.3pends$cluster_id)
model.predictions <- predict(rf_classifier, data = features_in_assembly.3pends)

# overlapping segments to clusters
features_in_assembly.3pends$probability_TRUE <- model.predictions$predictions[, 2]

assembly_clusters.annotated <- reduced_gr %>% 
  as.data.frame(.) %>%
  left_join(
    .,
    features_in_assembly.3pends %>%
      dplyr::select(polyAsignalClass, 
                    polyA.Signal, 
                    cluster_id, 
                    probability_TRUE,
                    inRef), by = "cluster_id") %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

# This step is necessary if the 3' ends were collapsed based on proximity at the
# start of the process, to make sure we go back to the original ranges of the input
# 3' ends
hits_to_clusters <- findOverlaps(tx_ends_gr, assembly_clusters.annotated)
annotated_clusters <- tx_ends_gr[queryHits(hits_to_clusters)]
mcols(annotated_clusters) <- mcols(assembly_clusters.annotated[subjectHits(hits_to_clusters)])

#add thresholds-based class
annotated_clusters$class <- ifelse(annotated_clusters$probability_TRUE >= true_threshold, "significant poly(A)", "unknown / non-significant")
annotated_clusters[annotated_clusters$probability_TRUE <= false_threshold,]$class <- "significant non-poly(A)"

### Fraction of annotated positive based on ref per cut-off

probability_cut_off <- c(0.001,0.01,0.05,0.1,0.25,0.5,0.75,0.95)
names(probability_cut_off) <- paste0(">", probability_cut_off)
t4 <- lapply(probability_cut_off, function(x){
  a1 <- annotated_clusters %>% as.data.frame(.) %>%
    filter(inRef==TRUE) %>% 
    mutate(isFound_inRef_and_Forest = ifelse(inRef==TRUE & probability_TRUE>x, TRUE, FALSE)) %>%
    group_by(isFound_inRef_and_Forest) %>% 
    tally() %>%
    mutate(percentage = n/sum(n)) %>%
    as.data.frame(.) 
  a1$probability <- x
  return(a1)
}) %>% do.call(rbind,.) %>%
  ggplot(., aes(fill=isFound_inRef_and_Forest, y=percentage, x=as.factor(probability))) + 
  geom_bar(position="stack", stat="identity") + 
  theme_classic()+ggtitle("Test ends in reference")+
  xlab("Cutoff")+ ylab("Fraction of 3' ends") +
  scale_fill_manual(values=c("dodgerblue","firebrick"))+theme(legend.title=element_blank())

t6 <- lapply(probability_cut_off, function(x){
  a1 <- annotated_clusters %>% as.data.frame(.) %>%
    filter(inRef==FALSE) %>% 
    mutate(isFound_in_Forest = ifelse(inRef==FALSE & probability_TRUE>x, TRUE, FALSE)) %>%
    group_by(isFound_in_Forest) %>% 
    tally() %>%
    mutate(percentage = n/sum(n)) %>%
    as.data.frame(.) 
  a1$probability <- x
  return(a1)
}) %>% do.call(rbind,.) %>%
  ggplot(., aes(fill=isFound_in_Forest, y=percentage, x=as.factor(probability))) + 
  geom_bar(position="stack", stat="identity") + 
  theme_classic()+ggtitle("All test dataset")+
  xlab("Cutoff")+ ylab("Fraction of 3' ends") +
  scale_fill_manual(values=c("dodgerblue","firebrick"))+theme(legend.title=element_blank())

### Recovered poly(A) signals per cut-off

t8 <- lapply(cut_offs, function(x){
  as.data.frame(annotated_clusters) %>% 
    mutate(predicted_true = ifelse(probability_TRUE>x, TRUE,FALSE)) %>% 
    group_by(polyAsignalClass, predicted_true) %>% 
    tally() %>% 
    group_by(predicted_true) %>%
    mutate(percentage = n / sum(n)) %>% 
    as.data.frame(.)}) %>% 
  do.call(rbind,.) %>% 
  mutate(probability = paste0("0." , 
                              str_split_fixed(rownames(.), 
                                              "\\.", 
                                              n=3)[,2])) %>%
  ggplot(., aes(fill=polyAsignalClass, 
                y=percentage, 
                x=predicted_true)) + 
  geom_bar(position="stack", 
           stat="identity") +
  theme_classic() + 
  facet_grid(cols=vars(probability)) + 
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust=1))+
  scale_fill_viridis_d()+ggtitle("Test dataset") +
  xlab("Prediction")+ ylab("Fraction of 3' ends")

qc_cutoff <- t4 | t6 | t8
ggsave(filename = paste0(outfig, "/qc_bycutoff.pdf"),
       plot = qc_cutoff, device = "pdf",
       width = 12, height = 6, dpi = 300)

nt.probs <- lapply(cut_offs,
                   function(x){
  res <- plotNucleotides(annotated_clusters[annotated_clusters$probability_TRUE > x],
                         reference_genome) 
  return(res)
})

ntp5 <- nt.probs$`>0.001` + ggtitle("Prob > 0.001") + ylim(0, 1) + theme_bw()
ntp6 <- nt.probs$`>0.01` + ggtitle("Prob > 0.01") + ylim(0, 1) + theme_bw()
ntp7 <- nt.probs$`>0.1` + ggtitle("Prob > 0.1") + ylim(0, 1) + theme_bw()
ntp8 <- nt.probs$`>0.25` + ggtitle("Prob > 0.25") + ylim(0, 1) + theme_bw()

nt_cutoff <- (ntp5 | ntp6 | ntp7 | ntp8)
ggsave(filename = paste0(outfig, "/ntprob_bycutoff.pdf"),
       plot = nt_cutoff, device = "pdf",
       width = 16, height = 4, dpi = 300)

## Export files
rtracklayer::export.gff2(annotated_clusters, output_gtf)

