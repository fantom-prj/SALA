# Helper functions 

require(LATER)
require(dplyr)
require(GenomicRanges)
library(caTools)
library(dplyr)
library(randomForest)
library(stringr)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(randomForest)
library(ranger)

#' Extract features from clusters
#'
#' @param clusters3seq
#' @param genomeSeq
#' @param dwindow
#' @param upwindow
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors Biostrings
#' @examples
#' 

extract_Seqfeatures <- function(clusters3seq, genomeSeq, dwindow, upwindow){
  # extend 3'coordinates for feature extraction
  message("Extending clusters in flanking windows")
  downstreamSite <-
    extend(clusters3seq, upstream = 0, downstream = dwindow)
  upstreamSite <-
    extend(clusters3seq, upstream = upwindow, downstream = 0)
  # make sequences
  seqsUpstream <- Biostrings::getSeq(genomeSeq, upstreamSite)
  names(seqsUpstream) <- upstreamSite$cluster_id
  seqsDownstream <- Biostrings::getSeq(genomeSeq, downstreamSite)
  names(seqsDownstream) <- downstreamSite$cluster_id
  message("retrieve nucleotide frequencies")
  seqsUpstream.dt <-
    as.data.frame(alphabetFrequency(seqsUpstream, 
                                    baseOnly = TRUE, 
                                    as.prob =TRUE)) %>% 
    dplyr::select(-other)
  colnames(seqsUpstream.dt) <-
    paste0("nt.prob.upstream", colnames(seqsUpstream.dt))
  seqsUpstream.dt$cluster_id <- names(seqsUpstream)
  seqsDownstream.dt <-
    as.data.frame(alphabetFrequency(seqsDownstream, 
                                    baseOnly = TRUE, 
                                    as.prob =
                                      TRUE)) %>% 
    dplyr::select(-other)
  colnames(seqsDownstream.dt) <-
    paste0("nt.prob.downstream", colnames(seqsDownstream.dt))
  seqsDownstream.dt$cluster_id <- names(seqsDownstream)
  clusterAContent <-
    left_join(seqsUpstream.dt, 
              seqsDownstream.dt, 
              by = "cluster_id")
  ## annotate poly(A) signals in upstream
  message("Retrieve poly(A) signals")
  message("Retrieve poly(A) signals Upstream")
  signalsUpStream <-
    getPolyAsignal(APAmotifs, 
                   hexamers,
                   seqsUpstream)
  signalsUpStream <- lapply(signalsUpStream,
                            as.data.frame)
  signalsUpStream <-
    bind_rows(signalsUpStream, .id = "column_label")
  message("Annotate poly(A) signals to clusters")
  signalsUpStream <-
    left_join(clusterAContent,
              signalsUpStream %>%
                dplyr::rename(cluster_id = names)) %>%
    dplyr::rename(polyA.Signal = column_label) %>%
    mutate(
      polyA.Signal = ifelse(is.na(polyA.Signal), 
                            "not detected", 
                            polyA.Signal),
      signalDistance = end - 46
    ) %>% 
    dplyr::select(-c(width, start, end))
  ## remove duplicated clusters with multiple assignments keeping signal closer to cleavage site
  signalsUpStream <- signalsUpStream %>%
    arrange(abs(signalDistance), 
            polyA.Signal) %>%
    distinct(cluster_id, 
             .keep_all = TRUE)
  # retrieve counts
  signalsUpStream <- left_join(
    signalsUpStream ,
    mcols(clusters3seq) %>% 
      as.data.frame(.) %>%
      dplyr::rename(cluster_id = cluster_id),
    by = "cluster_id"
  )
  # add oligo kmers
  message("Search for kmers")
  oligoDT.kmers <-
    lapply(c(10, 15, 20), function (x) {
      strrep("A", x)
    })
  kmerClustersUpstream <-
    lapply(oligoDT.kmers, 
           kmer_scan, 
           seqsUpstream)
  kmerClustersDownstream <-
    lapply(oligoDT.kmers, 
           kmer_scan, 
           seqsDownstream)
  message("Prepare table")
  signalsUpStream <-
    signalsUpStream %>% mutate(
      Ant.8mer.upstream = ifelse(cluster_id %in% kmerClustersUpstream[1], TRUE, FALSE),
      Ant.10mer.upstream = ifelse(cluster_id %in% kmerClustersUpstream[2], TRUE, FALSE),
      Ant.10mer.upstream = ifelse(cluster_id %in% kmerClustersUpstream[3], TRUE, FALSE),
      Ant.8mer.downstream = ifelse(cluster_id %in% kmerClustersDownstream[1], TRUE, FALSE),
      Ant.10mer.downstream = ifelse(cluster_id %in% kmerClustersDownstream[2], TRUE, FALSE),
      Ant.10mer.downstream = ifelse(cluster_id %in% kmerClustersDownstream[3], TRUE, FALSE)
    )  %>%
    mutate(polyAsignalClass =
             ifelse(polyA.Signal %in% c("AATATA", "ATTAAA", "AATAAA"),
                    polyA.Signal,
                    ifelse(polyA.Signal == "not_detected",
                           "not_detected",
                           "non_canonical")))
  return(signalsUpStream)
}

APAmotifs <-
  c(
    "AATAAA",
    "ATTAAA",
    "AATATA",
    "AAGAAA",
    "AATACA",
    "AATAGA",
    "AATGAA",
    "ACTAAA",
    "CATAAA",
    "GATAAA",
    "TATAAA",
    "TTTAAA"
  )
hexamers <- list()

#' Title
#'
#' @param x
#' @param upstream
#' @param downstream
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors Biostrings
#' @examples
extend <- function(x, upstream=0, downstream=0)   {
  if (any(BiocGenerics::strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- BiocGenerics::strand(x) == "+" | BiocGenerics::strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

#' Get PAS for 3' ends
#'
#' @param APAmotifs
#' @param hexamers
#' @param seqs
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors Biostrings
#' @examples
getPolyAsignal <- function(APAmotifs, hexamers, seqs){
  if(rlang::is_empty(APAmotifs)){ # no motifs remaining
    return(hexamers)
  } else {
    apaMotif = APAmotifs[1] #get first entry in list
    dnaseq <- DNAString(apaMotif) #canonical 1
    hexamers[[apaMotif]] <-unlist(vmatchPattern(dnaseq, seqs))
    seqs <- seqs[!names(seqs)  %in% names( hexamers[[apaMotif]] ),]
    APAmotifs <- APAmotifs[!APAmotifs %in%  apaMotif] # remove motif from list
    getPolyAsignal(APAmotifs, hexamers, seqs)
  }
}

#' Prepare clusters for processing
#'
#' @param clusters3seq
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors Biostrings
#' @examples
prepare_clusters <- function(clusters3seq){
  clusters3seq <-
    clusters3seq %>% makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  end(clusters3seq) <- start(clusters3seq)
  seqlevelsStyle(clusters3seq) <- "UCSC"
  clusters3seq$cluster_name <-
    paste0(seqnames(clusters3seq),
           ":",
           start(clusters3seq),
           ":",
           BiocGenerics::strand(clusters3seq))
  names(clusters3seq) <- clusters3seq$cluster_name
  clusters3seq <- keepStandardChromosomes(clusters3seq, pruning.mode = "coarse")
  return(clusters3seq)
}

#' Scan kmers for features
#'
#' @param pattern
#' @param seqsUp
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors
#' @examples
kmer_scan <- function(pattern,seqsUp){
  upstreamA <- unlist(vmatchPattern("AAAAAAAA", seqsUp)) %>% as.data.frame(.) %>% dplyr::pull(names)
  return(upstreamA)
}

## EXPORT TABLE function 
#' Prepare Table
#'
#' @param tes_from_assembly: TES from LATER prepareLinksDatabase Function 
#' @param polya_database: polyA database used for the analysis 
#' @param features3p_data : Features used in random forest classification. 
#'
#' @return
#' @export
#'
#' @examples
prepareTable <- function(tes_from_assembly, polya_database, features3p_data){ 
  exportObject <- list()
  x1 <- subsetByOverlaps(tes_from_assembly, polya_database)
  clusters_found_in_reference <- x1$cluster_id
  clusters_predicted_true <- features3p_data %>% 
    filter(predicted == TRUE ) %>% 
    dplyr::pull(cluster_id)
  exportObject$export_coordinates <- tes_from_assembly %>% 
    as.data.frame(.) %>% 
    mutate(tes_inReference = ifelse(cluster_id %in% clusters_found_in_reference, TRUE, FALSE),
           tes_predicted = ifelse(cluster_id %in% clusters_predicted_true, TRUE, FALSE)) %>% 
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE) 
  # prepare table with features 
  exportObject$feature_table <- features3p_data %>% 
    mutate(tes_inReference = ifelse(cluster_id %in% clusters_found_in_reference, TRUE, FALSE),
           tes_predicted = ifelse(cluster_id %in% clusters_predicted_true, TRUE, FALSE)) %>% 
    dplyr::select(cluster_id,
                  polyA.Signal,
                  signalDistance,
                  tes_predicted,
                  tes_inReference
    )
  return(exportObject)
}

getNucleotideProb <- function(cluster_regions, genome){
  flank = 50
  pas.flank <- flank(cluster_regions , width = flank, both = TRUE)
  pas.flank_seq <- Biostrings::getSeq(genome, pas.flank)
  prob <- consensusMatrix(pas.flank_seq,as.prob = TRUE)
  df <- reshape2::melt(prob)
  prob1 <- prob
  prop3 <- list(prob1, prob)
  df <- reshape2::melt(prop3)
  levels(df$Var1)[4] <- 'U'
  df = subset(df, grepl("A|C|G|U", Var1))
  return(df)
}
# plot
#' Prepare nucleotide frequencies for plotting
#'
#' @param x
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors ggplot2
#' @examples
nucleotidePlot <- function(x) {
  nt.plot <- ggplot(x, aes(Var2 - 50, value)) +
    geom_line(aes(color = Var1), alpha = 1, size = 1) +
    scale_color_manual(
      values =
        c(
          "A" = "red",
          "C" = "blue",
          "G" = "#ECC54E",
          "U" = "#A76BCF"
        )
    ) +
    xlab("Distance to Poly(A) site") +
    ylab("Relative frequency") +
    theme_classic() +
    ggtitle("Nucleotide probabilities") +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    )
  return(nt.plot)
}


#' Plot nucleotide frequencies
#'
#' @param three_prime_ends.gr
#' @param genome
#'
#' @return
#' @export
#' @import dplyr GenomicFeatures GenomicAlignments S4Vectors
#' @examples
plotNucleotides <- function(three_prime_ends.gr, genome) {
  nucleotide.probs <- getNucleotideProb(three_prime_ends.gr, genome)
  nucleotidePlot(nucleotide.probs)
  
}

#' Trim to most 3p nucleotide
#'
#' @param GenomicRanges object 
#'
#' @return GenomicRanges object single nt
#' @export
#'
#' @examples
trim_single_nucleotide <- function(reads.gr) {
  pos <- reads.gr[ BiocGenerics::strand(reads.gr)== "+",]
  start( pos) <-  end(pos)
  neg <- reads.gr[BiocGenerics::strand(reads.gr)== "-",]
  end( neg) <-  start(neg)
  trim.gr <- c(pos,neg )
  return(trim.gr)
}

#' Get sequence for motif enrichment
#'
#' @param granges_object
#' @param genome
#' @param output_file
#'
#' @return
#' @export
#' @examples
getPASeq <- function(granges_object, genome, output_file) {
  # Copy the input GRanges object to avoid modifying the original
  obj <- granges_object
  # Adjust start and end positions
  ends <- ifelse(strand(obj) == "+", end(obj), start(obj))
  start(obj) <- ends - 50
  end(obj) <- ends + 50
  # Extract sequences from the genome
  obj <- getSeq(genome, obj)
  # Set names based on cluster.ID
  names(obj) <- seq(1,length(obj),1)
  # Write to a FASTA file
  writeXStringSet(obj, output_file, format = "fasta")
  return(obj)
}


#' Plot motif frequencies
#'
#' @param motif
#' @param seq_objects_list
#' @param categories
#'
#' @return
#' @export
#' @examples
motif_plot_multi <- function(motif, seq_objects_list, categories) {
  motif_frequencies_list <- vector("list", length = length(seq_objects_list))
  
  for (j in seq_along(seq_objects_list)) {
    seq_object <- seq_objects_list[[j]]
    motif_frequencies <- integer(length(seq_object[[1]]))
    
    for (i in 1:(length(seq_object[[1]]) - length(motif) + 1)) {
      motif_counts <- vcountPattern(
        motif,
        subseq(seq_object, start = i, end = i + length(motif) - 1),
        fixed = FALSE
      )
      motif_frequencies[i] <- sum(motif_counts) / length(seq_object)
    }
    
    motif_frequencies_list[[j]] <- motif_frequencies
  }
  
  # Create a data frame for plotting
  motif_data <- data.frame(
    Position = 1:length(motif_frequencies_list[[1]]),
    Frequency = unlist(motif_frequencies_list),
    Category = rep(categories, each = length(motif_frequencies_list[[1]]))
  )
  
  # Plot the motif frequencies
  motif_plot <- ggplot(motif_data, aes(x = Position - 50, y = Frequency, fill = Category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("orange","dodgerblue","firebrick","seagreen", "darkblue", "lightgreen", "darkorange4", "darkgreen"))+
    labs(
      title = paste(as.character(motif), " Frequency"),
      x = "Position relative to cluster end",
      y = "Motif frequency"
    )+
    theme_bw()
  
  return(motif_plot)
}


#' Find top enriched kmers
#'
#' @param X
#' @param Y
#' @param seq_object_list
#'
#' @return
#' @export
#' @examples
find_top_kmers <- function(X, Y, seq_object_list) {
  # Function to calculate kmer frequencies within a single sequence
  calculate_kmer_frequencies <- function(seq, X) {
    all_kmers <- apply(expand.grid(replicate(X, c("A", "T", "C", "G"), simplify = FALSE), stringsAsFactors = FALSE), 1, paste, collapse = "")
    kmer_frequencies <- numeric(length(all_kmers))
    
    for (i in seq_along(all_kmers)) {
      kmer <- all_kmers[i]
      kmer_counts <- vcountPattern(
        as.character(kmer),
        seq,
        fixed = FALSE
      )
      kmer_frequencies[i] <- sum(kmer_counts)
    }
    
    return(data.frame(
      Kmer = all_kmers,
      Frequency = kmer_frequencies,
      Object = rep(as.character(substitute(seq)), length(all_kmers))
    ))
  }
  
  result_list <- lapply(seq_object_list, function(seq_object) {
    kmer_frequencies_list <- calculate_kmer_frequencies(seq_object, X)
    return(bind_rows(kmer_frequencies_list))
  })
  
  # Combine the data frames for all objects
  result <- bind_rows(result_list)
  
  # Group by Object and Kmer, then summarize to get the total frequency of each kmer per object
  top_kmers <- result %>%
    group_by(Object, Kmer) %>%
    summarize(TotalFrequency = sum(Frequency)) %>%
    group_by(Object) %>%
    top_n(Y, wt = TotalFrequency) %>%
    ungroup()
  
  return(top_kmers)
}


