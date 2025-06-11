# Random Forest classification and sequence analysis of transcriptome assemblies 3’ ends

Authors: Rodrigo Pracana, Carlos Alfonso-Gonzalez, Ivano Legnini

In order to call bona fide 3’ ends for the transcripts detected by CFC-seq, we use a random forest classifier. We give an example that can be ran with the datasets available [here](../../resources/tx_end_model/). Specifically, we use a training set with the 3’ ends generated in the SALA demo (the file [Neuron_series_demo.end3.bed](../../resources/tx_end_model/Neuron_series_demo.end3.bed.gz)). The classifier trained using a curated database (hereon the "reference" dataset) of 3’ ends identified by FLAM-seq and 3p-seq of iPS and neuron organoids (Legnini et al., 2019 and Alfonso-Gonzalez et al., 2023), given in the file [`reference_3prime.bed`](../../resources/tx_end_model/reference_3prime.bed).

The  script used to train and evaluate the model is [`RunForest.R`](RunForest.R), and it relies on functions defined in the the file [`code/rf_helpers.R`](code/rf_helpers.R). The libraries used by the script can be installed with the script [`code/install_packages.R`](code/install_packages.R).

## Model training and useage

We use the genome assembly [GRCh38 GCA_000001405.15](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/), which we subset to the standard chromosomes using [seqtk](https://github.com/lh3/seqtk):

```sh

seqtk subseq -l 80 GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz dataset/chr.txt \
  | bgzip -c > GRCh38_no_alt_analysis_set_GCA_000001405.std_chr.fa.gz

```

The first step of the `RunForest.R` script is to annotate the sequence features around each 3’ end, and define whether or not each end is supported by the reference dataset. The features included the presence, type and position of polyadenylation signals and the nucleotide composition of a 50 nt window upstream and downstream of each 3’ end.

The second step is to train a random forest classifier using the ranger package (v. 0.16.0), considering the 3’ ends found in the reference as true positives.

The last step is to use the trained model to assign a probability score to each 3’ end in the CFC-seq dataset. We define an upper threshold at which 95% of 3’ ends found in the reference and containing a polyadenylation signal were classified as significantly true 3’ ends, and a lower threshold at which 95% of all 3’ ends classified significantly false 3’ ends. 

The script can be called with the following parameters below. Note that due to the significantly larger number of non-overlapping sets compared to the reference, we balance the sampling to a fraction 0.8 of TRUE positives and 0.2 of the set of 3’ ends not found in reference. The seed parameter allows the random sampling to be replicable between runs.

```sh

mkdir output/fig
Rscript RunForest.R \
  --genome GRCh38_no_alt_analysis_set_GCA_000001405.std_chr.fa.gz \
  --reference_set reference_3prime.bed \
  --test_set Neuron_series_demo.end3.bed.gz \
  --output_gtf output/Neuron_series_demo.3pclusters.gtf \
  --output_model output/Neuron_series_demo.model.RDS \
  --output_figure output/fig \
  --threads 2 \
  --seed 123

```

Each 3’ end is given a probability that is a true 3’ end. An upper threshold is set such that 95% of the 3' ends in reference that contain a PAS (for increased stringency) are predicted to be true and a lower threshold such that 95% of the 3' ends that are not in reference and do not contain a PAS are predicted to be false.

The main output of this script is a GTF file (set in `--output_gtf`) with the annotated 3’ end, which include a field indicating this probability (`probability_TRUE`), and the field `class` assigned to `"significant poly(A)"` if above the upper threshold, `"significant non-poly(A)"` if below the lower threshold and `"unknown / non-significant"` if in between the two.

The script also outputs a directory (set in `--output_figure`) with figures exploring the features of the 3' ends and the accuracy of the model.

The model is itself saved as an RDS file (set in `--output_model`). If a file of the chosen name already exists, the script will simply load it instead of re-running the training step. 
