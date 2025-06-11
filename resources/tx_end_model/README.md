# Datasets for random forest classification and sequence analysis of transcriptome assemblies 3’ ends

Authors: Rodrigo Pracana, Carlos Alfonso-Gonzalez, Ivano Legnini

This directory includes datasets necessary to train a random forest classifier to identify true 3' ends from CFC-seq presented [here](../../code/tx_end_model/). Specifically, we use a training set with the 3’ ends generated in the SALA demo (the file [Neuron_series_demo.end3.bed](Neuron_series_demo.end3.bed.gz)). The classifier is trained using a curated database (hereon the "reference" dataset) of 3’ ends identified by FLAM-seq and 3p-seq of iPS and neuron organoids (Legnini et al., 2019 and Alfonso-Gonzalez et al., 2023), given in the file [`reference_3prime.bed`](reference_3prime.bed).

## Citation

Legnini, I., Alles, J., Karaiskos, N. et al. FLAM-seq: full-length mRNA sequencing reveals principles of poly(A) tail length control. Nat Methods 16, 879–886 (2019). https://doi.org/10.1038/s41592-019-0503-y

Alfonso-Gonzalez, C., Legnini, I., Holec, S., et al. Sites of transcription initiation drive mRNA isoform selection. Cell 186, 2438-2455.e22 (2023). https://doi.org/10.1016/j.cell.2023.04.012
