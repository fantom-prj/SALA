# Tail-trimmer v1.4

Authors: Chung Chau Hon

Tail trimmer remove the poly(A) tail from the 3' end of reads from FASTQ file. It removes the polyA together with some base-calling error at the 3' end. The last 5 As, if present, are added back to the reads for softclipping during mapping. Besides the trimmed FASTQ, also returns the information of the trimmed region including the number of As. This step is optional but useful if the length of As signal is too long after in vitro poly(A) tailing. This step also allows adding a prefix to the reads.

The script can be called with the following parameters below. 

```sh

  ./tail_trimmer_v1.4.pl \
    --rescue_tail_nt=5 \
    --revise_read_prefix=iPSC1 \
    --fastq_path=./iPSC1.fq \
    --out_dir=./iPSC1

```

