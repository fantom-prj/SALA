#!/usr/bin/env bash

set -euo pipefail

# Load conda
source ~/miniconda3/etc/profile.d/conda.sh

conda activate snakemake-slurm

snakemake \
    --executor slurm \
    --snakefile snakemake/Snakefile \
    --workflow-profile profile/test \
    --configfile config/test/config.test.yml \
    --jobs 6 \
    --cores 12 \
    --keep-going \
    --retries 2 \
    --rerun-incomplete
