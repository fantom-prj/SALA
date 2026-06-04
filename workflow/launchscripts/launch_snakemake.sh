#!/bin/bash

set -e  # Stop execution if any command fails

usage() {
    cat <<EOF
Usage: $(basename "$0") <config_file> [profile] [n_jobs] [n_threads] 

Description:
    config_file (required): path of the configuration file.
    profile (optional): path of the profile configuration for the execution. (default: profile)
    n_jobs (optional): max number of jobs that will be orchestrated by snakemake.
    n_threads (optional): max number of threads that will be orchestrated by snakemake.

Options:
  -h, --help    Show this help message and exit

EOF
}

# Handle help option
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    usage
    exit 0
fi

if [[ $# -lt 1 || $# -gt 4 ]]; then
    usage >&2
    exit 1
fi

# Load and activate conda
source ~/miniconda3/etc/profile.d/conda.sh
export CONDA_CHANNEL_PRIORITY=flexible
conda activate snakemake-slurm

config_file=$1
profile=${2-"profile"}
jobs=${3-8}
cores=${4-16}

echo "$(date): Snakemake installs conda enviroments if needed"

logs_dir=logs/snakemake
mkdir -p $logs_dir

#Launch snakemake
echo "$(date): Launch snakemake pipeline"
# MAX: --jobs 100 --cores 100
snakemake --executor slurm \
    --snakefile snakemake/Snakefile \
    --workflow-profile $profile \
    --slurm-logdir $logs_dir \
    --configfile $config_file \
    --jobs $jobs --cores $cores --keep-going \
    --retries 2 --rerun-incomplete

echo "$(date): Pipeline finished"
echo "$(date): DONE"
