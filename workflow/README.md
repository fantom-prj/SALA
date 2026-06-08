# SALA SLURM-Snakemake workflow

This workflow integrates SCAFE[^scafe] (v1.0.2) and SALA in a SLURM-ready snakemake pipeline to assemble a transcript and gene model starting from long-read alignments in BAM format. In particular, it goes through the steps 3., 4. and 5. of the SALA wiki[^sala-wiki]. 
This workflow allows to process independent libraries in parallel, through separated jobs orchestrated by snakemake. For each library, independent sub-jobs can be run in parallel as well. This is also handled by snakemake, depending on the available resources on your cluster.

## Dependencies
The workflow is designed to run on a SLURM-based HPC environment. The only prerequisite is a Conda-compatible package manager[^conda]. We recommend using Conda ≥ 23.0. All required channels and package dependencies are specified in the workflow environment files (see below).

## Configuration

### General config

The configuration files needed for the workflow are found in the config folder. the ```config.yaml``` file is the entrypoint for all information regarding the location of input, results, library configuration and parameters. 

```yaml
samplesheet: Path of the library table (see below)
results_dir: Path to the folder where results will be stored
reference_dir: Path where processed files from reference annotation will be stored
logs_dir: Path to the folder where logs will be stored
fasta_genome: Path of the reference genome (.fa.gz). It needs to be the same the BAM reads are aligned to.
mask_bed: Path of the promoter/enhancer annotation on the reference genome (.bed)
chr_sizes: Path of a tab-separated file containing the chromsome size in format of chrom\tsize (see example in test/resources/chrom.sizes_major.tsv)
gtf_annotation: Path of the reference transcriptome annotation (.gtf.gz)
genome_anno_id: Id for the genome+annotation folder (for SCAFE and SALA internal usage)
genome_code: Code name for the genome (for SALA internal usage)
library_config: Path of the library configuration file (see below)
params_set_file: Path of the parameters set file
params_set: Choice of parameters set (must be column in params_set_file)
```

A basic set of resources, based on GRCh38 and GENCODE 47, can be found in the ```test/resources``` folder.
Other genome-transcript annotation sets, including cCRE annotation from SCREEN, are available for human and mouse and can be downloaded at the following URLs:

| Genome | Transcript model | Mask BED | Download |
|----------|------------------|----------|----------|
| hg38 | GENCODE_v39 | SCREEN cCRE hg38 | [hg38.gencode_v39](https://drive.google.com/file/d/1Dx77UntpJZCZpffapQ33nXksXVivbMMX/view?usp=drive_link) |
| hg38 | GENCODE_v47 | SCREEN cCRE hg38 | [hg38.gencode_v47](https://drive.google.com/file/d/1QTOvFqD_yPMZtur5EP3pNrteFKudTjrc/view?usp=drive_link) |
| mm10 | GENCODE_vM25 | SCREEN cCRE mm10 | [mm10.gencode_vM25](https://drive.google.com/file/d/1WLF4cCQTcKbI8JzjhMcwI92daGQYzSBn/view?usp=drive_link) |
| mm39 | GENCODE_vM36 | SCREEN cCRE mm10_liftover_to_mm39 | [mm39.gencode_vM36](https://drive.google.com/file/d/1WOO3Sx33VqfqMU5_6C3-uvLzTQ7xTqAa/view?usp=drive_link) |

You may also use ```gdown --fuzzy```[^gdown] to download these files.

The ```params.csv``` file contains pre-defined parameters set (default and sensitive) for the various steps of the workflow. To define your own set of parameters, it is possible to add an additional column to the csv, and to change the parameters set choice accordingly in the config.yml.

The ```profile/config.yml``` file is related to resources configuration. Here it is possible to specify appropriate RAM requirements for the various steps. To handle out of memory issues, by default, the workflow will perform two retries for each step, each time doubling the allocated memory. The mem_mb default values have been defined according to tests on libraries with 60 to 250 million reads. We recommend to change them accordingly to your library size needs.

### Data config
The ```sample_table.csv``` contains the list of all libraries to be processed in parallel, and has to be provided in the sample table specified in the configuration (samplesheet), in this format (header included):
```csv
library_id,bam_list,rep_list
neuronset_id,/path/of/neuset/run/file.tsv,/path/of/neuset/replicates/file.tsv
thp1_id,/path/of/thp1/run/file.tsv,/path/of/thp1/replicates/file.tsv
```
The granularity of the sample set can be decided according to the user needs. Independent transcript and gene models will be produced for each of the libraries listed in the samplesheet.

For each library in the samplesheet, a replicate file and a run file are needed, both in TSV format.

```replicates.tsv```
The replicate file contains information on how to group the replicates in the library. The format is as follows (header included): 
```tsv
# library_prefix: prefix for the replicate, not containing underscores ("_")
# sample_ID: identifier for a sample in the library
# sample_rep: replicate of the sample
library_prefix	sample_ID	sample_rep
ipsc11	iPSC	iPSC_rep1
ipsc21	iPSC	iPSC_rep2
nsc11	NSC	NSC_rep1
nsc21	NSC	NSC_rep2
neu11	Neuron	Neuron_rep1
neu21	Neuron	Neuron_rep2
```

```runs.tsv```
The run file contains information on where the BAM files are stored, for each replicate in the library defined in the replicates table. It has to be headerless, and contain the following columns:
```tsv
# col1: prefix for the replicate
# col2: identifier for the replicate in the library
# col3: location of the run's alignment BAM file
ipsc11	iPSC_rep1_run1	test/iPSC_rep1_run1_subset.sorted.bam
ipsc21	iPSC_rep2_run1	test/iPSC_rep2_run1_subset.sorted.bam
nsc11	NSC_rep1_run1	test/NSC_rep1_run1_subset.sorted.bam
nsc21	NSC_rep2_run1	test/NSC_rep2_run1_subset.sorted.bam
neu11	Neuron_rep1_run1	test/Neuron_rep1_run1_subset.sorted.bam
neu21	Neuron_rep2_run1	test/Neuron_rep2_run1_subset.sorted.bam
```
For more details, refer to the SALA wiki.

```library.config.yml```
In this file, the following information about the libraries have to be specified:
* `internal_priming`: true if reads could be subjected to internal priming due to library preparation, false otherwise
* `5pr_confident`: true if the reads have confident 5' ends, false otherwise. In this latter case, external TSS clusters need to be provided as a .bed.gz file. Pre-set cluster coordinates for human (hg38) and mouse (mm10) can be found in the workflow/resources/CTSS_cluster folder.
* `5pr_cluster`: location of the externally provided TSS clusters file. Mandatory if 5pr_confident is false.

An example of a config file for a not confident 5' library can be found under config/test.

## How to run

### 1. Environment setup
* Make sure SALA scripts are executable
```sh
chmod 755 -R ./code/
```
* Move into the workflow directory
```sh
cd workflow
```
* Create the snakemake-slurm environment from the environment.yml file in the snakemake folder. This is the minimal environment that is needed for snakemake to run (with the slurm plugin) and to independently create and manage the actual environment needed for the run (in snakemake/envs/sala_wf_env.yml).
```sh
conda env create -f snakemake/environment.yml
```

### 2. File setup
You may copy the test files and fill the fields according to your resources and data.

* Prepare the library tables (sample_table.csv, replicates.tsv, runs.tsv)
```sh
cp config/test/all_replicates.tsv config/replicates.tsv 
# Edit your replicates.tsv

cp config/test/all_runs.tsv config/runs.tsv 
# Edit your runs.tsv

cp config/test/sample_table_all.csv config/sample_table.csv 
# Edit your sample_table.csv
```

* Set up library 5' confidence and internal priming details in the libraries configuration file (library.config.yml). 
```sh
# In case your libraries are 5' confident
cp config/test/library.config.yml config/library.config.yml

# In case your libraries are not 5' confident
cp config/test/library_notconf.config.yml config/library.config.yml

# Edit your library.config.yml
```

* Set up the configuration file (config.yml)
```sh
cp config/test/config.test.yml config/config.yml
# Edit your config.yml
```

### 3. Launch the workflow
* Activate the snakemake-slurm environment
```sh
conda activate snakemake-slurm
```

* Launch the workflow as follows, by setting appropriate number of jobs and cores depending on your libraries' replicates. We recommend as many jobs as the total number of replicates, and 4 cores per library. Note that these are the max resources that snakemake will be able to allocate as needed by the workflow steps; resources will be optimized to use the minimum number of jobs and cores as possible, depending on the required parallelization.
```sh
snakemake --executor slurm \
    --snakefile snakemake/Snakefile \
    --workflow-profile profile \
    --configfile config/config.yml \
    --jobs 12 --cores 16 --keep-going \
    --retries 2 --rerun-incomplete
```

### Test run
An example dataset is provided in the test folder. The set of configuration files can be found under config/test and profile/test. A test run can be launched as follows:
```sh
snakemake --executor slurm \
    --snakefile snakemake/Snakefile \
    --workflow-profile profile/test \
    --configfile config/test/config.test.yml \
    --jobs 6 --cores 12 --keep-going \
    --retries 2 --rerun-incomplete
```

[^scafe]: https://github.com/chung-lab/scafe
[^sala-wiki]: https://github.com/fantom-prj/SALA/wiki
[^conda]: https://www.anaconda.com/download
[^gdown]: https://github.com/wkentaro/gdown

## Rule diagram
![Snakemake rule graph](rulegraph.png)