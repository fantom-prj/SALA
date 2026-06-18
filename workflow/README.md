# SALA SLURM-Snakemake workflow

This workflow integrates core steps of SALA in a SLURM-ready snakemake pipeline to assemble a transcript and gene model starting from long-read alignments in BAM format. In particular, it goes through the steps 2, 3 and 4 of the SALA wiki[^sala-wiki]. 
This workflow allows to process independent libraries in parallel, through separated jobs orchestrated by snakemake. For each library, independent sub-jobs can be run in parallel as well. This is also handled by snakemake, depending on the available resources on your cluster.

## Dependencies
The workflow is designed to run on a SLURM-based HPC environment. The only prerequisite is a Conda-compatible package manager[^conda]. We recommend using Conda ≥ 23.0. All required channels and package dependencies are specified in the workflow environment files (see below).

## How to run

### 1. Environment setup
* Make sure SALA scripts are executable
* for snakemake pipeline, go to workflow directory and create snakemake environment with conda
```sh
chmod 755 -R ./code/
cd /your/SALA/directory/workflow
conda env create -f snakemake/environment.yml
```

### 2. Test run
* Activate the snakemake-slurm environment, download fasta / gtf for test run to test/resources
* Run the demo with the name of your slurm partition, you may identify it with sinfo

```sh
cd /your/SALA/directory/workflow/test/resources
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
gunzip GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz

```

* Update partition name of SLURM [slurm_partition: "cpuq" -> your partition name], you can identify it by sinfo
```sh
nano /your/SALA/directory/workflow/profile/test/config.yaml
```

* The set of configuration files can be found under config/test and profile/test.
* Test results located in test/results
```sh
conda activate snakemake-slurm

cd /your/SALA/directory/workflow

snakemake --executor slurm \
    --snakefile snakemake/Snakefile \
    --workflow-profile profile/test \
    --configfile config/test/config.test.yml \
    --jobs 6 --cores 12 --keep-going \
    --retries 2 --rerun-incomplete
```

### 3. Run your sample 
#### File setup
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

* Update partition name of SLURM [slurm_partition: "cpuq" -> your partition name], you can identify it by sinfo
```sh
nano /your/SALA/directory/workflow/profile/config.yaml
```

* Launch the workflow as follows, by setting appropriate number of jobs and cores depending on your libraries' replicates. We recommend as many jobs as the total number of replicates, and 4 cores per library. Note that these are the max resources that snakemake will be able to allocate as needed by the workflow steps; resources will be optimized to use the minimum number of jobs and cores as possible, depending on the required parallelization.
```sh
conda activate snakemake-slurm
cd /your/SALA/directory/workflow

snakemake --executor slurm \
    --snakefile snakemake/Snakefile \
    --workflow-profile profile \
    --configfile config/config.yml \
    --jobs 12 --cores 16 --keep-going \
    --retries 2 --rerun-incomplete
```

## Configuration

### General config

The configuration files needed for the workflow are found in the config folder. the ```config.yaml``` file is the entrypoint for all information regarding the location of input, results, library configuration and parameters. 

```yaml
samplesheet: config/test/sample_table_all.csv [Path of the library table (see below), *need update]
results_dir: test/results [Path to the folder where results will be stored, *need update]
reference_dir: test/reference [Path where processed files from reference annotation will be stored, *need update]
logs_dir: test/logs [Path to the folder where logs will be stored, *need update]
fasta_genome: test/resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta [Path of the reference genome (.fasta). It needs to be unzipped and the same the BAM reads are aligned to]
mask_bed: test/resources/GRCh38-cCREs.sorted.bed [For SCAFE alone Path of the promoter/enhancer annotation on the reference genome (.bed)]
gtf_annotation: test/resources/gencode.v47.annotation.gtf.gz [Path of the reference transcriptome annotation (.gtf.gz)]
genome_anno_id: hg38_gencode47 [Id for the genome+annotation folder (for SCAFE and SALA internal usage)]
genome_code: hg38 [Code name for the genome (for SALA internal usage)]
library_config: config/test/library.config.yml [Path of the library configuration file (see below), *need update]
params_set_file: config/params.csv [Path of the parameters set file]
params_set: default [Choice of parameters set (must be column in params_set_file), can be customized by adding column in config/params.csv]
```

The ```params.csv``` file contains pre-defined parameters set (default and sensitive) for the various steps of the workflow. To define your own set of parameters, it is possible to add an additional column to the csv, and to change the parameters set choice accordingly in the config.yml.

##Parameter description

| Parameter | Default | Description |
|------------|---------|-------------|
| `bamtoctss_TSS_mode` | `softclip` | Method used to define the transcription start site (TSS). `softclip`: TSS is defined by the 5′ soft-clipped position. `first_match`: TSS is defined by the first nucleotide aligned to the genome. |
| `bamtoctss_unencoded_G_upstrm_nt` | `3` | Number of upstream nucleotides examined for unencoded G residues. If set to `2`, a maximum of two unencoded Gs can be detected per read. |
| `bamtoctss_max_softclip_length` | `3` | Maximum allowed 5′ soft-clip length for a read to be considered a valid TSS. |
| `ctss_cluster_min_summit_count` | `0` | Minimum summit read count required for a CTSS cluster to be considered confident. |
| `ctss_cluster_min_cluster_count` | `1` | Minimum total read count required for a CTSS cluster to be considered confident. |
| `ctss_annotate_cre_min_cre_count` | `3` | Minimum read count required for a transcribed CRE (tCRE) to be considered confident. |
| `tes_cluster_min_summit_count` | `3` | Minimum summit read count required for a TES cluster to be considered confident. |
| `tes_cluster_min_cluster_count` | `5` | Minimum total read count required for a TES cluster to be considered confident. |
| `tes_cluster_min_nt_count` | `3` | Minimum cluster width (nt) required for a TES cluster to be considered confident. |
| `extract_junctions_min_nt_qual` | `10` | Minimum base quality required at all six splice junction flanking positions (-3,-2,-1,+1,+2,+3) for a junction to be considered high quality. |
| `extract_junctions_min_mapq` | `20` | Minimum read MAPQ required for splice junction extraction. |
| `assemble_tx_model_min_output_qry_count` | `1` | Output transcript models supported by at least this number of reads. |
| `assemble_tx_model_conf_end3_merge_flank` | `150` | Flanking distance (bp) used to merge nearby confident 3′ end clusters. Use `-1` to disable merging. |
| `assemble_tx_model_conf_end5_merge_flank` | `75` | Flanking distance (bp) used to merge nearby confident 5′ end clusters. Use `-1` to disable merging. |
| `assemble_tx_model_conf_end3_add_ref` | `yes` | Include annotated reference 3′ ends when constructing confident 3′ end regions. |
| `assemble_tx_model_conf_end5_add_ref` | `yes` | Include annotated reference 5′ ends when constructing confident 5′ end regions. |
| `assemble_tx_model_min_exon_length` | `1` | Minimum exon length required for a transcript model to be retained. |
| `assemble_tx_model_min_transcript_length` | `15` | Minimum transcript length required for a transcript model to be retained. |
| `assemble_tx_model_print_trnscrptID` | `no` | Report transcript IDs in the output. |
| `assemble_tx_model_trnscpt_set_end_priority` | `summit:commonest:longest` | Priority used to determine transcript set boundaries. `summit`: strongest end signal. `commonest`: most frequently observed end. `longest`: most distal observed end. |
| `assemble_tx_model_doubtful_end_merge_dist` | `150` | Distance used to group incomplete transcript ends. |
| `assemble_tx_model_doubtful_end_avoid_summit` | `yes` | Ignore summit-based end selection when resolving doubtful transcript ends. |
| `assemble_tx_model_min_summit_dist_split` | `50` | Minimum distance between two summits when splitting an end cluster. |
| `assemble_tx_model_retain_no_qry_ref_bound_set` | `no` | Retain transcript boundary sets supported only by the reference annotation. |
| `assemble_tx_model_min_size_split` | `100` | Minimum cluster size required to split an end cluster. |
| `assemble_tx_model_min_frac_split` | `0.2` | Minimum signal fraction required at both summits when splitting an end cluster. |
| `assemble_tx_model_min_qry_score` | `0` | Minimum query score (typically MAPQ) required for transcript assembly. |
| `assemble_gene_model_min_ref_exon_overlap_pct` | `10` | Minimum exon overlap percentage required to assign transcripts to the same gene. |
| `assemble_gene_model_exon_overlap_dist` | `-1` | Minimum exon overlap distance. Can be used together with `assemble_gene_model_min_ref_exon_overlap_pct`. |
| `assemble_gene_model_locus_merge_dist` | `100000` | Genomic search window used during gene model assembly. |
| `assemble_gene_model_exclude_t_type` | `retained_intron` | Transcript types excluded from gene model assembly. |
| `filter_tx_model_read_per_rep_ref_novel_tx` | `5` | Minimum full-length reads per replicate required for novel isoforms of annotated genes. |
| `filter_tx_model_read_per_rep_noref_novel_tx` | `1` | Minimum full-length reads per replicate required for transcripts from novel genes. |
| `filter_tx_model_isoform_ratio` | `0.1` | Minimum isoform abundance ratio required for novel isoforms of annotated genes. |
| `filter_tx_model_require_5pr_confidence` | `Yes` | If enabled, transcript models without supported 5′ ends are removed. |
| `filter_tx_model_n3_confid_ref_novel_tx` | `Yes` | If enabled, novel isoforms of annotated genes without supported 3′ ends are removed. |
| `filter_tx_model_n3_confid_noref_novel_tx` | `No` | If enabled, transcripts from novel genes without supported 3′ ends are removed. |



The ```profile/config.yml``` file is related to resources configuration. Here it is possible to specify appropriate RAM requirements for the various steps. To handle out of memory issues, by default, the workflow will perform two retries for each step, each time doubling the allocated memory. The mem_mb default values have been defined according to tests on libraries with 60 to 250 million reads. We recommend to change them accordingly to your library size needs.

### Data config
The ```sample_table.csv``` contains the list of all libraries to be processed in parallel, and has to be provided in the sample table specified in the configuration (samplesheet), in this format (header included):
```csv
library_id,bam_list,bam_tc_list,rep_list
neuronset_id,/path/of/neuset/run/file.tsv,/path/of/neuset/run_tc/file.tsv,/path/of/neuset/replicates/file.tsv
thp1_id,/path/of/thp1/run/file.tsv,/path/of/thp1/run_tc/file.tsv,/path/of/thp1/replicates/file.tsv
```
The second column, ```bam_list```, must contain the path of the runs file of the BAMs prior to running TranscriptClean[^tc], that are needed specifically for splicing junctions extraction. The third column, ```bam_tc_list```, can optionally contain a path for the runs file for the splicing-corrected BAMs (if TranscriptClean has been run for the library), or be left empty otherwise. In this case, the run file in bam_list will be used for all the steps. You can refer to the SALA wiki for more details about running TranscriptClean[^sala-wiki-tc] and splicing junction extraction[^sala-wiki-sj].

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
As explained in the sample table specification, in case you want to use TranscriptClean BAMs, two run files in this format must be produced: one with the paths of the pre-TranscriptClean BAMs and one with the paths of the post-TranscriptClean BAMs. For more details, refer to the SALA wiki.

```library.config.yml```
In this file, the following information about the libraries have to be specified:
* `internal_priming`: true to remove internal priming affected transcripts due to library preparation (poly-dT priming), false otherwise.
* `5pr_confident`: true if the reads have confident 5' ends, false otherwise. In this latter case, external TSS clusters need to be provided as a .bed.gz file. Pre-set cluster coordinates for human (hg38) and mouse (mm10) can be found in the workflow/resources/CTSS_cluster folder.
* `5pr_cluster`: location of the externally provided TSS clusters file. Mandatory if 5pr_confident is false.

An example of a config file for a not confident 5' library can be found under config/test.
```library_notconf.config.yml
internal_priming: false
5pr_confident: false
5pr_cluster: /your/SALA/directory/resources/CTSS_cluster/hg38_fair+new_CAGE_peaks_phase1and2.bed.bgz
```
5pr_cluster files (Human and mouse TSS regions) are provided in /your/SALA/directory/resources/CTSS_cluster


Specifically for SCAFE prepare genome, besides providing fasta files, cCRE annotation (mask_bed) downloaded from SCREEN should be applied.
SCREEN v3 human and mouse cCRE files were stored in /your/SALA/directory/workflow/test/resources/


[^sala-wiki]: https://github.com/fantom-prj/SALA/wiki
[^tc]: https://github.com/mortazavilab/TranscriptClean
[^sala-wiki-tc]: https://github.com/fantom-prj/SALA/wiki/1.-TranscriptClean
[^sala-wiki-sj]: https://github.com/fantom-prj/SALA/wiki/3.-Features-preparation#34-splicing-junction-collection
[^conda]: https://www.anaconda.com/download
[^gdown]: https://github.com/wkentaro/gdown

## Rule diagram
![Snakemake rule graph](rulegraph.png) <!--TODO generate again -->