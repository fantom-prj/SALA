rule prepare_genome:
    input:
        gtf=config["gtf_annotation"],
        fasta=config["fasta_genome"],
        chr_sizes=f"{REFERENCE}/chrom.sizes_major.tsv",
        mask=config["mask_bed"]
    output:
        touch(f"{RESULTS}/prep_genome.done")
    params:
        genome_anno_id=config["genome_anno_id"],
        scafe_path=SCAFE_PATH
    log:
        f"{LOGS}/prepare_genome.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    shell:
        """
        outdir=$CONDA_PREFIX/{params.scafe_path}/resources/genome
        scafe.tool.cm.prep_genome \
            --gtf_path={input.gtf} \
            --fasta_path={input.fasta} \
            --chrom_list_path={input.chr_sizes} \
            --mask_bed_path={input.mask} \
            --outputPrefix={params.genome_anno_id} \
            --outDir=$outdir \
            > {log} 2>&1
        """

rule copy_glm_folder:
    input:
        f"{RESULTS}/prep_genome.done"
    output:
        touch(f"{RESULTS}/copy_glm.done")
    params:
        scafe_path=SCAFE_PATH,
        genome_anno_id=config["genome_anno_id"]
    log:
        f"{LOGS}/copy_glm_folder.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    shell:
        """
        copy_from={SCAFE_RESOURCES}
        copy_to=$CONDA_PREFIX/{params.scafe_path}/resources/genome/{params.genome_anno_id}
        cp -r $copy_from/glm $copy_to \
        > {log} 2>&1
        """
