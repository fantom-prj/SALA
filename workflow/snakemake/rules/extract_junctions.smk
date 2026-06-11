rule extract_junctions:
    input:
        chrom_sizes=CHR_SIZES,
        bam_file=lambda wc: BAM_TO_PATH[wc.bam]
    output:
        f"{RESULTS}/junction_extractor/tmp/{{bam}}/log/{{bam}}.junct.info.tsv.gz"
    params:
        genome=config["fasta_genome"],
        min_nt_qual=PARAMS["extract_junctions_min_nt_qual"],
        min_mapq=PARAMS["extract_junctions_min_mapq"] 
    log:
        f"{LOGS}/extract_junctions/{{bam}}.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml" 
    shell:
        """
        perl {SALA_SCRIPTS}/junction_extractor/junction_extractor.pl \
        --in_bam={input.bam_file} \
        --chrom_size_path={input.chrom_sizes} \
        --chrom_fasta_path={params.genome} \
        --out_prefix={wildcards.bam} \
        --out_dir={RESULTS}/junction_extractor/tmp \
        --max_thread=1 \
        --min_nt_qual={params.min_nt_qual} \
        --min_MAPQ={params.min_mapq} \
        > {log} 2>&1
        """

rule group_library_bams:
    input:
        lambda wc: [
            f"{RESULTS}/junction_extractor/tmp/{bam[0]}/log/{bam[0]}.junct.info.tsv.gz"
            for bam in LIB_BAMS[wc.library]
        ],
        lib_runs=lambda wc: LIBRARY_TO_BAM[wc.library]
    output:
        touch(f"{RESULTS}/junction_extractor/{{library}}_group.done")      
    params:
        lib_dir=f"{RESULTS}/junction_extractor/{{library}}",
        tmp_dir=f"{RESULTS}/junction_extractor/tmp"
    log:
        f"{LOGS}/extract_junctions/{{library}}.log"
    run:
        with open(log[0], "w") as logfile:
            with redirect_stdout(logfile), redirect_stderr(logfile):
                with open(input.lib_runs) as f:
                    lib_samples =  [line.strip().split("\t")[1] for line in f if line.strip()]
                os.makedirs(params.lib_dir, exist_ok=True)
                for sample in lib_samples:
                    sample_file=f"{sample}.junct.info.tsv.gz"
                    move_from=os.path.join(params.tmp_dir, sample, "log", sample_file)
                    move_to=os.path.join(params.lib_dir, sample_file)
                    os.rename(move_from, move_to)

rule pool_junctions:
    input:
        f"{RESULTS}/junction_extractor/{{library}}_group.done",
    output:
        f"{RESULTS}/junction_extractor/{{library}}/{{library}}.hi_qual.junct.bed"
    log:
        f"{LOGS}/extract_junctions/{{library}}.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml" 
    shell:
        """
        {SALA_SCRIPTS}/junction_extractor/junction_pool.pl \
        --outDir={RESULTS}/junction_extractor/{wildcards.library} \
        --outTag={wildcards.library} \
        --findStr={RESULTS}/junction_extractor/{wildcards.library}/*.junct.info.tsv.gz \
        > {log} 2>&1
        """