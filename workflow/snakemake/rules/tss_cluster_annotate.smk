rule bam2ctss:
    input:
        f"{RESULTS}/copy_glm.done",
        lib_runs=lambda wc: LIBRARY_TO_BAM[wc.library]
    output:
        touch(f"{RESULTS}/ctss/{{library}}/.done")
    params:
        genome_anno_id=config["genome_anno_id"],
        scafe_path=SCAFE_PATH,
        tss_mode=PARAMS["bamtoctss_TSS_mode"],
        unencoded_G_upstrm_nt=PARAMS["bamtoctss_unencoded_G_upstrm_nt"],
        max_softclip_length=PARAMS["bamtoctss_max_softclip_length"]
    log:
        f"{LOGS}/bam2ctss/{{library}}.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    threads: 5
    shell:
        """
        scafe_bam2ctss="$CONDA_PREFIX/{params.scafe_path}/scripts/scafe.tool.bk.bam_to_ctss"
        perl {OTHER_SCRIPTS}/SCAFE.step1.bamtoctss.pl \
            --run_scafe_script_path=$scafe_bam2ctss \
            --in_lib_list_path={input.lib_runs} \
            --outDir={RESULTS}/ctss/{wildcards.library} \
            --genome={params.genome_anno_id} \
            --max_thread=5 \
            --TSS_mode={params.tss_mode} \
            --unencoded_G_upstrm_nt={params.unencoded_G_upstrm_nt} \
            --max_softclip_length={params.max_softclip_length} \
            > {log} 2>&1
        """

rule ctss_file_list:
    input:
        f"{RESULTS}/ctss/{{library}}/.done",
        lib_runs=lambda wc: LIBRARY_TO_BAM[wc.library]
    output:
        ctss_list=f"{RESULTS}/ctss/{{library}}/CTSS.list.txt"
    log:
        f"{LOGS}/list_files/{{library}}.log"
    run:
        with open(log[0], "w") as logfile:
            with redirect_stdout(logfile), redirect_stderr(logfile):
                with open(input.lib_runs) as f:
                    lib_samples =  [line.strip().split("\t")[1] for line in f if line.strip()]
                with open(output.ctss_list, "w") as out:
                    for sample in lib_samples:
                        collapse_ctss_file=f"{RESULTS}/ctss/{wildcards.library}/{sample}/bed/{sample}.collapse.ctss.bed.gz"
                        unencoded_G_ctss_file=f"{RESULTS}/ctss/{wildcards.library}/{sample}/bed/{sample}.unencoded_G.collapse.ctss.bed.gz"
                        out.write(f"{sample}\t{collapse_ctss_file}\t{unencoded_G_ctss_file}\n")

rule ctss_aggregate:
    input:
        ctss_list=f"{RESULTS}/ctss/{{library}}/CTSS.list.txt"
    output:
        f"{RESULTS}/aggregate/{{library}}/bed/{{library}}.aggregate.unencoded_G.collapse.ctss.bed.gz",
        f"{RESULTS}/5pr_data/{{library}}/5pr_signal.bed.gz",
        ctss_agg_file=f"{RESULTS}/aggregate/{{library}}/bed/{{library}}.aggregate.collapse.ctss.bed.gz"
    params:
        genome_anno_id=config["genome_anno_id"],
        signal_5pr_file=f"{RESULTS}/5pr_data/{{library}}/5pr_signal.bed.gz"
    log:
        f"{LOGS}/ctss_aggregate/{{library}}.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    threads: 5
    shell:
        """
        scafe.tool.cm.aggregate \
        --lib_list_path={input.ctss_list} \
        --max_thread=5 \
        --genome={params.genome_anno_id} \
        --outputPrefix={wildcards.library} \
        --outDir={RESULTS}/aggregate \
        > {log} 2>&1

        cp {output.ctss_agg_file} {params.signal_5pr_file}
        """

rule bed_to_bigwig:
    input:
        ctss_aggregate=f"{RESULTS}/aggregate/{{library}}/bed/{{library}}.aggregate.collapse.ctss.bed.gz",
        unencoded_G_aggregate=f"{RESULTS}/aggregate/{{library}}/bed/{{library}}.aggregate.unencoded_G.collapse.ctss.bed.gz"
    output:
        f"{RESULTS}/aggregate_bigwig/{{library}}.all/wig/{{library}}.all.count.fwd.bw",
        f"{RESULTS}/aggregate_bigwig/{{library}}.all/wig/{{library}}.all.count.rev.bw",
        f"{RESULTS}/aggregate_bigwig/{{library}}.ung/wig/{{library}}.ung.count.fwd.bw",
        f"{RESULTS}/aggregate_bigwig/{{library}}.ung/wig/{{library}}.ung.count.rev.bw"
    params:
        genome_anno_id=config["genome_anno_id"]
    log:
        f"{LOGS}/aggregate/{{library}}_tobw.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    shell:
        """
        (
        scafe.tool.cm.ctss_to_bigwig \
        --genome={params.genome_anno_id} \
        --ctss_bed_path={input.ctss_aggregate} \
        --outputPrefix={wildcards.library}.all \
        --outDir={RESULTS}/aggregate_bigwig

        scafe.tool.cm.ctss_to_bigwig \
        --genome={params.genome_anno_id} \
        --ctss_bed_path={input.unencoded_G_aggregate} \
        --outputPrefix={wildcards.library}.ung \
        --outDir={RESULTS}/aggregate_bigwig
        ) > {log} 2>&1
        """

rule ctss_cluster:
    input:
        ctss_aggregate=f"{RESULTS}/aggregate/{{library}}/bed/{{library}}.aggregate.collapse.ctss.bed.gz",
        unencoded_G_aggregate=f"{RESULTS}/aggregate/{{library}}/bed/{{library}}.aggregate.unencoded_G.collapse.ctss.bed.gz"
    output:
        f"{RESULTS}/cluster/{{library}}/bed/{{library}}.tssCluster.bed.gz"
    params:
        min_summit_count=PARAMS["ctss_cluster_min_summit_count"],
        min_cluster_count=PARAMS["ctss_cluster_min_cluster_count"]
    log:
        f"{LOGS}/ctss_cluster/{{library}}.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    shell:
        """
        scafe.tool.cm.cluster \
        --overwrite=yes \
        --cluster_ctss_bed_path={input.ctss_aggregate} \
        --count_ctss_bed_path={input.unencoded_G_aggregate} \
        --min_summit_count={params.min_summit_count} \
        --min_cluster_count={params.min_cluster_count} \
        --outputPrefix={wildcards.library} \
        --outDir={RESULTS}/cluster \
        > {log} 2>&1
        """
        
rule ctss_cluster_filter:
    input:
        ctss_aggregate=f"{RESULTS}/aggregate/{{library}}/bed/{{library}}.aggregate.collapse.ctss.bed.gz",
        unencoded_G_aggregate=f"{RESULTS}/aggregate/{{library}}/bed/{{library}}.aggregate.unencoded_G.collapse.ctss.bed.gz",
        ctss_cluster=f"{RESULTS}/cluster/{{library}}/bed/{{library}}.tssCluster.bed.gz"
    output:
        f"{RESULTS}/filter/{{library}}/bed/{{library}}.tssCluster.default.filtered.bed.gz",
        f"{RESULTS}/filter/{{library}}/log/{{library}}.tssCluster.log.tsv"
    params:
        genome_anno_id=config["genome_anno_id"]
    log:
        f"{LOGS}/ctss_cluster_filter/{{library}}.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    shell: 
        """
        scafe.tool.cm.filter \
        --overwrite=yes \
        --ctss_bed_path={input.ctss_aggregate} \
        --ung_ctss_bed_path={input.unencoded_G_aggregate} \
        --tssCluster_bed_path={input.ctss_cluster} \
        --genome={params.genome_anno_id} \
        --outputPrefix={wildcards.library} \
        --outDir={RESULTS}/filter \
        > {log} 2>&1
        """

rule ctss_annotate_cre:
    input:
        ctss_f_cluster=f"{RESULTS}/filter/{{library}}/bed/{{library}}.tssCluster.default.filtered.bed.gz",
        ctss_f_cluster_log=f"{RESULTS}/filter/{{library}}/log/{{library}}.tssCluster.log.tsv"
    output:
        f"{RESULTS}/annotate/{{library}}/bed/{{library}}.CRE.coord.bed.gz",
        f"{RESULTS}/annotate/{{library}}/log/{{library}}.CRE.info.tsv.gz",
        f"{RESULTS}/5pr_data/{{library}}/5pr_cluster.bed.gz",
        ctss_cluster_file=f"{RESULTS}/annotate/{{library}}/bed/{{library}}.cluster.coord.bed.gz"
    params:
        genome_anno_id=config["genome_anno_id"],
        cluster_5pr_file=f"{RESULTS}/5pr_data/{{library}}/5pr_cluster.bed.gz",
        min_cre_count=PARAMS["ctss_annotate_cre_min_cre_count"]
    log:
        f"{LOGS}/ctss_annotate_cre/{{library}}.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    shell:
        """
        scafe.tool.cm.annotate \
        --overwrite=yes \
        --tssCluster_bed_path={input.ctss_f_cluster} \
        --tssCluster_info_path={input.ctss_f_cluster_log} \
        --min_CRE_count={params.min_cre_count} \
        --genome={params.genome_anno_id} \
        --outputPrefix={wildcards.library} \
        --outDir={RESULTS}/annotate \
        > {log} 2>&1

        cp {output.ctss_cluster_file} {params.cluster_5pr_file}
        """

rule ctss_cre_directionality:
    input:
        cre_coord=f"{RESULTS}/annotate/{{library}}/bed/{{library}}.CRE.coord.bed.gz",
        cre_info=f"{RESULTS}/annotate/{{library}}/log/{{library}}.CRE.info.tsv.gz",
        ctss=f"{RESULTS}/aggregate/{{library}}/bed/{{library}}.aggregate.collapse.ctss.bed.gz"
    output:
        f"{RESULTS}/directionality/{{library}}/bed/{{library}}.CRE.directionality.summit_center.bed.gz",
        f"{RESULTS}/directionality/{{library}}/log/{{library}}.directionality.log.tsv.gz"
    log:
        f"{LOGS}/cre_directionality/{{library}}.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    shell:
        """
        scafe.tool.cm.directionality \
        --overwrite=yes \
        --CRE_bed_path={input.cre_coord} \
        --CRE_info_path={input.cre_info} \
        --ctss_bed_path={input.ctss} \
        --outputPrefix={wildcards.library} \
        --outDir={RESULTS}/directionality \
        > {log} 2>&1
        """

