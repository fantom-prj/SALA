rule bam_to_bed:
    input:
        lib_runs=lambda wc: LIBRARY_TO_BAM[wc.library]
    output:
        f"{RESULTS}/bamtobed/{{library}}/combined.bed.bgz"
    log:
        f"{LOGS}/bamtobed/{{library}}.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    threads: 4
    shell:
        """
        {OTHER_SCRIPTS}/SALA.input.bamtobed.sh \
        {input.lib_runs} \
        {RESULTS}/bamtobed/{wildcards.library} \
        > {log} 2>&1
        """

rule extract_ctes_signal:
    input:
        combined_bed=f"{RESULTS}/bamtobed/{{library}}/combined.bed.bgz",
        chrom_sizes=f"{REFERENCE}/chrom.sizes_major.tsv"
    output:
        f"{RESULTS}/ctes_signal/{{library}}/{{library}}.end3.bed",
        f"{RESULTS}/ctes_signal/{{library}}/{{library}}.end5.bed"          
    log:
        f"{LOGS}/ctes_signal/{{library}}.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    threads: 10
    shell:
        """
        perl {SALA_SCRIPTS}/3n_cluster/transcript_bed_to_end_bed_bigwig.pl \
        {input.combined_bed} \
        {input.chrom_sizes} \
        {wildcards.library} \
        {RESULTS}/ctes_signal/{wildcards.library} \
        > {log} 2>&1
        """

rule tes_cluster:
    input:
        tes_signal=f"{RESULTS}/ctes_signal/{{library}}/{{library}}.end3.bed"
    output:
        f"{RESULTS}/ctes_cluster/{{library}}/bed/{{library}}.tssCluster.bed.gz"
    params:
        min_summit_count=PARAMS["tes_cluster_min_summit_count"],
        min_nt_count=PARAMS["tes_cluster_min_nt_count"],
        min_cluster_count=PARAMS["tes_cluster_min_cluster_count"]
    log:
        f"{LOGS}/ctes_cluster/{{library}}.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    shell: 
        """
        scafe.tool.cm.cluster \
        --overwrite=yes \
        --cluster_ctss_bed_path={input.tes_signal} \
        --count_ctss_bed_path={input.tes_signal} \
        --min_summit_count={params.min_summit_count} \
        --min_nt_count={params.min_nt_count} \
        --min_cluster_count={params.min_cluster_count} \
        --outputPrefix={wildcards.library} \
        --outDir={RESULTS}/ctes_cluster \
        > {log} 2>&1
        """
        