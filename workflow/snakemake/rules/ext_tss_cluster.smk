rule handle_ext_tss_cluster:
    input: 
        f"{RESULTS}/copy_glm.done",
        ext_ctss_cluster=lib_config["5pr_cluster"],
        signal_5pr=f"{RESULTS}/ctes_signal/{{library}}/{{library}}.end5.bed"
    output:
        f"{RESULTS}/5pr_data/{{library}}/5pr_cluster.bed.gz",
        f"{RESULTS}/5pr_data/{{library}}/5pr_signal.bed.gz"
    params:
        cluster_5pr_file=f"{RESULTS}/5pr_data/{{library}}/5pr_cluster.bed.gz",
        signal_5pr_file=f"{RESULTS}/5pr_data/{{library}}/5pr_signal.bed.gz",
    log:
        f"{LOGS}/ext_ctss_cluster/{{library}}.log"
    shell:
        """
            gzipped_signal={input.signal_5pr}.gz
            gzip -c {input.signal_5pr} > $gzipped_signal
            cp {input.ext_ctss_cluster} {params.cluster_5pr_file}
            cp $gzipped_signal {params.signal_5pr_file}
        """
