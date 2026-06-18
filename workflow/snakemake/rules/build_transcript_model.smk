rule index_bed_files:
    input:
        cluster_5pr=f"{RESULTS}/5pr_data/{{library}}/5pr_cluster.bed.gz",
        signal_5pr=f"{RESULTS}/5pr_data/{{library}}/5pr_signal.bed.gz",
        cluster_3pr=f"{RESULTS}/ctes_cluster/{{library}}/bed/{{library}}.tssCluster.bed.gz",
        signal_3pr=f"{RESULTS}/ctes_signal/{{library}}/{{library}}.end3.bed"
    output:
        f"{RESULTS}/5pr_data/{{library}}/5pr_cluster.bed.bgz",
        f"{RESULTS}/5pr_data/{{library}}/5pr_cluster.bed.bgz.tbi",
        f"{RESULTS}/5pr_data/{{library}}/5pr_signal.bed.bgz",
        f"{RESULTS}/5pr_data/{{library}}/5pr_signal.bed.bgz.tbi",
        f"{RESULTS}/ctes_cluster/{{library}}/bed/{{library}}.tssCluster.bed.bgz",
        f"{RESULTS}/ctes_cluster/{{library}}/bed/{{library}}.tssCluster.bed.bgz.tbi",
        f"{RESULTS}/ctes_signal/{{library}}/{{library}}.end3.bed.bgz",
        f"{RESULTS}/ctes_signal/{{library}}/{{library}}.end3.bed.bgz.tbi"
    log:
        f"{LOGS}/index_bed_{{library}}.log"  
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml" 
    shell:
        """
        {OTHER_SCRIPTS}/SALA.input.bed_index.sh \
        {RESULTS} \
        {wildcards.library} \
        {input.cluster_5pr} \
        {input.signal_5pr} \
        {input.cluster_3pr} \
        {input.signal_3pr} \
        > {log} 2>&1
        """

rule assemble_transcript_model:
    input:
        chr_sizes=CHR_SIZES,
        reads_bed=f"{RESULTS}/bamtobed/{{library}}/combined.bed.bgz",
        transcr_ref=f"{REFERENCE}/transcript.bed.bgz",
        end5_cluster=f"{RESULTS}/5pr_data/{{library}}/5pr_cluster.bed.bgz",
        end3_cluster=f"{RESULTS}/ctes_cluster/{{library}}/bed/{{library}}.tssCluster.bed.bgz",
        end5_signal=f"{RESULTS}/5pr_data/{{library}}/5pr_signal.bed.bgz",
        end3_signal=f"{RESULTS}/ctes_signal/{{library}}/{{library}}.end3.bed.bgz",
        conf_junct=f"{RESULTS}/junction_extractor/{{library}}/{{library}}.hi_qual.junct.bed",
        ref_junct=f"{REFERENCE}/junction.bed"
    output:
        f"{RESULTS}/transcript_model/{{library}}/bed/{{library}}.model.bed.bgz",
        f"{RESULTS}/transcript_model/{{library}}/log/{{library}}.model.info.tsv.gz"
    params:
        genome=config["fasta_genome"],
        min_output_qry_count=PARAMS["assemble_tx_model_min_output_qry_count"],
        conf_end3_merge_flank=PARAMS["assemble_tx_model_conf_end3_merge_flank"],
        conf_end5_merge_flank=PARAMS["assemble_tx_model_conf_end5_merge_flank"],
        conf_end3_add_ref=PARAMS["assemble_tx_model_conf_end3_add_ref"],
        conf_end5_add_ref=PARAMS["assemble_tx_model_conf_end5_add_ref"],
        min_exon_length=PARAMS["assemble_tx_model_min_exon_length"],
        min_transcript_length=PARAMS["assemble_tx_model_min_transcript_length"],
        print_trnscrptID=PARAMS["assemble_tx_model_print_trnscrptID"],
        trnscpt_set_end_priority=PARAMS["assemble_tx_model_trnscpt_set_end_priority"],
        doubtful_end_merge_dist=PARAMS["assemble_tx_model_doubtful_end_merge_dist"],
        doubtful_end_avoid_summit=PARAMS["assemble_tx_model_doubtful_end_avoid_summit"],
        min_summit_dist_split=PARAMS["assemble_tx_model_min_summit_dist_split"],
        retain_no_qry_ref_bound_set=PARAMS["assemble_tx_model_retain_no_qry_ref_bound_set"],
        min_size_split=PARAMS["assemble_tx_model_min_size_split"],
        min_frac_split=PARAMS["assemble_tx_model_min_frac_split"],
        min_qry_score=PARAMS["assemble_tx_model_min_qry_score"]
    log:
        f"{LOGS}/assemble_tx_model_{{library}}.log"  
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    threads: 2
    shell: 
        """
        perl {SALA_SCRIPTS}/end5_guided_assembler.pl \
        --qry_bed_bgz={input.reads_bed} \
        --ref_bed_bgz={input.transcr_ref} \
        --out_dir={RESULTS}/transcript_model \
        --chrom_size_path={input.chr_sizes} \
        --chrom_fasta_path={params.genome} \
        --conf_end5_bed_bgz={input.end5_cluster} \
        --conf_end3_bed_bgz={input.end3_cluster} \
        --signal_end5_bed_bgz={input.end5_signal} \
        --signal_end3_bed_bgz={input.end3_signal} \
        --out_prefix={wildcards.library} \
        --novel_model_prefix=NOVT \
        --min_output_qry_count={params.min_output_qry_count} \
        --conf_end3_merge_flank={params.conf_end3_merge_flank} \
        --conf_end5_merge_flank={params.conf_end5_merge_flank} \
        --conf_end3_add_ref={params.conf_end3_add_ref} \
        --conf_end5_add_ref={params.conf_end5_add_ref} \
        --min_exon_length={params.min_exon_length} \
        --min_transcript_length={params.min_transcript_length} \
        --print_trnscrptID={params.print_trnscrptID} \
        --trnscpt_set_end_priority={params.trnscpt_set_end_priority} \
        --doubtful_end_merge_dist={params.doubtful_end_merge_dist} \
        --doubtful_end_avoid_summit={params.doubtful_end_avoid_summit} \
        --min_summit_dist_split={params.min_summit_dist_split} \
        --retain_no_qry_ref_bound_set={params.retain_no_qry_ref_bound_set} \
        --min_size_split={params.min_size_split} \
        --min_frac_split={params.min_frac_split} \
        --max_thread=2 \
        --conf_junction_bed={input.ref_junct},{input.conf_junct} \
        --min_qry_score={params.min_qry_score} \
        > {log} 2>&1
        """

rule assemble_gene_model:
    input:
        chr_sizes=CHR_SIZES,
        transcr_model=f"{RESULTS}/transcript_model/{{library}}/bed/{{library}}.model.bed.bgz",
        transcr_model_info=f"{RESULTS}/transcript_model/{{library}}/log/{{library}}.model.info.tsv.gz",
        transcr_ref=f"{REFERENCE}/transcript.bed.bgz",
        transcr_to_gene=f"{REFERENCE}/transcript_to_gene.tsv"
    output:
        f"{RESULTS}/gene_model/{{library}}/bed/{{library}}.model.bed.bgz",
        f"{RESULTS}/gene_model/{{library}}/log/{{library}}.model.info.tsv.gz"
    params:
        min_ref_exon_overlap_pct=PARAMS["assemble_gene_model_min_ref_exon_overlap_pct"],
        exon_overlap_dist=PARAMS["assemble_gene_model_exon_overlap_dist"],
        locus_merge_dist=PARAMS["assemble_gene_model_locus_merge_dist"],
        exclude_t_type=PARAMS["assemble_gene_model_exclude_t_type"]
    log:
        f"{LOGS}/assemble_gene_model_{{library}}.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    threads: 2
    shell:
        """
        perl {SALA_SCRIPTS}/assemble_gene_annotator.pl \
        --chrom_size_path={input.chr_sizes} \
        --model_bed_bgz={input.transcr_model} \
        --model_info_gz={input.transcr_model_info} \
        --revert_ref_model_bed_bgz={input.transcr_ref} \
        --ref_model_gene_link={input.transcr_to_gene} \
        --out_prefix={wildcards.library} \
        --out_dir={RESULTS}/gene_model \
        --novel_gene_prefix=NOVG \
        --disable_ref_chain_bound_gene_anno=yes \
        --min_ref_exon_overlap_pct={params.min_ref_exon_overlap_pct} \
        --exon_overlap_dist={params.exon_overlap_dist} \
        --locus_merge_dist={params.locus_merge_dist} \
        --exclude_t_type={params.exclude_t_type} \
        --max_thread=2 \
        > {log} 2>&1
        """

rule compute_count_matrix:
    input:
        transcr_to_gene=f"{REFERENCE}/transcript_to_gene.tsv",
        transcr_model=f"{RESULTS}/transcript_model/{{library}}/bed/{{library}}.model.bed.bgz"
    output:
        f"{RESULTS}/transcript_model/{{library}}/log/full_length_support_matrix.tsv.gz",
        f"{RESULTS}/transcript_model/{{library}}/log/partial_length_support_matrix.tsv.gz"
    params:
        transcr_model_fd=f"{RESULTS}/transcript_model/{{library}}",
        out_dir=f"{RESULTS}/transcript_model/{{library}}/log"
    log:
        f"{LOGS}/count_matrix_{{library}}.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    shell:
        """
            Rscript {OTHER_SCRIPTS}/SALA.count_matrix.R \
            {params.transcr_model_fd} \
            {params.out_dir} \
            {input.transcr_to_gene} \
            > {log} 2>&1
        """

rule filter_transcript_model:
    input:
        genome=config["fasta_genome"],
        lib_rep_file=lambda wc: LIBRARY_TO_REP[wc.library],
        transcr_model_counts=f"{RESULTS}/transcript_model/{{library}}/log/full_length_support_matrix.tsv.gz",
        gene_model=f"{RESULTS}/gene_model/{{library}}/log/{{library}}.model.info.tsv.gz"
    output:
        f"{RESULTS}/transcript_model/{{library}}/bed/{{library}}.table4.bed12.bed.bgz",
        f"{RESULTS}/transcript_model/{{library}}/log/{{library}}.table4_filtered.All_Ref.tsv.gz"
    params:
        transcr_model_fd=f"{RESULTS}/transcript_model/{{library}}",
        gene_model_fd=f"{RESULTS}/gene_model/{{library}}",
        genome_code=config["genome_code"],
        keep_IP=(not lib_config["internal_priming"]),
        read_per_rep_ref_novel_tx=PARAMS["filter_tx_model_read_per_rep_ref_novel_tx"],
        read_per_rep_noref_novel_tx=PARAMS["filter_tx_model_read_per_rep_noref_novel_tx"],
        isoform_ratio=PARAMS["filter_tx_model_isoform_ratio"],
        require_5pr_confidence=PARAMS["filter_tx_model_require_5pr_confidence"],
        require_3pr_confidence_ref_novel_tx=PARAMS["filter_tx_model_n3_confid_ref_novel_tx"],
        require_3pr_confidence_noref_novel_tx=PARAMS["filter_tx_model_n3_confid_noref_novel_tx"],
        scafe_dir=RESULTS if lib_config["5pr_confident"] else "NA"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    resources:
        retry_log=lambda wildcards, attempt:
            f"{LOGS}/filter_tx_model_{wildcards.library}_at{attempt}.log"
    shell:
        """
            Rscript {OTHER_SCRIPTS}/SALA.filter.R \
            {params.transcr_model_fd} \
            {wildcards.library} \
            {RESOURCES} \
            {REFERENCE} \
            {input.genome} \
            {params.read_per_rep_ref_novel_tx} \
            {params.read_per_rep_noref_novel_tx} \
            {params.isoform_ratio} \
            {params.require_5pr_confidence} \
            {params.gene_model_fd} \
            {input.lib_rep_file} \
            {params.scafe_dir} \
            {params.genome_code} \
            {params.keep_IP} \
            {params.require_3pr_confidence_ref_novel_tx} \
            {params.require_3pr_confidence_noref_novel_tx} \
            {RESULTS}/transcript_model/{wildcards.library}/cpat \
            > {resources.retry_log} 2>&1
        """

rule post_filter_assemble_gene_model:
    input:
        chr_sizes=CHR_SIZES,
        transcr_model=f"{RESULTS}/transcript_model/{{library}}/bed/{{library}}.table4.bed12.bed.bgz",
        transcr_model_info=f"{RESULTS}/transcript_model/{{library}}/log/{{library}}.table4_filtered.All_Ref.tsv.gz", #here too
        transcr_ref=f"{REFERENCE}/transcript.bed.bgz",
        transcr_to_gene=f"{REFERENCE}/transcript_to_gene.tsv"
    output:
        f"{RESULTS}/gene_model_filtered/{{library}}/bed/{{library}}.model.bed.bgz",
        f"{RESULTS}/gene_model_filtered/{{library}}/log/{{library}}.model.info.tsv.gz"
    params:
        min_ref_exon_overlap_pct=PARAMS["assemble_gene_model_min_ref_exon_overlap_pct"],
        exon_overlap_dist=PARAMS["assemble_gene_model_exon_overlap_dist"],
        locus_merge_dist=PARAMS["assemble_gene_model_locus_merge_dist"],
        exclude_t_type=PARAMS["assemble_gene_model_exclude_t_type"]
    log:
        f"{LOGS}/assemble_gene_model_filtered_{{library}}.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    threads: 2
    shell:
        """
        perl {SALA_SCRIPTS}/assemble_gene_annotator.pl \
        --chrom_size_path={input.chr_sizes} \
        --model_bed_bgz={input.transcr_model} \
        --model_info_gz={input.transcr_model_info} \
        --revert_ref_model_bed_bgz={input.transcr_ref} \
        --ref_model_gene_link={input.transcr_to_gene} \
        --out_prefix={wildcards.library} \
        --out_dir={RESULTS}/gene_model_filtered \
        --novel_gene_prefix=NOVG \
        --disable_ref_chain_bound_gene_anno=no \
        --min_ref_exon_overlap_pct={params.min_ref_exon_overlap_pct} \
        --exon_overlap_dist={params.exon_overlap_dist} \
        --locus_merge_dist={params.locus_merge_dist} \
        --exclude_t_type={params.exclude_t_type} \
        --max_thread=2 \
        > {log} 2>&1
        """

rule build_gtf_annotation:
    input:
        f"{RESULTS}/copy_gtf.done",
        f"{RESULTS}/gene_model_filtered/{{library}}/bed/{{library}}.model.bed.bgz",
        f"{RESULTS}/gene_model_filtered/{{library}}/log/{{library}}.model.info.tsv.gz",
        transcr_model=f"{RESULTS}/transcript_model/{{library}}/bed/{{library}}.table4.bed12.bed.bgz",
    output:
        f"{RESULTS}/transcript_model/{{library}}/log/table4.All_Ref.gtf.gz",
        f"{RESULTS}/transcript_model/{{library}}/log/table4.Detected_Ref.gtf.gz"
    params:
        transcr_model_fd=f"{RESULTS}/transcript_model/{{library}}",
        gene_model_fd=f"{RESULTS}/gene_model_filtered/{{library}}"
    log:
        f"{LOGS}/build_gtf_annotation_{{library}}.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    shell:
        """
        Rscript {OTHER_SCRIPTS}/SALA.gene_gtf_annotation.R \
        {params.transcr_model_fd} \
        {wildcards.library} \
        {RESOURCES} \
        {REFERENCE} \
        {params.gene_model_fd} \
        > {log} 2>&1
        """