rule gtf_to_bed12:
    input:
        transcriptome_gtf=config["gtf_annotation"]
    output:
        f"{REFERENCE}/transcript.bed.gz",
        f"{REFERENCE}/transcript.bed.bgz"
    log:
        f"{LOGS}/tx_gtf_to_bed.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    shell:
        """
        (
        transcript_bed={REFERENCE}/transcript.bed.gz
        transcript_bed_bgz={REFERENCE}/transcript.bed.bgz
        gzip -d -c {input.transcriptome_gtf} | bedparse gtf2bed | gzip > $transcript_bed
        gzip -d -c $transcript_bed | sort -k1,1 -k2,2n | bgzip > $transcript_bed_bgz
        tabix -p bed $transcript_bed_bgz
        ) > {log} 2>&1
        """

rule copy_gtf:
    input:
        transcriptome_gtf=config["gtf_annotation"]
    output:
        touch(f"{RESULTS}/copy_gtf.done")
    log:
        f"{LOGS}/copy_gtf.log"
    shell:
        """
        cp {input.transcriptome_gtf} {REFERENCE} > {log} 2>&1
        """

rule extract_transcriptome_junctions:
    input:
        transcriptome_bed=f"{REFERENCE}/transcript.bed.gz"
    output:
        f"{REFERENCE}/junction.bed"
    log:
        f"{LOGS}/tx_extract_junctions.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    shell:
        """
        (
        junction_bed={REFERENCE}/junction.bed
        tmp={REFERENCE}/temp.bed
        gzip -d -c {input.transcriptome_bed} | bedparse introns | bed12ToBed6 -i stdin > $junction_bed
        sort -k1,1 -k2,2n -k3,3n -k6,6 $junction_bed | awk '!seen[$1, $2, $3, $6]++' > $tmp && mv $tmp $junction_bed
        ) > {log} 2>&1
        """

rule extract_transcript_to_gene:
    input:
        genome=config["fasta_genome"],
        transcriptome_gtf=config["gtf_annotation"]
    output:
        f"{REFERENCE}/transcript_to_gene.tsv",
        CHR_SIZES
    log:
        f"{LOGS}/tx_to_gene.log"
    conda:
        f"{workflow.basedir}/envs/sala_wf_env.yml"
    shell:
        """
           Rscript {OTHER_SCRIPTS}/SALA.gtf_extract.R \
            {input.transcriptome_gtf} \
            {REFERENCE} \
            {input.genome} \
            > {log} 2>&1
        """
