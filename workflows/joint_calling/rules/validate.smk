rule gatk_concordance_NA12878:
    input:
        ref = rules.download_reference_genome.output,
        intervals = config["exon_bed"],
        truth_vcf = config["NA12878_validation_vcf"],
        truth_idx = config["NA12878_validation_vcf"] + ".tbi",
        eval_vcf = rules.merge_snp_indel_vcf.output.vcf,
        eval_idx = rules.merge_snp_indel_vcf.output.idx
    output:
        os.path.join(out_path, "concordance_NA12878/gatk/{group}.tsv")
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(out_path, "log/concordance_NA12878/gatk/{group}.log")
    shell:
        """
        gatk Concordance \
            -R {input.ref} \
            -eval {input.eval_vcf} \
            --truth {input.truth_vcf} \
            --intervals {input.intervals} \
            --summary {output} 2> {log}
        """

rule snpsift_concordance_NA12878:
    input:
        truth_vcf = config["NA12878_validation_vcf"],
        truth_idx = config["NA12878_validation_vcf"] + ".tbi",
        eval_vcf = rules.merge_snp_indel_vcf.output.vcf,
        eval_idx = rules.merge_snp_indel_vcf.output.idx
    output:
        truth_vcf = os.path.join(out_path, "concordance_NA12878/snpsift/truth_{group}.vcf"),
        eval_vcf = os.path.join(out_path, "concordance_NA12878/snpsift/eval_{group}.vcf"),
        summary = os.path.join(out_path, "concordance_NA12878/snpsift/{group}.tsv")
    conda:
        "../envs/snpsift.yaml"
    log:
        os.path.join(out_path, "log/concordance_NA12878/snpsift/{group}.log")
    shell:
        """
        zcat {input.truth_vcf} > {output.truth_vcf}
        zcat {input.eval_vcf} > {output.eval_vcf}

        SnpSift concordance \
            -v {output.truth_vcf} \
            {output.eval_vcf} > {output.summary} 2> {log}
        """
