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
        ref = rules.download_reference_genome.output,
        intervals = config["exon_bed"],
        truth_vcf = config["NA12878_validation_vcf"],
        truth_idx = config["NA12878_validation_vcf"] + ".tbi",
        eval_vcf = rules.merge_snp_indel_vcf.output.vcf,
        eval_idx = rules.merge_snp_indel_vcf.output.idx
    output:
        os.path.join(out_path, "concordance_NA12878/snpsift/{group}.tsv")
    conda:
        "../envs/snpsift.yaml"
    log:
        os.path.join(out_path, "log/concordance_NA12878/snpsift/{group}.log")
    shell:
        """
        SnpSift concordance \
            -v {input.truth_vcf} \
            {input.eval_vcf} > {output} 2> {log}
        """
