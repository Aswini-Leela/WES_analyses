rule gatk_combine_gvcfs:
	input:
		gvcfs = lambda wildcards: expand(rules.gatk_haplotypecaller.output, sample = get_samples_given_group(wildcards.group)),
		ref = rules.download_reference_genome.output,
		exon_bed = config["exon_bed"]
	output:
		vcf = temp(os.path.join(out_path, "gatk_combine/{group}.g.vcf.gz")),
		idx = temp(os.path.join(out_path, "gatk_combine/{group}.g.vcf.gz.tbi"))
	conda:
		"../envs/gatk.yaml"
	log:
		os.path.join(out_path, "log/gatk_combine_gvcfs/{group}.log")
	shell:
		"""
		gvcfs=$(for i in {input.gvcfs} ; do echo "-V "$i ; done | tr "\n" " ")

		gatk CombineGVCFs \
		        -R {input.ref} \
			$gvcfs \
		        -L {input.exon_bed} \
		        -O {output.vcf} 2> {log}
		"""

rule gatk_genotype_combined_gvcf:
	input:
		gvcf = rules.gatk_combine_gvcfs.output,
		ref = rules.download_reference_genome.output,
		exon_bed = config["exon_bed"]
	output:
		vcf = temp(os.path.join(out_path, "gatk_genotype/{group}.vcf.gz")),
		idx = temp(os.path.join(out_path, "gatk_genotype/{group}.vcf.gz.tbi"))
	conda:
		"../envs/gatk.yaml"
	log:
		os.path.join(out_path, "log/gatk_genotype/{group}.log")
	shell:
		"""
		gatk GenotypeGVCFs \
		        -R {input.ref} \
			-V {input.gvcf} \
		        -L {input.exon_bed} \
		        -O {output.vcf} 2> {log}
		"""
			
rule gatk_selectvariants_snp:
	input:
		vcf = rules.gatk_genotype_combined_gvcf.output,
		ref = rules.download_reference_genome.output,
	output:
		vcf = temp(os.path.join(out_path, "gatk_selectvariants_snp/{group}.vcf.gz")),
		idx = temp(os.path.join(out_path, "gatk_selectvariants_snp/{group}.vcf.gz.tbi"))
	conda:
		"../envs/gatk.yaml"
	log:
		os.path.join(out_path, "log/gatk_selectvariants_snp/{group}.log")
	params:
		select_type_to_include = "SNP"
	shell:
		"""
		gatk SelectVariants \
		        -R {input.ref} \
			-V {input.vcf} \
			--select-type-to-include {params.select_type_to_include} \
		        -O {output.vcf} 2> {log}
		"""

rule gatk_variant_recalibrator_snp:
	input:
		vcf = rules.gatk_selectvariants_snp.output,
		resource1 = config["hapmap_resource"],
		resource2 = config["omni_resource"],
		resource3 = config["1000G_resource"],
		resource4 = config["dbsnp_resource"]
	output:
		recal = temp(os.path.join(out_path, "gatk_variant_recalibrator_snp/{group}.recal")),
		tranche_file = temp(os.path.join(out_path, "gatk_variant_recalibrator_snp/{group}.tranches"))
	conda:
		"../envs/gatk.yaml"
	log:
		os.path.join(out_path, "log/gatk_variant_recalibrator_snp/{group}.log")
	resources:
		mem_mb = 50000
	params:
		java_options = "-Xmx50G",
		tranches = " 100.0 99.95 99.9 99.8 99.6 99.5 99.4 99.3 99.0 98.0 97.0 90.0",
		annots = " QD MQRankSum ReadPosRankSum FS MQ SOR DP",
		mode = "SNP",
		max_gaussians = 4,
		max_attempts = 5,
		resource1 = "hapmap,known=false,training=true,truth=true,prior=15",
		resource2 = "omni,known=false,training=true,truth=true,prior=12",
		resource3 = "1000G,known=false,training=true,truth=false,prior=10",
		resource4 = "dbsnp,known=true,training=false,truth=false,prior=7"
	shell:
		"""
		tranches=$(echo {params.tranches} | sed "s/ / -tranche /g" | sed "s/^/ -tranche /")
		annots=$(echo {params.annots} | sed "s/ / -an /g" | sed "s/^/ -an /")

		gatk VariantRecalibrator \
			--java-options \"{params.java_options}\" \
			-V {input.vcf} \
			--trust-all-polymorphic \
			$tranches \
			$annots \
			--mode {params.mode} \
			--max-attempts {params.max_attempts} \
			--max-gaussians {params.max_gaussians} \
			--resource:{params.resource1} {input.resource1} \
			--resource:{params.resource2} {input.resource2} \
			--resource:{params.resource3} {input.resource3} \
			--resource:{params.resource4} {input.resource4} \
			-O {output.recal} \
			--tranches-file {output.tranche_file} 2> {log}
		"""

rule gatk_selectvariants_indel:
	input:
		vcf = rules.gatk_genotype_combined_gvcf.output,
		ref = rules.download_reference_genome.output,
	output:
		vcf = temp(os.path.join(out_path, "gatk_selectvariants_indel/{group}.vcf.gz")),
		idx = temp(os.path.join(out_path, "gatk_selectvariants_indel/{group}.vcf.gz.tbi"))
	conda:
		"../envs/gatk.yaml"
	log:
		os.path.join(out_path, "log/gatk_selectvariants_indel/{group}.log")
	params:
		select_type_to_include = "INDEL"
	shell:
		"""
		gatk SelectVariants \
		        -R {input.ref} \
			-V {input.vcf} \
			--select-type-to-include {params.select_type_to_include} \
		        -O {output.vcf} 2> {log}
		"""

rule gatk_variant_recalibrator_indel:
	input:
		vcf = rules.gatk_selectvariants_indel.output,
		resource1 = config["mills_resource"],
		resource2 = config["axiompoly_resource"],
		resource3 = config["dbsnp_resource"]
	output:
		recal = temp(os.path.join(out_path, "gatk_variant_recalibrator_indel/{group}.recal")),
		tranche_file = temp(os.path.join(out_path, "gatk_variant_recalibrator_indel/{group}.tranches"))
	conda:
		"../envs/gatk.yaml"
	log:
		os.path.join(out_path, "log/gatk_variant_recalibrator_indel/{group}.log")
	params:
		java_options = "-Xmx50G",
		tranches = " 100.0 99.95 99.9 99.8 99.6 99.5 99.4 99.3 99.0 98.0 97.0 90.0",
		annots = " QD MQRankSum ReadPosRankSum FS MQ SOR DP",
		mode = "INDEL",
		max_attempts = 20,
		resource1 = "mills,known=false,training=true,truth=true,prior=12",
		resource2 = "axiomPoly,known=false,training=true,truth=false,prior=10",
		resource3 = "dbsnp,known=true,training=false,truth=false,prior=2"
	resources:
		mem_mb = 50000
	shell:
		"""
		tranches=$(echo {params.tranches} | sed "s/ / -tranche /g" | sed "s/^/ -tranche /")
		annots=$(echo {params.annots} | sed "s/ / -an /g" | sed "s/^/ -an /")

		gatk VariantRecalibrator \
			--java-options \"{params.java_options}\" \
			-V {input.vcf} \
			--trust-all-polymorphic \
			$tranches \
			$annots \
			--mode {params.mode} \
			--max-attempts {params.max_attempts} \
			--resource:{params.resource1} {input.resource1} \
			--resource:{params.resource2} {input.resource2} \
			--resource:{params.resource3} {input.resource3} \
			-O {output.recal} \
			--tranches-file {output.tranche_file} 2> {log}
		"""


rule gatk_applyVQSR_indel:
	input:
		vcf = rules.gatk_genotype_combined_gvcf.output,
		recal = rules.gatk_variant_recalibrator_indel.output.recal,
		tranch = rules.gatk_variant_recalibrator_indel.output.tranche_file,
	output:
		vcf = temp(os.path.join(out_path, "gatk_VQSR_indel/{group}.vcf.gz")),
		idx = temp(os.path.join(out_path, "gatk_VQSR_indel/{group}.vcf.gz.tbi"))
	conda:
		"../envs/gatk.yaml"
	log:
		os.path.join(out_path, "log/gatk_applyVQSR_indel/{group}.log")
	params:
		truth_sensitivity_filter_level = 99.7,
		create_output_variant_index = "true",
		mode = "INDEL"
	shell:
		"""
		gatk ApplyVQSR \
			-V {input.vcf} \
			--recal-file {input.recal} \
			--tranches-file {input.tranch} \
			--truth-sensitivity-filter-level {params.truth_sensitivity_filter_level} \
			--create-output-variant-index {params.create_output_variant_index} \
			-mode {params.mode}
			-O {output.vcf}
		"""

rule gatk_applyVQSR_snp:
	input:
		vcf = rules.gatk_genotype_combined_gvcf.output,
		recal = rules.gatk_variant_recalibrator_snp.output.recal,
		tranch = rules.gatk_variant_recalibrator_snp.output.tranche_file,
	output:
		vcf = temp(os.path.join(out_path, "gatk_VQSR_snp/{group}.vcf.gz")),
		idx = temp(os.path.join(out_path, "gatk_VQSR_snp/{group}.vcf.gz.tbi"))
	conda:
		"../envs/gatk.yaml"
	log:
		os.path.join(out_path, "log/gatk_applyVQSR_snp/{group}.log")
	params:
		truth_sensitivity_filter_level = 99.7,
		create_output_variant_index = "true",
		mode = "SNP"
	shell:
		"""
		gatk ApplyVQSR \
			-V {input.vcf} \
			--recal-file {input.recal} \
			--tranches-file {input.tranch} \
			--truth-sensitivity-filter-level {params.truth_sensitivity_filter_level} \
			--create-output-variant-index {params.create_output_variant_index} \
			-mode {params.mode} \
			-O {output.vcf}
		"""

