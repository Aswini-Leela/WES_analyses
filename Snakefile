import glob
import pandas as pd

configfile: "config.yaml"

metadata = pd.read_csv("input_folder/metadata.tsv", sep = "\t", index_col = 0)
samples = metadata.index.tolist()
groups = metadata["group"].unique().tolist()

include: "rules/downloads.smk"
include: "rules/basics.smk"
include: "rules/bwa.smk"
include: "rules/gatk_preprocessing.smk"
include: "rules/gatk_single_calling.smk"
# include: "rules/gatk_joint_calling.smk"
include: "rules/validate.smk"

rule all:
	input:
		expand(rules.gatk_haplotypecaller.output.vcf, sample = samples),
		expand(rules.funcotator.output.vcf, sample = samples),
		expand(rules.concordance_NA12878.output, sample = ["NIST7035", "NIST7086"])

def get_samples_given_group(group):
	return metadata.loc[metadata["group"] == group, :].index.tolist()

