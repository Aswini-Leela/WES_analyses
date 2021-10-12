# rule concatenate_fastq:
# 	input:
# 		dir = "input_folder/datasets/raw_data/{sample}"
# 	output:
# 		tmp_r1 = temp(os.path.join(config["intermediate_path"], "concatenate_fastq/_{sample}_1.fq.gz")),
# 		tmp_r2 = temp(os.path.join(config["intermediate_path"], "concatenate_fastq/_{sample}_2.fq.gz")),
# 		r1 = os.path.join(config["intermediate_path"], "concatenate_fastq/{sample}_1.fq.gz"),
# 		r2 = os.path.join(config["intermediate_path"], "concatenate_fastq/{sample}_2.fq.gz")
# 	conda:
# 		"../envs/bbmap.yaml"
# 	resources:
# 		ram="100G"
# 	log:
# 		"results/log/concatenate_fastq/{sample}.log"
# 	shell:
# 		"""
# 		for i in $(find {input.dir} | grep "_1.fq.gz") ; do
# 			cat $i >> {output.tmp_r1}
# 		done
# 		
# 		for i in $(find {input.dir} | grep "_2.fq.gz") ; do
# 			cat $i >> {output.tmp_r2}
# 		done
# 
# 		repair.sh in={output.tmp_r1} out={output.r1} repair -Xmx{resources.ram} 2> {log}
# 		repair.sh in={output.tmp_r2} out={output.r2} repair -Xmx{resources.ram} 2>> {log}
# 		"""


def get_samples_given_group(group):
	return metadata.loc[metadata["group"] == group, :].index.tolist()

