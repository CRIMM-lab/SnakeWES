rule merge_tumor_nocontrol_mode:
	input:
		vcf=expand(f"results/{{sample}}_tumor/{{sample}}.{{caller}}.filt.vcf.gz", sample=config["samples"].values(), caller=config["callers"].values()),
		tbi=expand(f"results/{{sample}}_tumor/{{sample}}.{{caller}}.filt.vcf.gz.tbi", sample=config["samples"].values(), caller=config["callers"].values())
	output:
		vcf="results/multisample.{{caller}}.nocontrol.vcf.gz",
		tbi="results/multisample.{{caller}}.nocontrol.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/merge_varscan_nocontrol_mode.log"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2>{log} &&"
		"tabix -p vcf {output.vcf}"

