rule ParseAnnotationVepMutect2SingleSample:
	input:
		"results/{sample}.mutect2.paired.filtered.vep.vcf.gz"
	output:
		"results/{sample}.mutect2.paired.filtered.vep.tsv"
	threads: 1
	log:
		"logs/{sample}.ParseAnnotationVepMutect2SingleSample.log"
	conda:
		"../envs/bcftools.yaml"
	params:
		sample="{sample}_tumor",
		header="resources/header.vep.txt"
	shell:
		"""
		cat {params.header} <(bcftools view -s {params.sample} {input} |bcftools +split-vep -f "%CHROM\t%POS\t%REF\t%ALT\t%CSQ\t[%GT\t%AD\t%AF]\n" -d -A tab |tr "," "\t") > {output} 2>{log}
		"""

rule ParseAnnotationVepVarScanSingleSample:
	input:
		"results/{sample}.varscan.paired.vep.vcf.gz"
	output:
		"results/{sample}.varscan.paired.vep.tsv"
	threads: 1
	log:
		"logs/{sample}.ParseAnnotationVepVarScanSingleSample.log"
	conda:
		"../envs/bcftools.yaml"
	params:
		sample="{sample}_tumor",
		header="resources/header.vep.txt"
	shell:
		"""
		cat {params.header} <(bcftools view -s {params.sample} {input} |bcftools +split-vep -f "%CHROM\t%POS\t%REF\t%ALT\t%CSQ\t[%GT\t%RD\t%AD\t%FREQ]\n" -d -A tab) > {output} 2>{log}
		"""
