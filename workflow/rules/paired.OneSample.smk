rule ParseAnnotationVepMutect2SingleSample:
	input:
		"results/{sample}.mutect2.paired.filtered.vep.vcf.gz"
	output:
		"results/{sample}.mutect2.paired.filtered.vep.tmp.tsv"
	threads: 1
	log:
		"logs/{sample}.ParseAnnotationVepMutect2SingleSample.log"
	conda:
		"../envs/bcftools.yaml"
	params:
		""
	shell:
		"""
		bcftools +split-vep {input} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%GT\t%RD\t%AD\t%FREQ]\n" -d -A tab > {output} 2>{log}
		"""

rule ParseAnnotationVepVarScanSingleSample:
	input:
		"results/{sample}.varscan.paired.vep.vcf.gz"
	output:
		"results/{sample}.varscan.paired.vep.tmp.tsv"
	threads: 1
	log:
		"logs/{sample}.ParseAnnotationVepVarScanSingleSample.log"
	conda:
		"../envs/bcftools.yaml"
	params:
		""
	shell:
		"""
		bcftools +split-vep {input} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%GT\t%RD\t%AD\t%FREQ]\n" -d -A tab > {output} 2>{log}
		"""
