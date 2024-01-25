# Processing control samples

rule bcftools_call_on_controls:  
	input:
		bam="alignments/{sample}.control.dd.rec.bam",
		bai="alignments/{sample}.control.dd.rec.bai"
	output:
		"results/{sample}_control/{sample}.bcftools.vcf.gz"
	threads: 1
	params:
		ref=config["genome"],
		intervals=config["intervals"]
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.bcftools.call.log"
	shell:
		"bcftools mpileup -Ou -d 10000 -R {params.intervals} -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR -f {params.ref} {input.bam} | "
		"bcftools call -mv | "
		"bcftools norm -m - -f {params.ref} -Oz -o {output} 2> {log}"


rule bcftools_call_on_tumors:  
	input:
		bam="alignments/{sample}.tumor.dd.rec.bam",
		bai="alignments/{sample}.tumor.dd.rec.bai"
	output:
		"results/{sample}_tumor/{sample}.bcftools.vcf.gz"
	threads: 1
	params:
		ref=config["genome"],
		intervals=config["intervals"]
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.bcftools.call.log"
	shell:
		"bcftools mpileup -Ou -d 10000 -R {params.intervals} -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR -f {params.ref} {input.bam} | "
		"bcftools call -mv | "
		"bcftools norm -m - -f {params.ref} -Oz -o {output} 2> {log}"

rule bcftools_add_af_field_on_controls:  
	input:
		"results/{sample}_control/{sample}.bcftools.vcf.gz"
	output:
		vcf=temp("results/{sample}_control/{sample}.bcftools.addaf.vcf.gz"),
		tsv=temp("results/{sample}_control/{sample}.annot.tsv.gz"),
		tsv_tbi=temp("results/{sample}_control/{sample}.annot.tsv.gz.tbi")
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.bcftools.addaf.log"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%AD{{0}}\t%AD{{1}}]' {input} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {output.tsv} && 
		tabix -b2 -e2 {output.tsv} && 
		bcftools annotate -a {output.tsv} -h <(echo '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,AF {input} -Oz -o {output.vcf}
		"""

rule bcftools_add_af_field_on_tumors:  
	input:
		"results/{sample}_tumor/{sample}.bcftools.vcf.gz"
	output:
		vcf=temp("results/{sample}_tumor/{sample}.bcftools.addaf.vcf.gz"),
		tsv=temp("results/{sample}_tumor/{sample}.annot.tsv.gz"),
		tsv_tbi=temp("results/{sample}_tumor/{sample}.annot.tsv.gz.tbi")
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.tum.bcftools.addaf.log"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%AD{{0}}\t%AD{{1}}]' {input} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {output.tsv} && 
		tabix -b2 -e2 {output.tsv} && 
		bcftools annotate -a {output.tsv} -h <(echo '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,AF {input} -Oz -o {output.vcf}
		"""

rule bcftools_filter_controls: 
	input:
		"results/{sample}_control/{sample}.bcftools.addaf.vcf.gz"
	output:
		vcf=temp("results/{sample}_control/{sample}.bcftools.filt.vcf.gz"),
		tbi=temp("results/{sample}_control/{sample}.bcftools.filt.vcf.gz.tbi")
	threads: 1
	params:
		excl=config["chr_to_exclude"],
		depth=config['filtering_controls']['min_depth'],
		vaf=config['filtering_controls']['vaf'],
		alt=config['filtering_controls']['alt_depth']		
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.bcftools.filter.log"
	shell:
		"bcftools view -i'FORMAT/DP >= {params.depth} & INFO/AF >= {params.vaf} & FORMAT/AD[0:1] >= {params.alt}' {input} | "
		"grep -v -f {params.excl} | "
		"bcftools sort -Oz -o {output.vcf} 2> {log} && "
		"tabix -p vcf {output.vcf}"

rule bcftools_filter_tumors: 
	input:
		"results/{sample}_tumor/{sample}.bcftools.addaf.vcf.gz"
	output:
		vcf=temp("results/{sample}_tumor/{sample}.bcftools.filt.vcf.gz"),
		tbi=temp("results/{sample}_tumor/{sample}.bcftools.filt.vcf.gz.tbi")
	threads: 1
	params:
		config["chr_to_exclude"],
		excl=config["chr_to_exclude"],
		depth=config['filtering_tumors']['min_depth'],
		vaf=config['filtering_tumors']['vaf'],
		alt=config['filtering_tumors']['alt_depth']
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.bcftools.filter.log"
	shell:
		"bcftools view -i 'FORMAT/DP >= {params.depth} & INFO/AF >= {params.vaf} & FORMAT/AD[0:1] >= {params.alt}' {input} | "
		"grep -v -f {params.excl} | "
		"bcftools sort -Oz -o {output.vcf} 2> {log} && "
		"tabix -p vcf {output.vcf}"

rule bcftools_merge_bcftools_controls:
	input:
		vcf=expand(f"results/{{sample}}_control/{{sample}}.bcftools.filt.vcf.gz", sample=config['controls'].values()),
		tbi=expand(f"results/{{sample}}_control/{{sample}}.bcftools.filt.vcf.gz.tbi", sample=config['controls'].values())
	output:
		vcf="results/custom_pon.bcftools.vcf.gz",
		tbi="results/custom_pon.bcftools.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/bcftools.merge.normals.vars.log"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2> {log} && "
		"tabix -p vcf {output.vcf}"


rule bcftools_remove_controls: 
	input:
		vcf="results/{sample}_tumor/{sample}.bcftools.filt.vcf.gz",
		vcf_tbi="results/{sample}_tumor/{sample}.bcftools.filt.vcf.gz.tbi",
		custom_pon="results/custom_pon.bcftools.vcf.gz",
		pon_tbi="results/custom_pon.bcftools.vcf.gz.tbi"
	output:
		vcf="results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/{sample}.filtOnCtr.vars.log"
	shell:
		"bcftools isec -w1 -Oz -c none -n~10 -o {output.vcf} {input.vcf} {input.custom_pon} 2> {log} && "
		"tabix -p vcf {output.vcf}"	


rule merge_bcftools_samples_variants: 
	input:
		vcf=expand("results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vcf.gz", sample=config['tumors'].values()),
		tbi=expand("results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vcf.gz.tbi", sample=config['tumors'].values())
	output:
		vcf="results/multisample.bcftools.vcf.gz",
		tbi="results/multisample.bcftools.vcf.gz.tbi"
	threads: 1
	log: 
		"logs/multisample.bcftools.log"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2> {log} &&"
		"tabix -p vcf {output.vcf}"
