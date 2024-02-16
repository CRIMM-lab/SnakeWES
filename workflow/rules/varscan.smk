# true command line parameters "--min-avg-qual 15 --p-value 0.1 --min-var-freq 0.01" instead of "--min-avg-qual 1 --p-value 1 --min-var-freq 0 --min-coverage 1 --min-reads2 1"

rule varscan_call_snv_on_controls: 
	input:
		bam="alignments/{sample}.control.dd.rec.bam",
		bai="alignments/{sample}.control.dd.rec.bai"
	output:
		vcf_snv="results/{sample}_control/{sample}.varscan.snv.vcf.gz",
		tbi_snv="results/{sample}_control/{sample}.varscan.snv.vcf.gz.tbi"
	message:
		"varscan call snv - control {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.varscan.snv.callandNorm.txt"
	threads: 1
	params:
		intervals=config["intervals"],
		ref=config["genome"],
	conda:
		"../envs/varscan.yaml"
	log:
		"logs/{sample}.varscan.snv.callandNorm.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} | "
		"varscan mpileup2snp --output-vcf 1 --min-avg-qual 1 --p-value 1 --min-var-freq 0 --min-coverage 1 --min-reads2 1 | "
		"bgzip -c > {output.vcf_snv} 2> {log} &&"
		"tabix -p vcf {output.vcf_snv}"

rule varscan_call_snv_on_tumors: 
	input:
		bam="alignments/{sample}.tumor.dd.rec.bam",
		bai="alignments/{sample}.tumor.dd.rec.bai"
	output:
		vcf_snv="results/{sample}_tumor/{sample}.varscan.snv.vcf.gz",
		tbi_snv="results/{sample}_tumor/{sample}.varscan.snv.vcf.gz.tbi"
	threads: 1
	params:
		#varscan="/home/simone/programs/VarScan.v2.3.9.jar",
		intervals=config["intervals"],
		ref=config["genome"]
	conda:
		"../envs/varscan.yaml"
	log:
		"logs/{sample}.varscan.snv.callandNorm.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} |"
		"varscan mpileup2snp --output-vcf --min-avg-qual 1 --p-value 1 --min-var-freq 0 --min-coverage 1 --min-reads2 1 |"
		"bgzip -c > {output.vcf_snv} &&"
		"tabix -p vcf {output.vcf_snv}"


rule varscan_call_indel_on_controls: 
	input:
		bam="alignments/{sample}.control.dd.rec.bam",
		bai="alignments/{sample}.control.dd.rec.bai"
	output:
		vcf_indel="results/{sample}_control/{sample}.varscan.indel.vcf.gz",
		tbi_indel="results/{sample}_control/{sample}.varscan.indel.vcf.gz.tbi"
	message:
		"varscan call indel - control {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.varscan.indel.callandNorm.txt"
	threads: 1
	params:
		intervals=config["intervals"],
		ref=config["genome"]
	conda:
		"../envs/varscan.yaml"
	log:
		"logs/{sample}.varscan.indel.callandNorm.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} |"
		"varscan mpileup2indel --output-vcf --min-avg-qual 1 --p-value 1 --min-var-freq 0 --min-coverage 1 --min-reads2 1 |"
		"bgzip -c > {output.vcf_indel} 2> {log} &&"
		"tabix -p vcf {output.vcf_indel}"
		
rule varscan_call_indel_on_tumors: 
	input:
		bam="alignments/{sample}.tumor.dd.rec.bam",
		bai="alignments/{sample}.tumor.dd.rec.bai"
	output:
		vcf_indel="results/{sample}_tumor/{sample}.varscan.indel.vcf.gz",
		tbi_indel="results/{sample}_tumor/{sample}.varscan.indel.vcf.gz.tbi"
	message:
		"varscan call indel - tumor {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.varscan.indel.callandNorm.txt"
	threads: 1
	params:
		#varscan="/home/simone/programs/VarScan.v2.3.9.jar",
		intervals=config["intervals"],
		ref=config["genome"]
	conda:
		"../envs/varscan.yaml"
	log:
		"logs/{sample}.varscan.indel.callandNorm.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} | "
		"varscan mpileup2indel --output-vcf --min-avg-qual 1 --p-value 1 --min-var-freq 0 --min-coverage 1 --min-reads2 1 | "
		"bgzip -c > {output.vcf_indel} && "
		"tabix -p vcf {output.vcf_indel} "

rule varscan_concat_controls: 
	input:
		vcf_snv="results/{sample}_control/{sample}.varscan.snv.vcf.gz",
		vcf_indel="results/{sample}_control/{sample}.varscan.indel.vcf.gz",
		tbi_snv="results/{sample}_control/{sample}.varscan.snv.vcf.gz.tbi",
		tbi_indel="results/{sample}_control/{sample}.varscan.indel.vcf.gz.tbi"
	output:
		temp("results/{sample}_control/{sample}.varscan.concat.tmp.vcf.gz")
	message:
		"concat snv and indel - control {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.varscan.concat.txt"
	params: 
		config["genome"]
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.varscan.concat.log"
	shell: 
		"bcftools concat -a -Ov {input.vcf_snv} {input.vcf_indel} | "
		"bcftools norm -Oz -m - -f {params} -Oz -o {output} 2> {log}"				


rule varscan_concat_tumors: 
	input:
		vcf_snv="results/{sample}_tumor/{sample}.varscan.snv.vcf.gz",
		vcf_indel="results/{sample}_tumor/{sample}.varscan.indel.vcf.gz",
		tbi_snv="results/{sample}_tumor/{sample}.varscan.snv.vcf.gz.tbi",
		tbi_indel="results/{sample}_tumor/{sample}.varscan.indel.vcf.gz.tbi"
	output:
		temp("results/{sample}_tumor/{sample}.varscan.concat.tmp.vcf.gz")
	message:
		"concat snv and indel - tumor {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.varscan.concat.txt"
	params: 
		config["genome"]
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.varscan.concat.log"
	shell: 
		"bcftools concat -a -Ov {input.vcf_snv} {input.vcf_indel} | "
		"bcftools norm -Oz -m - -f {params} -Oz -o {output} 2> {log}"

rule varscan_mod_AF_controls:
	input:
		"results/{sample}_control/{sample}.varscan.concat.tmp.vcf.gz"
	output:
		vcf="results/{sample}_control/{sample}.varscan.concat.vcf.gz",
		tsv="results/{sample}_control/{sample}.varscan.annot.tsv.gz",
		tsv_tbi=temp("results/{sample}_control/{sample}.varscan.annot.tsv.gz.tbi")
	log:
		"logs/{sample}.varscan_mod_AF_controls.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%FREQ]' {input} | 
		awk -F'\t' 'BEGIN {{OFS="\t"}} {{ $5 = $5 / 100; print }}'| bgzip -c > {output.tsv} && 
		tabix -b2 -e2 {output.tsv} &&
		bcftools annotate -x FORMAT/FREQ -a {output.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF  {input} -Oz -o {output.vcf}
		"""

rule varscan_mod_AF_tumors:
	input:
		"results/{sample}_tumor/{sample}.varscan.concat.tmp.vcf.gz"
	output:
		vcf="results/{sample}_tumor/{sample}.varscan.concat.vcf.gz",
		tsv="results/{sample}_tumor/{sample}.varscan.annot.tsv.gz",
		tsv_tbi="results/{sample}_tumor/{sample}.varscan.annot.tsv.gz.tbi"
	log:
		"logs/{sample}.varscan_mod_AF_tumors.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%FREQ]' {input} | awk -F'\t' 'BEGIN {{OFS="\t"}} {{ $5 = $5 / 100; print }}'| bgzip -c > {output.tsv} && 
		tabix -b2 -e2 {output.tsv} &&
		bcftools annotate -x FORMAT/FREQ -a {output.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF  {input} -Oz -o {output.vcf}
		"""

rule varscan_reheader_controls:
	input:
		"results/{sample}_control/{sample}.varscan.concat.vcf.gz"
	output:
		"results/{sample}_control/{sample}.varscan.vcf.gz"
	message:
		"reheader - control {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.varscan.rh.txt"
	threads: 1
	params:
		txt="results/{sample}_control/{sample}.rh.txt"
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.varscan.rh.log"
	shell: 
		"echo {wildcards.sample} > {params.txt} && "
		"bcftools reheader -s {params.txt} -o {output} {input} && "
		"rm {params.txt}"


rule varscan_reheader_tumors:
	input:
		"results/{sample}_tumor/{sample}.varscan.concat.vcf.gz"
	output:
		"results/{sample}_tumor/{sample}.varscan.vcf.gz"
	message:
		"reheader - tumor {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.varscan.rh.txt"
	threads: 1
	params:
		txt="results/{sample}_tumor/{sample}.rh.txt"
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.varscan.rh.log"
	shell: 
		"echo {wildcards.sample} > {params.txt} && "
		"bcftools reheader -s {params.txt} -o {output} {input} && "
		"rm {params.txt}"


rule filter_varscan_controls_vars: 
	input:
		"results/{sample}_control/{sample}.varscan.vcf.gz"
	output:
		vcf="results/{sample}_control/{sample}.varscan.filt.vcf.gz",
		tbi="results/{sample}_control/{sample}.varscan.filt.vcf.gz.tbi"
	message:
		"filter varscan controls vars - control {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.filter.varscan.txt"
	threads: 1
	params:
		excl=config["chr_to_exclude"],
		alt=config["filtering_controls"]["alt_depth"],
		min_depth=config["filtering_controls"]["min_depth"],
		vaf=config["filtering_controls"]["vaf"],
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.filter.varscan.log"
	shell:
		"bcftools view -Ov -i 'FORMAT/DP >= {params.min_depth} & FORMAT/AD >= {params.alt} & FORMAT/AF >= {params.vaf}' {input} |"
		"grep -v -f {params.excl} |"
		"bcftools sort -Oz -o {output.vcf} 2> {log} &&"
		"tabix -f -p vcf {output.vcf}"

rule filter_varscan_tumors_vars: 
	input:
		"results/{sample}_tumor/{sample}.varscan.vcf.gz"
	output:
		vcf="results/{sample}_tumor/{sample}.varscan.filt.vcf.gz"
		tbi="results/{sample}_tumor/{sample}.varscan.filt.vcf.gz.tbi"
	message:
		"filter varscan tumor vars - tumor {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.filter.varscan.txt"
	threads: 1
	params:
		excl=config["chr_to_exclude"],
		alt=config["filtering_tumors"]["alt_depth"],
		min_depth=config["filtering_tumors"]["min_depth"],
		vaf=config["filtering_tumors"]["vaf"],
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.filter.varscan.log"
	shell:
		"bcftools view -Ov -i'FORMAT/DP >= {params.min_depth} & FORMAT/AD >= {params.alt} & FORMAT/AF >= {params.vaf}' {input} | "
		"grep -v -f {params.excl} | "
		"bcftools sort -Oz -o {output.vcf} && "
		"tabix -f -p vcf {output.vcf}"

rule merge_varscan_controls_variants:
	input:
		vcf=expand(f"results/{{sample}}_control/{{sample}}.varscan.filt.vcf.gz", sample=config['controls'].values()),
		tbi=expand(f"results/{{sample}}_control/{{sample}}.varscan.filt.vcf.gz.tbi", sample=config['controls'].values())
	output:
		vcf="results/custom_pon.varscan.vcf.gz",
		tbi="results/custom_pon.varscan.vcf.gz.tbi"
	message:
		"merge varscan controls vars"
	benchmark:
		"benchmarks/varscan.merge.ctr.vars.txt"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/varscan.merge.ctr.vars.log"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2> {log} &&"
		"tabix -p vcf {output.vcf} "


rule remove_varscan_controls_vars: 
	input:
		vcf="results/{sample}_tumor/{sample}.varscan.filt.vcf.gz",
		custom_pon="results/custom_pon.varscan.vcf.gz",
		tbi="results/custom_pon.varscan.vcf.gz.tbi"
	output:
		vcf="results/{sample}_tumor/{sample}.varscan.filtOnCtr.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.varscan.filtOnCtr.vcf.gz.tbi"
	message:
		"remove varscan controls vars"
	benchmark:
		"benchmarks/{sample}.varscan.rh.txt"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/{sample}.filtOnCtr.vars.log"
	shell:
		"bcftools isec -w1 -Oz -c none -n~10 -o {output.vcf} {input.vcf} {input.custom_pon} 2> {log} && "
		"tabix -p vcf {output.vcf}"


rule merge_varscan_tumors_variants: ## difference between single curly and doucle curly brackets ##
	input:
		vcf=expand(f"results/{{sample}}_tumor/{{sample}}.varscan.filtOnCtr.vcf.gz", sample=config['samples'].values()),
		tbi=expand(f"results/{{sample}}_tumor/{{sample}}.varscan.filtOnCtr.vcf.gz.tbi", sample=config['samples'].values())
	output:
		vcf="results/multisample.varscan.vcf.gz",
		tbi="results/multisample.varscan.vcf.gz.tbi"
	message:
		"merge varscan tumors vars"
	benchmark:
		"benchmarks/multisample.varscan.txt"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/multisample.varscan.log"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2> {log} &&"
		"tabix -p vcf {output.vcf}"


