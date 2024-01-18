rule varscan_call_snv_on_controls: 
	input:
		bam="data/bam/{control}.dedup.recal.bam",
		bai="data/bam/{control}.dedup.recal.bam.bai"
	output:
		vcf_snv=temp("results/varscan/custom_pon/{control}.ctr.varscan.snv.vcf.gz"),
		tbi_snv=temp("results/varscan/custom_pon/{control}.ctr.varscan.snv.vcf.gz.tbi")
	threads: 1
	params:
		varscan="/home/simone/programs/VarScan.v2.3.9.jar",
		intervals=config["intervals"],
		ref=config["genome"],
		#vcf_tmp=temp("results/varscan/custom_pon/{control}.ctr.varscan.snv.vcf")
	log:
		"logs/{control}.varscan.snv.callandNorm.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} | "
		"java -jar {params.varscan} mpileup2snp --output-vcf 1 --min-avg-qual 15 --p-value 0.1 --min-var-freq 0.01 | "
		"bgzip -c > {output.vcf_snv} 2> {log} &&"
		"tabix -p vcf {output.vcf_snv}"

rule varscan_call_indel_on_controls: 
	input:
		bam="data/bam/{control}.dedup.recal.bam",
		bai="data/bam/{control}.dedup.recal.bam.bai"
	output:
		vcf_indel=temp("results/varscan/custom_pon/{control}.ctr.varscan.indel.vcf.gz"),
		tbi_indel=temp("results/varscan/custom_pon/{control}.ctr.varscan.indel.vcf.gz.tbi")
	threads: 1
	params:
		varscan="/home/simone/programs/VarScan.v2.3.9.jar",
		intervals=config["intervals"],
		ref=config["genome"]
	log:
		"logs/{control}.varscan.indel.callandNorm.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} |"
		"java -jar {params.varscan} mpileup2indel --output-vcf --min-avg-qual 15 --p-value 0.1 --min-var-freq 0.01 |"
		"bgzip -c > {output.vcf_indel} 2> {log} &&"
		"tabix -p vcf {output.vcf_indel}"
		
rule varscan_concat_controls: 
	input:
		vcf_snv="results/varscan/custom_pon/{control}.ctr.varscan.snv.vcf.gz",
		vcf_indel="results/varscan/custom_pon/{control}.ctr.varscan.indel.vcf.gz",
		tbi_snv="results/varscan/custom_pon/{control}.ctr.varscan.snv.vcf.gz.tbi",
		tbi_indel="results/varscan/custom_pon/{control}.ctr.varscan.indel.vcf.gz.tbi"
	output:
		vcf_concat=temp("results/varscan/custom_pon/{control}.ctr.varscan.concat.vcf.gz")
	params: 
		vcf_tmp="results/varscan/custom_pon/{control}.ctr.varscan.tmp.vcf",
		ref=config["genome"]
	threads: 1
	log:
		"logs/{control}.varscan.concat.log"
	shell: 
		"bcftools concat -a -Ov -o {params.vcf_tmp} {input.vcf_snv} {input.vcf_indel} && "
		"cat <(bcftools view -h {params.vcf_tmp} | sed -e 's/<ID=FREQ,Number=1,Type=String/<ID=FREQ,Number=A,Type=Float/g') <(bcftools view -H {params.vcf_tmp} | sed 's/%//g' | sed 's/,/./g') | "
		"bcftools norm -Oz -m - -f {params.ref} -Oz -o {output.vcf_concat} 2> {log} && "
		"rm {params.vcf_tmp}"

rule varscan_reheader_controls:
	input:
		"results/varscan/custom_pon/{control}.ctr.varscan.concat.vcf.gz"
	output:
		"results/varscan/custom_pon/{control}.ctr.varscan.vcf.gz"
	threads: 1
	params:
		txt="results/varscan/custom_pon/{control}.rh.txt"
	log:
		"logs/{control}.varscan.rh.log"
	shell: 
		"echo {wildcards.control} > {params.txt} && "
		"bcftools reheader -s {params.txt} -o {output} {input} && "
		"rm {params.txt}"

rule filter_varscan_controls_vars: 
	input:
		"results/varscan/custom_pon/{control}.ctr.varscan.vcf.gz"
	output:
		vcf="results/varscan/custom_pon/{control}.ctr.varscan.filt.vcf.gz",
		tbi="results/varscan/custom_pon/{control}.ctr.varscan.filt.vcf.gz.tbi"
	threads: 1
	params:
		chr_to_exclude=config["chr_to_exclude"]
	log:
		"logs/{control}.filter.varscan.log"
	shell:
		"bcftools view -Ov -i 'FORMAT/ADF >= 1 && FORMAT/ADR >= 1 && FREQ >= 1' {input} |"
		"grep -v -f {params.chr_to_exclude} |"
		"bcftools sort -Oz -o {output.vcf} 2> {log} &&"
		"tabix -f -p vcf {output.vcf}"

rule merge_varscan_controls_variants:
	input:
		vcf=expand(f"results/varscan/custom_pon/{{control}}.ctr.varscan.filt.vcf.gz", control=config['controls'].values()),
		tbi=expand(f"results/varscan/custom_pon/{{control}}.ctr.varscan.filt.vcf.gz.tbi", control=config['controls'].values())
	output:
		vcf="results/varscan/custom_pon/custom_pon.varscan.vcf.gz",
		tbi="results/varscan/custom_pon/custom_pon.varscan.vcf.gz.tbi"
	threads: 1
	log: 
		"logs/varscan.merge.ctr.vars.log"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2> {log} &&"
		"tabix -p vcf {output.vcf} "

#---------------------------------------------------------------------------------------------------
rule varscan_call_snv_on_samples: 
	input:
		bam="data/bam/{sample}.dedup.recal.bam",
		bai="data/bam/{sample}.dedup.recal.bam.bai"
	output:
		vcf_snv=temp("results/varscan/{sample}.sample.varscan.snv.vcf.gz"),
		tbi_snv=temp("results/varscan/{sample}.sample.varscan.snv.vcf.gz.tbi")
	threads: 1
	params:
		varscan="/home/simone/programs/VarScan.v2.3.9.jar",
		intervals=config["intervals"],
		ref=config["genome"]
	log:
		"logs/{sample}.varscan.snv.callandNorm.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} |"
		"java -jar {params.varscan} mpileup2snp --output-vcf --min-avg-qual 15 --p-value 0.1 --min-var-freq 0.01 |"
		"bgzip -c > {output.vcf_snv} &&"
		"tabix -p vcf {output.vcf_snv}"

rule varscan_call_indel_on_samples: 
	input:
		bam="data/bam/{sample}.dedup.recal.bam",
		bai="data/bam/{sample}.dedup.recal.bam.bai"
	output:
		vcf_indel=temp("results/varscan/{sample}.sample.varscan.indel.vcf.gz"),
		tbi_indel=temp("results/varscan/{sample}.sample.varscan.indel.vcf.gz.tbi")
	threads: 1
	params:
		varscan="/home/simone/programs/VarScan.v2.3.9.jar",
		intervals=config["intervals"],
		ref=config["genome"]
	log:
		"logs/{sample}.varscan.indel.callandNorm.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} | "
		"java -jar {params.varscan} mpileup2indel --output-vcf --min-avg-qual 15 --p-value 0.1 --min-var-freq 0.01 | "
		"bgzip -c > {output.vcf_indel} && "
		"tabix -p vcf {output.vcf_indel} "

rule varscan_concat_samples: 
	input:
		vcf_snv="results/varscan/{sample}.sample.varscan.snv.vcf.gz",
		vcf_indel="results/varscan/{sample}.sample.varscan.indel.vcf.gz",
		tbi_snv="results/varscan/{sample}.sample.varscan.snv.vcf.gz.tbi",
		tbi_indel="results/varscan/{sample}.sample.varscan.indel.vcf.gz.tbi"
	output:
		vcf_concat="results/varscan/{sample}.sample.varscan.concat.vcf.gz"
	params: 
		vcf_tmp="results/varscan/{sample}.sample.varscan.tmp.vcf", 
		ref=config["genome"]
	threads: 1
	log:
		"logs/{sample}.varscan.concat.log"
	shell: 
		"bcftools concat -a -Ov -o {params.vcf_tmp} {input.vcf_snv} {input.vcf_indel} && "
		"cat <(bcftools view -h {params.vcf_tmp} | sed -e 's/<ID=FREQ,Number=1,Type=String/<ID=FREQ,Number=A,Type=Float/g') <(bcftools view -H {params.vcf_tmp} | sed 's/%//g' | sed 's/,/./g') | "
		"bcftools norm -Oz -m - -f {params.ref} -Oz -o {output.vcf_concat} 2> {log} && "
		"rm {params.vcf_tmp} "

rule varscan_reheader_samples:
	input:
		"results/varscan/{sample}.sample.varscan.concat.vcf.gz"
	output:
		"results/varscan/{sample}.sample.varscan.vcf.gz"
	threads: 1
	params:
		txt="results/varscan/{sample}.rh.txt"
	log:
		"logs/{sample}.varscan.rh.log"
	shell: 
		"echo {wildcards.sample} > {params.txt} && "
		"bcftools reheader -s {params.txt} -o {output} {input} && "
		"rm {params.txt}"

rule filter_varscan_samples_vars: 
	input:
		"results/varscan/{sample}.sample.varscan.vcf.gz"
	output:
		"results/varscan/{sample}.sample.varscan.filt.vcf.gz"
	threads: 1
	params:
		varscan="",
		ref=config["genome"],
		chr_to_exclude=config["chr_to_exclude"]
	log:
		"logs/{sample}.filter.varscan.log"
	shell:
		"bcftools view -Ov -i'FORMAT/DP >= 20 && FORMAT/ADF >= 1 && FORMAT/ADR >= 1 && FREQ >= 3' {input} | "
		"grep -v -f {params.chr_to_exclude} | "
		"bcftools sort -Oz -o {output} && "
		"tabix -f -p vcf {output}"

rule remove_varscan_controls_vars: 
	input:
		vcf="results/varscan/{sample}.sample.varscan.filt.vcf.gz",
		custom_pon="results/varscan/custom_pon/custom_pon.varscan.vcf.gz",
		tbi="results/varscan/custom_pon/custom_pon.varscan.vcf.gz.tbi"
	output:
		vcf="results/varscan/{sample}.sample.varscan.filtOnCtr.vcf.gz",
		tbi="results/varscan/{sample}.sample.varscan.filtOnCtr.vcf.gz.tbi"
	threads: 1
	params:
		gatk_pon=config["gatk_pon"]
	log: 
		"logs/{sample}.filtOnCtr.vars.log"
	shell:
		"bcftools isec -w1 -Oz -c none -n~100 -o {output.vcf} {input.vcf} {input.custom_pon} {params.gatk_pon} 2> {log} && "
		"tabix -p vcf {output.vcf}"

rule merge_varscan_samples_variants: ## difference between single curly and doucle curly brackets ##
    input:
        vcf=expand("results/varscan/{sample}.sample.varscan.filtOnCtr.vcf.gz", sample=config['samples'].values()),
        tbi=expand("results/varscan/{sample}.sample.varscan.filtOnCtr.vcf.gz.tbi", sample=config['samples'].values())
    output:
        vcf="results/varscan/multisample.varscan.vcf.gz",
        tbi="results/varscan/multisample.varscan.vcf.gz.tbi"
    threads: 1
    log: 
        "logs/multisample.varscan.log"
    shell:
        "bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2> {log} &&"
        "tabix -p vcf {output.vcf}"
