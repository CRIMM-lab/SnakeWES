rule freebayes_call_and_norm_on_samples: 
	input:
		bam="alignments/{sample}.tumor.dd.rec.bam",
		bai="alignments/{sample}.tumor.dd.rec.bai"
	output:
		"results/{sample}_tumor/{sample}.freebayes.vcf.gz"
	threads: 1
	conda:
		"../envs/freebayes.yaml"
	params:
		ref=config["genome"],
		intervals=config["intervals"]
	log:
		"logs/{sample}.freebayes_call_and_norm_on_samples.log"
	shell:
		"freebayes -f {params.ref} -F 0.01 -C 2 -t {params.intervals} --pooled-continuous {input.bam} | "
		"bcftools norm -m - -f {params.ref} -Oz -o {output} 2>{log}"


rule freebayes_call_and_norm_on_controls: 
	input:
		bam="alignments/{sample}.control.dd.rec.bam",
		bai="alignments/{sample}.control.dd.rec.bai"
	output:
		"results/{sample}_control/{sample}.freebayes.vcf.gz"
	threads: 1
	conda:
		"../envs/freebayes.yaml"
	params:
		ref=config["genome"],
		intervals=config["intervals"]
	log:
		"logs/{sample}.freebayes_call_and_norm_on_controls.log"
	shell:
		"freebayes -f {params.ref} -F 0.01 -C 2 -t {params.intervals} --pooled-continuous {input.bam} | "
		"bcftools norm -m - -f {params.ref} -Oz -o {output} 2>{log}"


rule freebayes_add_af_samples:
	input:
		"results/{sample}_tumor/{sample}.freebayes.vcf.gz"
	output:
		"results/{sample}_tumor/{sample}.freebayes.annotAF.vcf.gz"
	message:
		"freebayes add AF on tumor {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.freebayesAddAFtumor.txt"
	log:
		"logs/{sample}.freebayesAddAFtumor.log",
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	params:
		"results/{sample}_tumor/{sample}.freebayes.annotAF.tsv.gz"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%RO\t%AO\n' {input} 2> {log} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {params} 2>> {log} && 
		tabix -b2 -e2 {params} && bcftools annotate -a {params} --columns CHROM,POS,REF,ALT,AF {input} -Oz -o {output} 2>> {log} && 
		rm {params}
		"""

rule freebayes_add_af_control:
	input:
		"results/{sample}_control/{sample}.freebayes.vcf.gz"
	output:
		"results/{sample}_control/{sample}.freebayes.annotAF.vcf.gz"
	message:
		"freebayes add AF on control {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.freebayesAddAFcontrol.txt"
	log:
		"logs/{sample}.freebayesAddAFcontrol.log",
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	params:
		"results/{sample}_control/{sample}.freebayes.annotAF.tsv.gz"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%RO\t%AO\n' {input} 2>{log} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {params} 2>>{log} && 
		tabix -b2 -e2 {params} && bcftools annotate -a {params} --columns CHROM,POS,REF,ALT,AF {input} -Oz -o {output} 2>>{log} && 
		rm {params}
		"""

rule freebayes_filter_samples: 
	input:
		"results/{sample}_tumor/{sample}.freebayes.annotAF.vcf.gz"
	output:
		"results/{sample}_tumor/{sample}.freebayes.filt.vcf.gz"
	message:
		"freebayes filter on tumor {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.freebayesFilterTumor.txt"
	log:
		"logs/{sample}.freebayesFilterTumor.log",
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	params:
		excl=config["chr_to_exclude"],
		depth=config['filtering_tumors']['min_depth'],
		vaf=config['filtering_tumors']['vaf'],
		alt=config['filtering_tumors']['alt_depth']
	shell:
		"bcftools view -i 'QUAL > 1 & INFO/DP >= {params.depth} & INFO/AF >= {params.vaf} & INFO/AC >= {params.alt}' {input} 2> {log}| "
		"grep -v -f {params.excl} 2>> {log} | "
		"bcftools sort -Oz -o {output} 2>> {log} && "
		"tabix -p vcf {output} 2>> {log}"

rule freebayes_filter_control: 
	input:
		"results/{sample}_control/{sample}.freebayes.annotAF.vcf.gz"
	output:
		vcf="results/{sample}_control/{sample}.freebayes.filt.vcf.gz",
		tbi="results/{sample}_control/{sample}.freebayes.filt.vcf.gz.tbi"
	message:
		"freebayes filter on control {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.freebayesFilterControl.txt"
	log:
		"logs/{sample}.freebayesFilterControl.log",
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	params:
		excl=config["chr_to_exclude"],
		depth=config['filtering_controls']['min_depth'],
		vaf=config['filtering_controls']['vaf'],
		alt=config['filtering_controls']['alt_depth']
	shell:
		"bcftools view -i 'QUAL > 1 & INFO/DP >= {params.depth} & INFO/AF >= {params.vaf} & INFO/AC >= {params.alt}' {input} 2> {log}| "
		"grep -v -f {params.excl} 2>> {log}| "
		"bcftools sort -Oz -o {output.vcf} 2>> {log} && "
		"tabix -p vcf {output.vcf}"

rule freebayes_merge_control_vars:
	input:
		vcf=expand(f"results/{{sample}}_control/{{sample}}.freebayes.filt.vcf.gz", sample=config['controls'].values()),
		tbi=expand(f"results/{{sample}}_control/{{sample}}.freebayes.filt.vcf.gz.tbi", sample=config['controls'].values())
	output:
		vcf="results/customPoN.freebayes.vcf.gz",
		tbi="results/customPoN.freebayes.vcf.gz.tbi"
	message:
		"freebayes merge controls variants"
	benchmark:
		"benchmarks/freebayesMergeControl.txt"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/freebayesMergeControl.log"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2> {log} && "
		"tabix -p vcf {output.vcf} 2>> {log}"


rule freebayes_remove_control_vars: 
	input:
		vcf="results/{sample}_tumor/{sample}.freebayes.filt.vcf.gz",
		PoN="results/customPoN.freebayes.vcf.gz",
		tbi="results/customPoN.freebayes.vcf.gz.tbi"
	output:
		vcf="results/{sample}_tumor/{sample}.freebayes.filtOnCtr.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.freebayes.filtOnCtr.vcf.gz.tbi"
	message:
		"freebayes remove controls variants"
	benchmark:
		"benchmarks/{sample}.freebayesRemoveControl.txt"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/{sample}.freebayesRemoveControl.log"
	shell:
		"bcftools isec -w1 -Oz -c none -n~10 -o {output.vcf} {input.vcf} {input.PoN} 2> {log} && "
		"tabix -p vcf {output.vcf}"


rule merge_freebayes_samples_variants: 
    input:
        vcf=expand("results/{sample}_tumor/{sample}.freebayes.filtOnCtr.vcf.gz", sample=config['samples'].values()),
        tbi=expand("results/{sample}_tumor/{sample}.freebayes.filtOnCtr.vcf.gz.tbi", sample=config['samples'].values())
    output:
        vcf="results/multisample.freebayes.vcf.gz",
        tbi="results/multisample.freebayes.vcf.gz.tbi"
    threads: 1
    log: 
        "logs/multisample.freebayes.log"
    shell:
        "bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2>{log} && "
        "tabix -p vcf {output.vcf} 2>>{log}"
