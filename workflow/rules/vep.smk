rule get_vep_cache:
	output:
		directory("resources/vep/cache"),
	params:
		species="homo_sapiens_merged",
		build="GRCh38",
		release="110",
	log:
		"logs/vep_cache.log",
	cache: "omit-software"  # save space and time with between workflow caching (see docs)
	wrapper:
		"v3.3.5/bio/vep/cache"


rule download_vep_plugins:
	output:
		directory("resources/vep/plugins")
	params:
		release=110
	wrapper:
		"v3.3.5/bio/vep/plugins"

################################################# VEP ####################################################

rule vep_single_sample:
	input:
		vcf="results/{sample}_tumor/{sample}.{caller}.filtOnCtr.vcf.gz",
		cache=config["vepcache"]
	output:
		vcf="results/{sample}_tumor/{sample}.{caller}.filtOnCtr.vep.vcf.gz",  
		tbi="results/{sample}_tumor/{sample}.{caller}.filtOnCtr.vep.vcf.gz.tbi",
		html="results/{sample}_tumor/{sample}.{caller}.filtOnCtr.vep.html"
	threads: 
		config["threads"]
	log:
		"logs/{sample}.vep_{caller}_single_sample.log"
	conda:
		"../envs/vep.yaml"
	params:
		ref=config["genome"],
		clinvar=config["clinvar"],
		dbNSFP=config["dbNSFP"],
	shell:
		"vep -i {input.vcf} -o {output.vcf} --fork {threads} --compress_output bgzip --everything --offline --species homo_sapiens --stats_file {output.html} --assembly GRCh38 --cache --dir_cache {input.cache} --cache_version 110 --merged --fasta {params.ref} --format vcf --symbol --no_intergenic --merged --cache --pick --pick_allele --vcf --plugin dbNSFP,{params.dbNSFP},SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,FATHMM_converted_rankscore,FATHMM_pred --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN 2>{log} && tabix {output.vcf} 2>>{log}"


rule vep_multisample:
	input:
		vcf="results/multisample.{caller}.vcf.gz",
		cache=config["vepcache"]
	output:
		vcf="results/multisample.{caller}.vep.vcf.gz",
		tbi="results/multisample.{caller}.vep.vcf.gz.tbi",
		html="results/multisample.{caller}.vep.html"
	threads:
		config["threads"]
	log:
		"logs/vep_{caller}_multisample.log"
	conda:
		"../envs/vep.yaml"
	params:
		ref=config["genome"],
		clinvar=config["clinvar"],
		dbNSFP=config["dbNSFP"],
	shell:
		"vep -i {input.vcf} -o {output.vcf} --fork {threads} --compress_output bgzip --everything --offline --species homo_sapiens --stats_file {output.html} --assembly GRCh38 --cache --dir_cache {input.cache} --cache_version 110 --merged --fasta {params.ref} --format vcf --symbol --no_intergenic --merged --cache --pick --pick_allele --vcf --plugin dbNSFP,{params.dbNSFP},SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,FATHMM_converted_rankscore,FATHMM_pred --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN 2>{log} && tabix {output.vcf} 2>>{log}"

########################################## SPLIT VEP #########################################################

#MUTECT2

rule split_vep_mutect2_single_sample:
	input:
		vcf="results/{sample}_tumor/{sample}.mutect2.filtOnCtr.vep.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.mutect2.filtOnCtr.vep.vcf.gz.tbi",
		html="results/{sample}_tumor/{sample}.mutect2.filtOnCtr.vep.html",
	output:
		"results/{sample}_tumor/{sample}.mutect2.filtOnCtr.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/{sample}.split_vep_mutect2_single_sample.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule split_vep_mutect2_multisample:
	input:
		vcf="results/multisample.mutect2.vep.vcf.gz",
		tbi="results/multisample.mutect2.vep.vcf.gz.tbi",
		html="results/multisample.mutect2.vep.html"
	output:
		"results/multisample.mutect2.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/split_vep_mutect2_multisample.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

# FREEBAYES

rule split_vep_freebayes_single_sample:
	input:
		vcf="results/{sample}_tumor/{sample}.freebayes.filtOnCtr.vep.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.freebayes.filtOnCtr.vep.vcf.gz.tbi",
		html="results/{sample}_tumor/{sample}.freebayes.filtOnCtr.vep.html",
	output:
		"results/{sample}_tumor/{sample}.freebayes.filtOnCtr.vep.tmp01.tsv"
	threads:1
	log:
		"logs/{sample}.split_vep_freebayes_single_sample.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%RO\t%AO\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule split_vep_freebayes_multisample:
	input:
		vcf="results/multisample.freebayes.vep.vcf.gz",
		tbi="results/multisample.freebayes.vep.vcf.gz.tbi",
		html="results/multisample.freebayes.vep.html"
	output:
		"results/multisample.freebayes.vep.tmp01.tsv"
	threads: 1 
	log:
		"logs/split_vep_freebayes_multisample.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%RO\t%AO\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

# BCFTOOLS 

rule split_vep_bcftools_single_sample:
	input:
		vcf="results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vep.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vep.vcf.gz.tbi",
		html="results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vep.html",
	output:
		"results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/{sample}.split_vep_bcftools_single_sample.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule split_vep_bcftools_multisample:
	input:
		vcf="results/multisample.bcftools.vep.vcf.gz",
		tbi="results/multisample.bcftools.vep.vcf.gz.tbi",
		html="results/multisample.bcftools.vep.html"
	output:
		"results/multisample.bcftools.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/split_vep_bcftools_multisample.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

# VARSCAN 

rule split_vep_varscan_single_sample:
	input:
		vcf="results/{sample}_tumor/{sample}.varscan.filtOnCtr.vep.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.varscan.filtOnCtr.vep.vcf.gz.tbi",
		html="results/{sample}_tumor/{sample}.varscan.filtOnCtr.vep.html",
	output:
		"results/{sample}_tumor/{sample}.varscan.filtOnCtr.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/{sample}.split_vep_varscan_single_sample.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%RD\t%AD\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule split_vep_varscan_multisample:
	input:
		vcf="results/multisample.varscan.vep.vcf.gz",
		tbi="results/multisample.varscan.vep.vcf.gz.tbi",
		html="results/multisample.varscan.vep.html"
	output:
		"results/multisample.varscan.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/split_vep_varscan_multisample.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%RD\t%AD\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""
