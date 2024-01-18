# Processing control samples

rule bcftools_call_on_controls:  
	input:
		bam="data/bam/{ctr}.dd.rec.bam",
		bai="data/bam/{ctr}.dd.rec.bai"
	output:
		"results/{{ctr}}_control/{ctr}.bcftools.vcf.gz"
	threads: 1
	params:
		ref=config["genome"],
		intervals=config["intervals"]
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{ctr}.ctr.bcftools.call.log"
	shell:
		"bcftools mpileup -Ou -d 10000 -R {params.intervals} -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR -f {params.ref} {input.bam} | "
		"bcftools call -mv | "
		"bcftools norm -m - -f {params.ref} -Oz -o {output} 2> {log}"


rule bcftools_add_af_field_on_controls:  
	input:
		"results/{{ctr}}_control/{ctr}.bcftools.vcf.gz"
	output:
		vcf=temp("results/{{ctr}}_control/{ctr}.bcftools.addaf.vcf.gz"),
		tsv=temp("results/{{ctr}}_control/{ctr}.annot.tsv.gz"),
		tsv_tbi=temp("results/{{ctr}}_control/{ctr}.annot.tsv.gz.tbi")
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{ctr}.ctr.bcftools.addaf.log"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%AD{{0}}\t%AD{{1}}]' {input} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {output.tsv} && 
		tabix -b2 -e2 {output.tsv} && 
		bcftools annotate -a {output.tsv} -h <(echo '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,AF {input} -Oz -o {output.vcf}
		"""


rule bcftools_filter_controls: 
	input:
		"results/{{ctr}}_control/{ctr}.bcftools.addaf.vcf.gz"
	output:
		vcf=temp("results/{{ctr}}_control/{ctr}.bcftools.filt.vcf.gz"),
		tbi=temp("results/{{ctr}}_control/{ctr}.bcftools.filt.vcf.gz.tbi")
	threads: 1
	params:
		excl=config["chr_to_exclude"],
		depth=config['filtering_controls']['min_depth'],
		vaf=config['filtering_controls']['vaf'],
		alt=config['filtering_controls']['alt_depth']		
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{ctr}.bcftools.filter.log"
	shell:
		"bcftools view -i'FORMAT/DP >= {params.depth} & INFO/AF >= {params.vaf} & FORMAT/AD[0:1] >= {params.alt}' {input} | "
		"grep -v -f {params.excl} | "
		"bcftools sort -Oz -o {output.vcf} 2> {log} && "
		"tabix -p vcf {output.vcf}"


rule bcftools_merge_bcftools_controls:
	input:
		vcf=expand(f"results/{{ctr}}_control/{{ctr}}.bcftools.filt.vcf.gz", ctr=config['controls'].values()),
		tbi=expand(f"results/{{ctr}}_control/{{ctr}}.bcftools.filt.vcf.gz.tbi", ctr=config['controls'].values())
	output:
		vcf="results/{{ctr}}_control/custom_pon.bcftools.vcf.gz",
		tbi="results/{{ctr}}_control/custom_pon.bcftools.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/bcftools.merge.normals.vars.log"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2> {log} && "
		"tabix -p vcf {output.vcf}"


###########################################################################################


rule bcftools_call_on_tumors:  
	input:
		bam="data/bam/{tum}.dd.rec.bam",
		bai="data/bam/{tum}.dd.rec.bai"
	output:
		"results/{tum}/{tum}.bcftools.vcf.gz"
	threads: 1
	params:
		ref=config["genome"],
		intervals=config["intervals"]
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{tum}.bcftools.call.log"
	shell:
		"bcftools mpileup -Ou -d 10000 -R {params.intervals} -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR -f {params.ref} {input.bam} | "
		"bcftools call -mv | "
		"bcftools norm -m - -f {params.ref} -Oz -o {output} 2> {log}"


rule bcftools_add_af_field_on_tumors:  
	input:
		"results/{tum}/{tum}.bcftools.vcf.gz"
	output:
		vcf=temp("results/{tum}/{tum}.bcftools.addaf.vcf.gz"),
		tsv=temp("results/{tum}/{tum}.annot.tsv.gz"),
		tsv_tbi=temp("results/{tum}/{tum}.annot.tsv.gz.tbi")
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{tum}.tum.bcftools.addaf.log"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%AD{{0}}\t%AD{{1}}]' {input} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {output.tsv} && 
		tabix -b2 -e2 {output.tsv} && 
		bcftools annotate -a {output.tsv} -h <(echo '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,AF {input} -Oz -o {output.vcf}
		"""


rule bcftools_filter_tumors: 
	input:
		"results/{tum}/{tum}.bcftools.addaf.vcf.gz"
	output:
		vcf=temp("results/{tum}/{tum}.bcftools.filt.vcf.gz"),
		tbi=temp("results/{tum}/{tum}.bcftools.filt.vcf.gz.tbi")
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
		"logs/{tum}.tum.filt.bcftools.log"
	shell:
		"bcftools view -i 'FORMAT/DP >= {params.depth} & INFO/AF >= {params.vaf} & FORMAT/AD[0:1] >= {params.alt}' {input} | "
		"grep -v -f {params.excl} | "
		"bcftools sort -Oz -o {output.vcf} 2> {log} && "
		"tabix -p vcf {output.vcf}"	


rule bcftools_remove_controls: 
	input:
		vcf="results/{tum}/{tum}.bcftools.filt.vcf.gz",
		vcf_tbi="results/{tum}/{tum}.bcftools.filt.vcf.gz.tbi",
		custom_pon="results/controls/custom_pon.bcftools.vcf.gz",
		pon_tbi="results/controls/custom_pon.bcftools.vcf.gz.tbi"
	output:
		"results/{tum}/{tum}.bcftools.filtOnCtr.vcf.gz"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/{tum}.filtOnCtr.vars.log"
	shell:
		"bcftools isec -w1 -Oz -c none -n~10 -o {output} {input.vcf} {input.custom_pon} 2> {log}"


#rule bcftools_vep_annotation: 
#	input:
#		"results/variants/{sample}.bt.filtOnCtr.vcf.gz"
#	output:
#		vcf="results/variants/bcftools/{sample}.bt.filtOnCtr.vep.vcf.gz",
#		summ_html="results/bcftools/{sample}.vep.summary.html"
#	threads: config['ncores']
#	params:
#		vep="/home/alessio/programs/ensembl-vep/vep",
#		ref=config["genome"]
#	conda:
#		"../envs/bcftools.yaml"
#	log: 
#		"logs/bcftools/{sample}.bt.filtOnCtr.vep.log"
#	shell:
#		"perl {params.vep} -i {input} -o {output.vcf} --everything --species homo_sapiens --vcf "
#			"--compress_output gzip --stats_file {output.summ_html} --assembly GRCh38 --no_progress --buffer_size 5000 " 
#			"--ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot "
#			"--cache --cache_version 110 --merged --tsl --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number "
#			"--no_escape --xref_refseq --failed 1 --pick --dir_cache /home/alessio/.vep --fasta {params.ref} --format vcf --offline "
#			"--pubmed --fork {threads} --sift b --polyphen b --af --af_1kg --af_gnomade --af_gnomadg --regulatory"


#rule bcftools_splitvep:
#	input:
#		"results/bcftools/{sample}.bt.filtOnCtr.vep.vcf.gz"
#	output:
#		"results/bcftools/{sample}.bt.filtOnCtr.vep.tsv"
#	threads: 1
#	conda:
#		"../envs/bcftools.yaml"
#	log: 
#		"logs/bcftools/{sample}.bt.splitvep.log"
#	shell:
#		"""
#		bcftools +split-vep -o {output} {input} -f"%CHROM\t%POS\t%REF\t%ALT\t%IMPACT\t%Consequence\t%SYMBOL\t%EUR_AF\t%gnomADe_AF\t%gnomADe_NFE_AF\t%gnomADg_AF\t%gnomADg_NFE_AF\t%Existing_variation\t%CLIN_SIG\t%SIFT\t%PolyPhen[\t%AF][\t%AD]\n" -d -A tab 
#		"""
#
#rule bcftools_db_algorithm:
#	input:
#		"results/bcftools/{sample}.sample.bcftools.filtOnCtr.vep.tsv" 
#	output:
#		"results/bcftools/{sample}.sample.bcftools.filtOnDB.tsv"
#	threads: 1 
#	params:
#		"workflow/script/"
#	shell:
#		"{params} {input} {output}"



#rule merge_bcftools_samples_variants: 
#   input:
#       vcf=expand("results/bcftools/{sample}.sample.bcftools.filtOnCtr.vcf.gz", sample=config['samples'].values()),
#       tbi=expand("results/bcftools/{sample}.sample.bcftools.filtOnCtr.vcf.gz.tbi", sample=config['samples'].values())
#    output:
#        vcf="results/bcftools/multisample.bcftools.vcf.gz",
#        tbi="results/bcftools/multisample.bcftools.vcf.gz.tbi"
#    threads: 1
#    log: 
#        "logs/multisample.bcftools.log"
#    shell:
#        "bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2> {log} &&"
#        "tabix -p vcf {output.vcf}"
