rule freebayesCallTumor:
	input:
		alns="alignments/{sample}.tumor.dd.rec.bam",
		ref=config["genome"],
		regions=config["intervals"],
	output:
		vcf="results/{sample}_tumor/{sample}.freebayes.vcf.gz",
	message:
		"freebayes calling on tumor {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.freebayesCallTumor.txt"
	log:
		"logs/{sample}.freebayesCallTumor.log",
	params:
		chunksize=100000,
		normalize="-m - -f "+config["genome"],
		extra="-F 0.01 -C 2 --pooled-continuous"
	threads: config["threads"]
	resources:
		mem_mb=4096,
	wrapper:
		"v3.3.3/bio/freebayes"

rule freebayesCallControl:
	input:
		alns="alignments/{sample}.control.dd.rec.bam",
		ref=config["genome"],
		regions=config["intervals"],
	output:
		vcf="results/{sample}_control/{sample}.freebayes.vcf.gz",
	message:
		"freebayes calling on control {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.freebayesCallControl.txt"
	log:
		"logs/{sample}.freebayesCallControl.log",
	params:
		chunksize=100000,
		normalize="-m - -f "+config["genome"],
		extra="-F 0.01 -C 2 --pooled-continuous"
	threads: config["threads"]
	resources:
		mem_mb=4096,
	wrapper:
		"v3.3.3/bio/freebayes"


rule freebayesAddAFtumor:
	input:
		"results/{sample}_tumor/{sample}.freebayes.vcf.gz"
	output:
		temp("results/{sample}_tumor/{sample}.freebayes.annotAF.vcf.gz")
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

rule freebayesAddAFcontrol:
	input:
		"results/{sample}_control/{sample}.freebayes.vcf.gz"
	output:
		temp("results/{sample}_control/{sample}.freebayes.annotAF.vcf.gz")
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
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%RO\t%AO\n' {input} 2> {log} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {params} 2>> {log} && 
		tabix -b2 -e2 {params} && bcftools annotate -a {params} --columns CHROM,POS,REF,ALT,AF {input} -Oz -o {output} 2>> {log} && 
		rm {params}
		"""

rule freebayesFilterTumor: 
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

rule freebayesFilterControl: 
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

rule freebayesMergeControl:
	input:
		expand(f"results/{{sample}}_control/{{sample}}.freebayes.filt.vcf.gz", sample=config['controls'].values()),
		expand(f"results/{{sample}}_control/{{sample}}.freebayes.filt.vcf.gz.tbi", sample=config['controls'].values())
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
		"bcftools merge -m none -Oz -o {output.vcf} {input} 2> {log} && "
		"tabix -p vcf {output.vcf} 2>> {log}"


rule freebayesRemoveControl: 
	input:
		vcf="results/{sample}_tumor/{sample}.freebayes.filt.vcf.gz",
		PoN="results/customPoN.freebayes.vcf.gz",
		tbi="results/customPoN.freebayes.vcf.gz.tbi"
	output:
		temp("results/{sample}_tumor/{sample}.freebayes.filtOnCtr.vcf.gz")
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
		"bcftools isec -w1 -Oz -c none -n~10 -o {output} {input.vcf} {input.PoN} 2> {log}"


#rule freebayes_vep_annotation: 
#	input:
#		vcf="results/freebayes/{sample}.freebayes.filtOnCtr.vcf.gz"
#	output:
#		vcf="results/freebayes/{sample}.freebayes.filtOnCtr.vep.vcf.gz",
#		summ_html="results/freebayes/{sample}.freebayes.vep.summary.html"
#	threads: 10
#	params:
#		vep="/home/alessio/programs/ensembl-vep/vep",
#		ref=config["genome"]
#	log: 
#		"logs/{sample}.freebayes.filtOnCtr.vep.log"
#	shell:
#		"perl {params.vep} -i {input.vcf} -o {output.vcf} --everything --species homo_sapiens --vcf "
#			"--compress_output gzip --stats_file {output.summ_html} --assembly GRCh38 --no_progress --buffer_size 5000 " 
#			"--ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot "
#			"--cache --cache_version 109 --merged --tsl --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number "
#			"--no_escape --xref_refseq --failed 1 --pick --dir_cache /home/alessio/.vep --fasta {params.ref} --format vcf --offline "
#			"--pubmed --fork {threads} --sift b --polyphen b --af --af_1kg --af_gnomade --af_gnomadg --regulatory"

#rule freebayes_splitvep:
#	input:
#		vcf="results/freebayes/{sample}.freebayes.filtOnCtr.vep.vcf.gz"
#	output:
#		tsv="results/freebayes/{sample}.freebayes.filtOnCtr.vep.tsv"
#	threads: 1
#	log: 
#		"logs/{sample}.freebayes.splitvep.tsv"
#	shell:
#		"""
#		bcftools +split-vep -o {output.tsv} {input.vcf} -f "%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%IMPACT\t%SYMBOL\t%Gene\t%Feature_type\t%Feature\t%BIOTYPE\t%EXON\t%INTRON\t%HGVSc\t%HGVSp\t%CDS_position\t%Amino_acids\t%Codons\t%Existing_variation\t%DISTANCE\t%VARIANT_CLASS\t%CANONICAL\t%MANE_SELECT\t%MANE_PLUS_CLINICAL\t%TSL\t%ENSP\t%UNIPROT_ISOFORM\t%REFSEQ_MATCH\t%GENE_PHENO\t%SIFT\t%PolyPhen\t%DOMAINS\t%AF\t%EUR_AF\t%gnomADe_AF\t%gnomADe_NFE_AF\t%gnomADg_AF\t%gnomADg_NFE_AF\t%MAX_AF\t%CLIN_SIG\t%SOMATIC\t%PHENO\t%PUBMED\t%MOTIF_NAME\t%HIGH_INF_POS\t%MOTIF_SCORE_CHANGE\t%TRANSCRIPTION_FACTORS[\t%AD][\t%AF][\t%GT]\n" -d -A tab 
#		"""



#rule merge_freebayes_samples_variants: 
#    input:
#        vcf=expand("results/freebayes/{sample}.freebayes.filtOnCtr.vcf.gz", sample=config['samples'].values()),
#        tbi=expand("results/freebayes/{sample}.freebayes.filtOnCtr.vcf.gz.tbi", sample=config['samples'].values())
#    output:
#        vcf="results/freebayes/multisample.freebayes.vcf.gz",
#        tbi="results/freebayes/multisample.freebayes.vcf.gz.tbi"
#    threads: 1
#    log: 
#        "logs/multisample.freebayes.log"
#    shell:
#        "bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2> {log} &&"
#        "tabix -p vcf {output.vcf}"
