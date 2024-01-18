# Varscan

rule varscan_vep_annotation: 
	input:
		vcf="results/varscan/multisample.varscan.vcf.gz",
		tbi="results/varscan/multisample.varscan.vcf.gz.tbi"
	output:
		vcf="results/varscan/multisample.varscan.vep.vcf.gz",
		summ_html="results/varscan/multisample.varscan.vep.summary.html"
	threads: 10
	params:
		vep="/home/alessio/programs/ensembl-vep/vep",
		ref=config["genome"]
	log: 
		"logs/multisample.varscan.vep.log"
	shell:
		"perl {params.vep} -i {input.vcf} -o {output.vcf} --everything --species homo_sapiens --vcf "
			"--compress_output gzip --stats_file {output.summ_html} --assembly GRCh38 --no_progress --buffer_size 5000 " 
			"--ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot "
			"--cache --cache_version 110 --merged --tsl --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number "
			"--no_escape --xref_refseq --failed 1 --pick --dir_cache /home/alessio/.vep --fasta {params.ref} --format vcf --offline "
			"--pubmed --fork {threads} --sift b --polyphen b --af --af_1kg --af_gnomade --af_gnomadg --regulatory"

rule varscan_splitvep:
	input:
		"results/varscan/multisample.varscan.vep.vcf.gz"
	output:
		temp("results/varscan/multisample.varscan.vep.tmp.tsv")
	threads: 1
	log: 
		"logs/multisample.varscan.splitvep.tsv"
	shell:
		"""
		bcftools +split-vep -o {output} {input} -f "%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%IMPACT\t%SYMBOL\t%Gene\t%Feature_type\t%Feature\t%BIOTYPE\t%EXON\t%INTRON\t%HGVSc\t%HGVSp\t%CDS_position\t%Amino_acids\t%Codons\t%Existing_variation\t%DISTANCE\t%VARIANT_CLASS\t%CANONICAL\t%MANE_SELECT\t%MANE_PLUS_CLINICAL\t%TSL\t%ENSP\t%UNIPROT_ISOFORM\t%REFSEQ_MATCH\t%GENE_PHENO\t%SIFT\t%PolyPhen\t%DOMAINS\t%AF\t%EUR_AF\t%gnomADe_AF\t%gnomADe_NFE_AF\t%gnomADg_AF\t%gnomADg_NFE_AF\t%MAX_AF\t%CLIN_SIG\t%SOMATIC\t%PHENO\t%PUBMED\t%MOTIF_NAME\t%HIGH_INF_POS\t%MOTIF_SCORE_CHANGE\t%TRANSCRIPTION_FACTORS[\t%AD][\t%AF][\t%GT]\n" -d -A tab 
		"""

rule varscan_add_sample_name:
	input:
		vcf="results/varscan/multisample.varscan.vcf.gz",
		fi_tsv="results/varscan/multisample.varscan.vep.tmp.tsv"
	output:
		fo_tsv="results/varscan/multisample.varscan.vep.tsv"
	threads: 1
	params: 
		txt="results/varscan/sample_list.txt", 
		script="workflow/scripts/add_sample_names_wes.py"
	log:
	shell:
		"bcftools query -l {input.vcf} > {params.txt} && "
		"python {params.script} {input.fi_tsv} {output.fo_tsv} {params.txt} && "
		"rm {params.txt}"

# freebayes

#rule freebayes_vep_annotation: 
#	input:
#		vcf="results/freebayes/multisample.freebayes.vcf.gz",
#		tbi="results/freebayes/multisample.freebayes.vcf.gz.tbi"
#	output:
#		vcf="results/freebayes/multisample.freebayes.vep.vcf.gz",
#		summ_html="results/freebayes/multisample.freebayes.vep.summary.html"
#	threads: 10
#	params:
#		vep="/home/alessio/programs/ensembl-vep/vep",
#		ref=config["genome"]
#	log: 
#		"logs/multisample.freebayes.vep.log"
#	shell:
#		"perl {params.vep} -i {input.vcf} -o {output.vcf} --everything --species homo_sapiens --vcf "
#			"--compress_output gzip --stats_file {output.summ_html} --assembly GRCh38 --no_progress --buffer_size 5000 " 
#			"--ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot "
#			"--cache --cache_version 109 --merged --tsl --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number "
#			"--no_escape --xref_refseq --failed 1 --pick --dir_cache /home/alessio/.vep --fasta {params.ref} --format vcf --offline "
#			"--pubmed --fork {threads} --sift b --polyphen b --af --af_1kg --af_gnomade --af_gnomadg --regulatory"
#
#rule freebayes_splitvep:
#	input:
#		"results/freebayes/multisample.freebayes.vep.vcf.gz"
#	output:
#		temp("results/freebayes/multisample.freebayes.vep.tmp.tsv")
#	threads: 1
#	log: 
#		"logs/multisample.freebayes.splitvep.tsv"
#	shell:
#		"""
#		bcftools +split-vep -o {output} {input} -f "%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%IMPACT\t%SYMBOL\t%Gene\t%Feature_type\t%Feature\t%BIOTYPE\t%EXON\t%INTRON\t%HGVSc\t%HGVSp\t%CDS_position\t%Amino_acids\t%Codons\t%Existing_variation\t%DISTANCE\t%VARIANT_CLASS\t%CANONICAL\t%MANE_SELECT\t%MANE_PLUS_CLINICAL\t%TSL\t%ENSP\t%UNIPROT_ISOFORM\t%REFSEQ_MATCH\t%GENE_PHENO\t%SIFT\t%PolyPhen\t%DOMAINS\t%AF\t%EUR_AF\t%gnomADe_AF\t%gnomADe_NFE_AF\t%gnomADg_AF\t%gnomADg_NFE_AF\t%MAX_AF\t%CLIN_SIG\t%SOMATIC\t%PHENO\t%PUBMED\t%MOTIF_NAME\t%HIGH_INF_POS\t%MOTIF_SCORE_CHANGE\t%TRANSCRIPTION_FACTORS[\t%AD][\t%AF][\t%GT]\n" -d -A tab 
#		"""

#rule freebayes_add_sample_name:
#	input:
#		vcf="results/freebayes/multisample.freebayes.vcf.gz",
#		fi_tsv="results/freebayes/multisample.freebayes.vep.tmp.tsv"
#	output:
#		fo_tsv="results/freebayes/multisample.freebayes.vep.tsv"
#	threads: 1
#	params: 
#		txt="results/freebayes/sample_list.txt", 
#		script="workflow/scripts/add_sample_names_wes.py"
#	log:
#	shell:
#		"bcftools query -l {input.vcf} > {params.txt} && "
#		"python {params.script} {input.fi_tsv} {output.fo_tsv} {params.txt} && "
#		"rm {params.txt}"
#
# mutect2

rule mutect2_vep_annotation: 
	input:
		vcf="results/mutect2/multisample.mutect2.vcf.gz",
		tbi="results/mutect2/multisample.mutect2.vcf.gz.tbi"
	output:
		vcf="results/mutect2/multisample.mutect2.vep.vcf.gz",
		summ_html="results/mutect2/multisample.mutect2.vep.summary.html"
	threads: 10
	params:
		vep="/home/alessio/programs/ensembl-vep/vep",
		ref=config["genome"]
	log: 
		"logs/multisample.freebayes.vep.log"
	shell:
		"perl {params.vep} -i {input.vcf} -o {output.vcf} --everything --species homo_sapiens --vcf "
			"--compress_output gzip --stats_file {output.summ_html} --assembly GRCh38 --no_progress --buffer_size 5000 " 
			"--ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot "
			"--cache --cache_version 109 --merged --tsl --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number "
			"--no_escape --xref_refseq --failed 1 --pick --dir_cache /home/alessio/.vep --fasta {params.ref} --format vcf --offline "
			"--pubmed --fork {threads} --sift b --polyphen b --af --af_1kg --af_gnomade --af_gnomadg --regulatory"

#rule freebayes_splitvep:
#	input:
#		"results/freebayes/multisample.freebayes.vep.vcf.gz"
#	output:
#		temp("results/freebayes/multisample.freebayes.vep.tmp.tsv")
#	threads: 1
#	log: 
#		"logs/multisample.freebayes.splitvep.tsv"
#	shell:
#		"""
#		bcftools +split-vep -o {output} {input} -f "%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%IMPACT\t%SYMBOL\t%Gene\t%Feature_type\t%Feature\t%BIOTYPE\t%EXON\t%INTRON\t%HGVSc\t%HGVSp\t%CDS_position\t%Amino_acids\t%Codons\t%Existing_variation\t%DISTANCE\t%VARIANT_CLASS\t%CANONICAL\t%MANE_SELECT\t%MANE_PLUS_CLINICAL\t%TSL\t%ENSP\t%UNIPROT_ISOFORM\t%REFSEQ_MATCH\t%GENE_PHENO\t%SIFT\t%PolyPhen\t%DOMAINS\t%AF\t%EUR_AF\t%gnomADe_AF\t%gnomADe_NFE_AF\t%gnomADg_AF\t%gnomADg_NFE_AF\t%MAX_AF\t%CLIN_SIG\t%SOMATIC\t%PHENO\t%PUBMED\t%MOTIF_NAME\t%HIGH_INF_POS\t%MOTIF_SCORE_CHANGE\t%TRANSCRIPTION_FACTORS[\t%AD][\t%AF][\t%GT]\n" -d -A tab 
#		"""
#
#rule freebayes_add_sample_name:
#	input:
#		vcf="results/freebayes/multisample.freebayes.vcf.gz",
#		fi_tsv="results/freebayes/multisample.freebayes.vep.tmp.tsv"
#	output:
#		fo_tsv="results/freebayes/multisample.freebayes.vep.tsv"
#	threads: 1
#	params: 
#		txt="results/freebayes/sample_list.txt", 
#		script="workflow/scripts/add_sample_names_wes.py"
#	log:
#	shell:
#		"bcftools query -l {input.vcf} > {params.txt} && "
#		"python {params.script} {input.fi_tsv} {output.fo_tsv} {params.txt} && "
#		"rm {params.txt}"
#
