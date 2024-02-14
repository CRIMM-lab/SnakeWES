rule Mutect2TumorPureCN: 
	input:
		fasta=config["genome"],
		map="alignments/{sample}.tumor.dd.rec.bam",
		intervals=config["intervals"],
		pon=config["gatk_pon"],
		germline=config["gnomAD"]
	output:
		vcf="results/{sample}_tumor/{sample}.mutect2.pureCN.vcf.gz",
		#bam="alignments/{sample}.tumor.mutect2.bam",        
		#f1r2="data/{sample}.tumor.f1r2.tar.gz"
	message:
		"Mutect2 calling for PureCN - tumor {wildcards.sample}"
	threads: config["threads"]
	resources:
		mem_mb=4096
	log:
		"logs/{sample}.Mutect2Tumor.pureCN.log",
	params:
		extra="--max-reads-per-alignment-start 0 --genotype-germline-sites true --genotype-pon-sites true --interval-padding 75"
	wrapper:
		"v3.3.3/bio/gatk/mutect"


rule Mutect2ControlPureCN: 
	input:
		fasta=config["genome"],
		map="alignments/{sample}.control.dd.rec.bam",
		intervals=config["intervals"],
		pon=config["gatk_pon"],
		germline=config["gnomAD"]
	output:
		vcf="results/{sample}_control/{sample}.mutect2.pureCN.vcf.gz",
		#bam="alignments/{sample}.control.mutect2.bam",        
		#f1r2="data/{sample}.control.f1r2.tar.gz"
	message:
		"Mutect2 calling for PureCN - control {wildcards.sample}"
	threads: config["threads"]
	resources:
		mem_mb=4096
	log:
		"logs/{sample}.Mutect2Control.pureCN.log",
	params:
		extra="--max-reads-per-alignment-start 0 --genotype-germline-sites true --genotype-pon-sites true --interval-padding 75"
	wrapper:
		"v3.3.3/bio/gatk/mutect"


rule Mutect2PoNpureCN:
	input:
		expand(f"results/{{sample}}_control/{{sample}}.mutect2.pureCN.vcf.gz", sample=config["controls"].values()),
	output:
		vcf="results/customPoN.PureCN.vcf.gz",
		tbi="results/customPoN.PureCN.vcf.gz.tbi"
	conda:
		"../envs/bcftools.yaml"
	params:
		alt=config["filtering_controls"]["alt_depth"],
		min_depth=config["filtering_controls"]["min_depth"],
		vaf=config["filtering_controls"]["vaf"],
	log:
		"logs/Mutect2PoNpureCN.log"
	shell:
		"bcftools merge -m none {input} | "
		"bcftools view -Ov -i'FORMAT/DP >= {params.min_depth} & FORMAT/AD >= {params.alt} & FORMAT/AF >= {params.vaf}' | "
		"bcftools sort -Oz -o {output} && "
		"tabix -p vcf {output}"


rule IntervalFilePureCN:
	input:
		bed="resources/SureSelectV8.wes.mod.bed",
		refgen=config["genome"],
		mappability="resources/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw"
	output:
		txt="resources/baits_hg38_intervals.txt",
	log:
		"logs/pureCN",
	conda:
		"../envs/pureCN.yaml"
	params:
		"workflow/scripts/interval_pureCN.R"
	shell:
		"Rscript {params} {input.bed} {input.refgen} {input.mappability} {output.txt}"


rule CoverageTumorPureCN:
	input:
		txt="resources/baits_hg38_intervals.txt",
		bam="alignments/{sample}.tumor.dd.rec.bam",
		bai="alignments/{sample}.tumor.dd.rec.bai"
	output:
		cov_raw="alignments/{sample}.tumor.coveragePureCN.txt",
		cov_loess="alignments/{sample}.tumor.coveragePureCN.loess.txt",
		qc="alignments/{sample}.control.coveragePureCN.loess.qc.txt"
	conda:
		"../envs/pureCN.yaml"
	log:
		"logs/{sample}.CoverageTumorPureCN.log"
	params:
		"workflow/scripts/coverage_pureCN.R"
	shell:
		"Rscript {params} {input.txt} {input.bam} {output.cov_raw} {output.cov_loess} {output.qc}"


rule CoverageControlPureCN:
	input:
		txt="resources/baits_hg38_intervals.txt",
		bam="alignments/{sample}.control.dd.rec.bam",
		bai="alignments/{sample}.control.dd.rec.bai"
	output:
		cov_raw="alignments/{sample}.control.coveragePureCN.txt",
		cov_loess="alignments/{sample}.control.coveragePureCN.loess.txt",
		qc="alignments/{sample}.control.coveragePureCN.loess.qc.txt"
	conda:
		"../envs/pureCN.yaml"
	log:
		"logs/{sample}.CoverageControlPureCN.log"
	params:
		"workflow/scripts/coverage_pureCN.R"
	shell:
		"Rscript {params} {input.txt} {input.bam} {output.cov_raw} {output.cov_loess} {output.qc}"


rule ControlDBpureCN:
	input:
		ctr_cov=expand(f"alignments/{{sample}}.control.coveragePureCN.loess.txt", sample=config["controls"].values()),
		pon="results/customPoN.PureCN.vcf.gz",
		tbi="results/customPoN.PureCN.vcf.gz.tbi"
	output:
		cov_list="alignments/control_coverages.list",
		rds="data/controlDB.rds",
		rds_bias="data/mapping_bias.rds"
	conda:
		"../envs/pureCN.yaml"
	params:
		"workflow/scripts/normalDB_pureCN.R"
	shell:
		"ls {input.ctr_cov} > {output.cov_list} && "
		"Rscript {params} {output.cov_list} {input.pon} {output.rds} {output.rds_bias}"

#rule PureCN:
#	input:
		
#	output:


