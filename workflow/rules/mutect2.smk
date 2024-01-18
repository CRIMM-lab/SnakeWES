rule Mutect2Tumor: 
	input:
		fasta=config["genome"],
		map="alignments/{sample}.tumor.dd.rec.bam",
		intervals=config["intervals"],
		pon=config["gatk_pon"],
		germline=config["gnomAD"]
	output:
		vcf="results/{sample}_tumor/{sample}.mutect2.vcf.gz",
		bam="alignments/{sample}.tumor.mutect2.bam",        
		f1r2="data/{sample}.tumor.f1r2.tar.gz"
	message:
		"Mutect2 calling with {wildcards.sample}"
	threads: 1
	resources:
		mem_mb=4096
	log:
		"logs/{sample}.Mutect2Tumor.log",
	params:
		extra="--max-reads-per-alignment-start 0" 
	wrapper:
		"v3.3.3/bio/gatk/mutect"

rule Mutect2Control: 
	input:
		fasta=config["genome"],
		map="alignments/{sample}.control.dd.rec.bam",
		intervals=config["intervals"],
		pon=config["gatk_pon"],
		germline=config["gnomAD"]
	output:
		vcf="results/{sample}_control/{sample}.mutect2.vcf.gz",
		bam="alignments/{sample}.control.mutect2.bam",        
		f1r2="data/{sample}.control.f1r2.tar.gz"
	message:
		"Mutect2 calling with {wildcards.sample}"
	threads: 1
	resources:
		mem_mb=4096
	log:
		"logs/{sample}.Mutect2Control.log",
	params:
		extra="--max-reads-per-alignment-start 0" 
	wrapper:
		"v3.3.3/bio/gatk/mutect"

rule OrientationModelTumor:
	input:
		f1r2="data/{sample}.tumor.f1r2.tar.gz",
	output:
		"data/{sample}.tumor.artifacts_prior.tar.gz",
	resources:
		mem_mb=4096,
	log:
		"logs/{sample}.OrientationModelTumor.log",
	wrapper:
		"v3.3.3/bio/gatk/learnreadorientationmodel"

rule OrientationModelControl:
	input:
		f1r2="data/{sample}.control.f1r2.tar.gz",
	output:
		"data/{sample}.control.artifacts_prior.tar.gz",
	resources:
		mem_mb=4096,
	log:
		"logs/{sample}.OrientationModelControl.log",
	wrapper:
		"v3.3.3/bio/gatk/learnreadorientationmodel"

rule GetPileupSummariesTumor:
	input:
		bam="alignments/{sample}.tumor.dd.rec.bam",
		intervals=config["intervals"],
		variants=config["gnomAD"]
	output:
		"data/{sample}.tumor.pileup.table"
	threads: 1
	resources:
		mem_mb=4096,
	params:
		extra="",
	log:
		"logs/{sample}.GetPileupSummariesTumor.log",
	wrapper:
		"v3.3.3/bio/gatk/getpileupsummaries"

rule GetPileupSummariesControl:
	input:
		bam="alignments/{sample}.control.dd.rec.bam",
		intervals=config["intervals"],
		variants=config["gnomAD"]
	output:
		"data/{sample}.control.pileup.table"
	threads: 1
	resources:
		mem_mb=4096,
	params:
		extra="",
	log:
		"logs/{sample}.GetPileupSummariesControl.log",
	wrapper:
		"v3.3.3/bio/gatk/getpileupsummaries"

rule CalculateContaminationTumor:
	input:
		"data/{sample}.tumor.pileup.table"
	output:
		contamination="data/{sample}.tumor.contamination.table",
		segmentation="data/{sample}.tumor.tumorseg.txt"
	threads: 1
	log:
		"logs/{sample}.CalculateContaminationTumor.log",
	conda:
		"../envs/gatk4.yaml"
	params:
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"]),
		mem_mb="-Xmx4G"
		#extra="--tumor-segmentation data/{wildcards.sample}.tumseg.txt"
	shell:
		"gatk --java-options {params.mem_mb} CalculateContamination -I {input} -O {output.contamination} --tumor-segmentation {output.segmentation} > {log} 2>&1"

rule CalculateContaminationControl:
	input:
		"data/{sample}.control.pileup.table"
	output:
		contamination="data/{sample}.control.contamination.table",
		segmentation="data/{sample}.control.tumorseg.txt"
	threads: 1
	log:
		"logs/{sample}.CalculateContaminationControl.log",
	conda:
		"../envs/gatk4.yaml"
	params:
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"]),
		mem_mb="-Xmx4G"
		#extra="--tumor-segmentation data/{wildcards.sample}.tumseg.txt"
	shell:
		"gatk --java-options {params.mem_mb} CalculateContamination -I {input} -O {output.contamination} --tumor-segmentation {output.segmentation} > {log} 2>&1"

rule FilterMutectCallsTumor:
	input:
		vcf="results/{sample}_tumor/{sample}.mutect2.vcf.gz",
		ref=config["genome"],
		bam="alignments/{sample}.tumor.dd.rec.bam",
		intervals=config["intervals"],
		contamination="data/{sample}.tumor.contamination.table", # from gatk CalculateContamination
		segmentation="data/{sample}.tumor.tumorseg.txt", # from gatk CalculateContamination
		f1r2="data/{sample}.tumor.artifacts_prior.tar.gz" # from gatk LearnReadOrientationBias
	output:
		vcf="results/{sample}_tumor/{sample}.mutect2.filtered.vcf.gz",
	log:
		"logs/{sample}.FilterMutectCallsTumor.log",
	params:
		#extra="--tumor-segmentation data/{wildcard.sample}.tumseg.txt",  # optional arguments, see GATK docs
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"])  # optional
	resources:
		mem_mb=4096,
	wrapper:
		"v3.3.3/bio/gatk/filtermutectcalls"

rule FilterMutectCallsControl:
	input:
		vcf="results/{sample}_control/{sample}.mutect2.vcf.gz",
		ref=config["genome"],
		bam="alignments/{sample}.control.dd.rec.bam",
		intervals=config["intervals"],
		contamination="data/{sample}.control.contamination.table", # from gatk CalculateContamination
		segmentation="data/{sample}.control.tumorseg.txt", # from gatk CalculateContamination
		f1r2="data/{sample}.control.artifacts_prior.tar.gz" # from gatk LearnReadOrientationBias
	output:
		vcf="results/{sample}_control/{sample}.mutect2.filtered.vcf.gz",
	log:
		"logs/{sample}.FilterMutectCallsControl.log",
	params:
		#extra="--tumor-segmentation data/{wildcard.sample}.tumseg.txt",  # optional arguments, see GATK docs
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"])  # optional
	resources:
		mem_mb=4096,
	wrapper:
		"v3.3.3/bio/gatk/filtermutectcalls"
