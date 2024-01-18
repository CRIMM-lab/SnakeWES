rule fastqcUntrimmedTumor:
	input:
		"data/fastq/{sample}.tumor.{strand}.fastq.gz"
	output:
		html="qc/{sample}_tumor_{strand}_fastqc.html",
		zip="qc/{sample}_tumor_{strand}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	params:
		extra = "--quiet"
	log:
		"logs/{sample}_{strand}.fastqcUntrimmedTumor.log"
	threads: config["threads"]
	resources:
		mem_mb = 1024
	wrapper:
		"v3.3.3/bio/fastqc"

rule fastqcUntrimmedControl:
	input:
		"data/fastq/{sample}.control.{strand}.fastq.gz"
	output:
		html="qc/{sample}_control_{strand}_fastqc.html",
		zip="qc/{sample}_control_{strand}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	params:
		extra = "--quiet"
	log:
		"logs/{sample}_{strand}.fastqcUntrimmedControl.log"
	threads: config["threads"]
	resources:
		mem_mb = 1024
	wrapper:
		"v3.3.3/bio/fastqc"

rule fastqcUntrimmedGermline:
	input:
		"data/fastq/{sample}.germline.{strand}.fastq.gz"
	output:
		html="qc/{sample}_germline_{strand}_fastqc.html",
		zip="qc/{sample}_germline_{strand}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	params:
		extra = "--quiet"
	log:
		"logs/{sample}_{strand}.fastqcUntrimmedGermline.log"
	threads: config["threads"]
	resources:
		mem_mb = 1024
	wrapper:
		"v3.3.3/bio/fastqc"

rule fastqcTrimmedTumor:
	input:
		"data/fastq/{sample}.tumor.{strand}.tr.fastq.gz"
	output:
		html="qc/{sample}_tumor_{strand}_tr_fastqc.html",
		zip="qc/{sample}_tumor_{strand}_tr_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	params:
		extra = "--quiet"
	log:
		"logs/{sample}_{strand}.fastqcUntrimmedTumor.log"
	threads: config["threads"]
	resources:
		mem_mb = 1024
	wrapper:
		"v3.3.3/bio/fastqc"

rule fastqcTrimmedControl:
	input:
		"data/fastq/{sample}.control.{strand}.tr.fastq.gz"
	output:
		html="qc/{sample}_control_{strand}_tr_fastqc.html",
		zip="qc/{sample}_control_{strand}_tr_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	params:
		extra = "--quiet"
	log:
		"logs/{sample}_{strand}.fastqcTrimmedControl.log"
	threads: config["threads"]
	resources:
		mem_mb = 1024
	wrapper:
		"v3.3.3/bio/fastqc"

rule fastqcTrimmedGermline:
	input:
		"data/fastq/{sample}.germline.{strand}.tr.fastq.gz"
	output:
		html="qc/{sample}_germline_{strand}_tr_fastqc.html",
		zip="qc/{sample}_germline_{strand}_tr_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	params:
		extra = "--quiet"
	log:
		"logs/{sample}_{strand}.fastqcTrimmedGermline.log"
	threads: config["threads"]
	resources:
		mem_mb = 1024
	wrapper:
		"v3.3.3/bio/fastqc"