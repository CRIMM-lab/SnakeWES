rule fastpTumor: ## aggiungere flag per i threads
	input:
		sample=["data/fastq/{sample}.tumor.R1.fastq.gz", "data/fastq/{sample}.tumor.R2.fastq.gz"]
	output:
		trimmed=["data/fastq/{sample}.tumor.R1.tr.fastq.gz", "data/fastq/{sample}.tumor.R2.tr.fastq.gz"],
		json="data/{sample}.tumor.json",
		failed="data/{sample}.tumor.failedreads.txt",
		html="data/{sample}.tumor.html",
		unpaired1="data/{sample}.tumor.u1.fq.gz",
		unpaired2="data/{sample}.tumor.u2.fq.gz"
	threads: 1
	log:
		"logs/{sample}.fastpTumor.log"
	message:
		"Trimming with fastp"	
	params:
		#adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
		extra="--length_required 30 --detect_adapter_for_pe --disable_quality_filtering"
	benchmark:
		"benchmarks/{sample}.fastpTumor.txt"
	wrapper:
		"v3.3.3/bio/fastp"

rule fastpControl: ## aggiungere flag per i threads
	input:
		sample=["data/fastq/{sample}.control.R1.fastq.gz", "data/fastq/{sample}.control.R2.fastq.gz"]
	output:
		trimmed=["data/fastq/{sample}.control.R1.tr.fastq.gz", "data/fastq/{sample}.control.R2.tr.fastq.gz"],
		json="data/{sample}.control.json",
		failed="data/{sample}.control.failedreads.txt",
		html="data/{sample}.control.html",
		unpaired1="data/{sample}.control.u1.fq.gz",
		unpaired2="data/{sample}.control.u2.fq.gz"
	threads: 1
	log:
		"logs/{sample}.fastpControl.log"
	message:
		"Trimming with fastp"	
	params:
		#adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
		extra="--length_required 30 --detect_adapter_for_pe --disable_quality_filtering"
	benchmark:
		"benchmarks/{sample}.fastpControl.txt"
	wrapper:
		"v3.3.3/bio/fastp"

rule fastpGermline: ## aggiungere flag per i threads
	input:
		sample=["data/fastq/{sample}.germline.R1.fastq.gz", "data/fastq/{sample}.germline.R2.fastq.gz"]
	output:
		trimmed=["data/fastq/{sample}.germline.R1.tr.fastq.gz", "data/fastq/{sample}.germline.R2.tr.fastq.gz"],
		json="data/{sample}.germline.json",
		failed="data/{sample}.germline.failedreads.txt",
		html="data/{sample}.germline.html",
		unpaired1="data/{sample}.germline.u1.fq.gz",
		unpaired2="data/{sample}.germline.u2.fq.gz"
	threads: 1
	log:
		"logs/{sample}.fastpGermline.log"
	message:
		"Trimming with fastp"	
	params:
		#adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
		extra="--length_required 30 --detect_adapter_for_pe --disable_quality_filtering"
	benchmark:
		"benchmarks/{sample}.fastpGermline.txt"
	wrapper:
		"v3.3.3/bio/fastp"