rule bwaTumor:
	input:
		reads=["data/fastq/{sample}.tumor.R1.tr.fastq.gz", "data/fastq/{sample}.tumor.R2.tr.fastq.gz"],
		idx=multiext(config['genome'], ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123")
	output:
		temp("alignments/{sample}.tumor.srt.bam")
	log:
		"logs/{sample}.bwaTumor.log"
	message:
		"Mapping reads with bwa-mem2 on tumor {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.bwaTumor.txt"
	params:
		extra=r"-R '@RG\tID:{sample}_wesMasto\tSM:{sample}_tumor\tDT:20230515\tPI:150\tLB:SureSelectV8\tPL:ILLUMINA\tCN:BIODIVERSA'",
		sort="samtools",  
		sort_order="coordinate"
	threads: config["threads"]
	wrapper:
		"v3.3.3/bio/bwa-mem2/mem"

rule bwaControl:    
	input: 
		reads=["data/fastq/{sample}.control.R1.tr.fastq.gz", "data/fastq/{sample}.control.R2.tr.fastq.gz"],
		idx=multiext(config['genome'], ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123")
	output: 
		temp("alignments/{sample}.control.srt.bam")
	threads: config["threads"]
	message:
		"Mapping reads with bwa-mem2 on control {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.bwaControl.txt"
	params:
		extra=r"-R '@RG\tID:{sample}_wesMasto\tSM:{sample}_control\tDT:20230515\tPI:150\tLB:SureSelectV8\tPL:ILLUMINA\tCN:BIODIVERSA'",
		sort="samtools",  
		sort_order="coordinate"
	log:
		"logs/{sample}.bwaControl.log"
	wrapper:
		"v3.3.3/bio/bwa-mem2/mem"

rule bwaGermline:    
	input: 
		reads=["data/fastq/{sample}.germline.R1.tr.fastq.gz", "data/fastq/{sample}.germline.R2.tr.fastq.gz"],
		idx=multiext(config['genome'], ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123")
	output: 
		temp("alignments/{sample}.germline.srt.bam")
	threads: config["threads"]
	message:
		"Mapping reads with bwa-mem2 on germline {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.bwaGermline.txt"
	params:
		extra=r"-R '@RG\tID:{sample}_wesMasto\tSM:{sample}_germline\tDT:20230515\tPI:150\tLB:SureSelectV8\tPL:ILLUMINA\tCN:BIODIVERSA'",
		sort="samtools",  
		sort_order="coordinate"
	log:
		"logs/{sample}.bwaGermline.log"
	wrapper:
		"v3.3.3/bio/bwa-mem2/mem"

rule markDuplicatesTumor:
	input:
		bams="alignments/{sample}.tumor.srt.bam"
	output:
		bam=temp("alignments/{sample}.tumor.dd.bam"),
		bai=temp("alignments/{sample}.tumor.dd.bai"),
		metrics="data/{sample}.tumor.metrics.txt"
	log:
		"logs/{sample}.markDuplicatesTumor.log"
	message:
		"Mark and remove PCR-optical duplicates"
	benchmark:
		"benchmarks/{sample}.markDuplicatesTumor.txt"
	params:
		extra="--REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT",
		java_opts="-XX:ParallelGCThreads="+str(config["threads"])
	resources:
		mem_mb=4096    
	log:
		"logs/{sample}.markDuplicatesTumor.log"
	wrapper:
		 "v3.3.3/bio/picard/markduplicates"

rule markDuplicatesControl:
	input:
		bams="alignments/{sample}.control.srt.bam"
	output:
		bam=temp("alignments/{sample}.control.dd.bam"),
		bai=temp("alignments/{sample}.control.dd.bai"),
		metrics="data/{sample}.control.metrics.txt"
	log:
		"logs/{sample}.markDuplicatesControl.log"
	message:
		"Mark and remove PCR-optical duplicates"
	benchmark:
		"benchmarks/{sample}.markDuplicatesControl.txt"
	params:
		extra="--REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT",
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"])
	resources:
		mem_mb=4096    
	log:
		"logs/{sample}.markDuplicatesControl.log"
	wrapper:
		 "v3.3.3/bio/picard/markduplicates"

rule markDuplicatesGermline:
	input:
		bams="alignments/{sample}.germline.srt.bam"
	output:
		bam=temp("alignments/{sample}.germline.dd.bam"),
		bai=temp("alignments/{sample}.germline.dd.bai"),
		metrics="data/{sample}.germline.metrics.txt"
	log:
		"logs/{sample}.markDuplicatesGermline.log"
	message:
		"Mark and remove PCR-optical duplicates"
	benchmark:
		"benchmarks/{sample}.markDuplicatesGermline.txt"
	params:
		extra="--REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT",
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"])
	resources:
		mem_mb=4096    
	log:
		"logs/{sample}.markDuplicatesGermline.log"
	wrapper:
		 "v3.3.3/bio/picard/markduplicates"

rule baseRecalibratorTumor:
	input:
		bam="alignments/{sample}.tumor.dd.bam",
		bai="alignments/{sample}.tumor.dd.bai",
		ref=config['genome'],
		dict="resources/GRCh38_full_analysis_set_plus_decoy_hla.dict",
		known=["resources/dbSNP.b156.filt.vcf.gz", 
		"resources/Homo_sapiens_assembly38.known_indels.vcf.gz",
		config["gnomAD"], "resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"],  # optional known sites - single or a list
	output:
		recal_table="data/{sample}.tumor.recal.table",
	benchmark:
		"benchmarks/{sample}.baseRecalibratorTumor.txt"
	log:
		"logs/{sample}.baseRecalibratorTumor.log"
	message:
		"Perform base recalibration model"
	params:
		extra="-L " + config['intervals'],  # optional
		java_opts="-XX:ParallelGCThreads=10"  # optional
	resources:
		mem_mb=4096
	wrapper:
		"v3.3.3/bio/gatk/baserecalibrator"

rule baseRecalibratorControl:
	input:
		bam="alignments/{sample}.control.dd.bam",
		bai="alignments/{sample}.control.dd.bai",
		ref=config['genome'],
		dict="resources/GRCh38_full_analysis_set_plus_decoy_hla.dict",
		known=["resources/dbSNP.b156.filt.vcf.gz", 
		"resources/Homo_sapiens_assembly38.known_indels.vcf.gz",
		config["gnomAD"], "resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"],  # optional known sites - single or a list
	output:
		recal_table="data/{sample}.control.recal.table",
	benchmark:
		"benchmarks/{sample}.baseRecalibratorControl.txt"
	log:
		"logs/{sample}.baseRecalibratorControl.log"
	message:
		"Perform base recalibration model"
	params:
		extra="-L " + config['intervals'],  # optional
		java_opts="-XX:ParallelGCThreads=10"  # optional
	resources:
		mem_mb=4096
	wrapper:
		"v3.3.3/bio/gatk/baserecalibrator"

rule baseRecalibratorGermline:
	input:
		bam="alignments/{sample}.germline.dd.bam",
		bai="alignments/{sample}.germline.dd.bai",
		ref=config['genome'],
		dict="resources/GRCh38_full_analysis_set_plus_decoy_hla.dict",
		known=["resources/dbSNP.b156.filt.vcf.gz", 
		"resources/Homo_sapiens_assembly38.known_indels.vcf.gz",
		config["gnomAD"], "resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"],  # optional known sites - single or a list
	output:
		recal_table="data/{sample}.germline.recal.table",
	benchmark:
		"benchmarks/{sample}.baseRecalibratorGermline.txt"
	log:
		"logs/{sample}.baseRecalibratorGermline.log"
	message:
		"Perform base recalibration model"
	params:
		extra="-L " + config['intervals'],  # optional
		java_opts="-XX:ParallelGCThreads=10"  # optional
	resources:
		mem_mb=4096
	wrapper:
		"v3.3.3/bio/gatk/baserecalibrator"

rule ApplyBQSRTumor:
	input:
		bam="alignments/{sample}.tumor.dd.bam",
		bai="alignments/{sample}.tumor.dd.bai",	
		ref=config['genome'],
		dict="resources/GRCh38_full_analysis_set_plus_decoy_hla.dict",
		recal_table="data/{sample}.tumor.recal.table",
	output:
		bam="alignments/{sample}.tumor.dd.rec.bam",
		bai="alignments/{sample}.tumor.dd.rec.bai"
	benchmark:
		"benchmarks/{sample}.ApplyBQSRTumor.txt"
	log:
		"logs/recal/{sample}.ApplyBQSRTumor.log",
	params:
		extra="-L " + config['intervals'],  # optional
		java_opts="-XX:ParallelGCThreads=10",  # optional
		#embed_ref=True,  # embed the reference in cram output
	resources:
		mem_mb=4096
	wrapper:
		"v3.3.3/bio/gatk/applybqsr"

rule ApplyBQSRControl:
	input:
		bam="alignments/{sample}.control.dd.bam",
		bai="alignments/{sample}.control.dd.bai",   
		ref=config['genome'],
		dict="resources/GRCh38_full_analysis_set_plus_decoy_hla.dict",
		recal_table="data/{sample}.control.recal.table",
	output:
		bam="alignments/{sample}.control.dd.rec.bam",
		bai="alignments/{sample}.control.dd.rec.bai"
	benchmark:
		"benchmarks/{sample}.ApplyBQSRControl.txt"
	log:
		"logs/recal/{sample}.ApplyBQSRControl.log",
	params:
		extra="-L " + config['intervals'],  # optional
		java_opts="-XX:ParallelGCThreads=10",  # optional
		#embed_ref=True,  # embed the reference in cram output
	resources:
		mem_mb=4096
	wrapper:
		"v3.3.3/bio/gatk/applybqsr"

rule ApplyBQSRGermline:
	input:
		bam="alignments/{sample}.germline.dd.bam",
		bai="alignments/{sample}.germline.dd.bai",   
		ref=config['genome'],
		dict="resources/GRCh38_full_analysis_set_plus_decoy_hla.dict",
		recal_table="data/{sample}.germline.recal.table",
	output:
		bam="alignments/{sample}.germline.dd.rec.bam",
		bai="alignments/{sample}.germline.dd.rec.bai"
	benchmark:
		"benchmarks/{sample}.ApplyBQSRGermline.txt"
	log:
		"logs/recal/{sample}.ApplyBQSRGermline.log",
	params:
		extra="-L " + config['intervals'],  # optional
		java_opts="-XX:ParallelGCThreads=10",  # optional
		#embed_ref=True,  # embed the reference in cram output
	resources:
		mem_mb=4096
	wrapper:
		"v3.3.3/bio/gatk/applybqsr"
