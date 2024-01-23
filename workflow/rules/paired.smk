rule GetPileupSummariesGermline:
	input:
		bam="alignments/{sample}.germline.dd.rec.bam",
		intervals=config["intervals"],
		variants=config["gnomAD"]
	output:
		"data/{sample}.germline.pileup.table"
	threads: 1
	resources:
		mem_mb=4096,
	params:
		extra="",
	log:
		"logs/{sample}.GetPileupSummariesGermlinePaired.log",
	wrapper:
		"v3.3.3/bio/gatk/getpileupsummaries"

rule CalculateContaminationPaired:
	input:
		tumor="data/{sample}.tumor.pileup.table",
		germline="data/{sample}.germline.pileup.table"
	output:
		contamination="data/{sample}.paired.contamination.table",
		segmentation="data/{sample}.paired.tumorseg.txt"
	threads: 1
	log:
		"logs/{sample}.CalculateContaminationPaired.log",
	conda:
		"../envs/gatk4.yaml"
	params:
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"]),
		mem_mb="-Xmx4G"
	shell:
		"gatk --java-options {params.mem_mb} CalculateContamination -I {input.tumor} -matched {input.germline} -O {output.contamination} --tumor-segmentation {output.segmentation} 2>&1> {log}"

rule Mutect2Paired:
	input:
		bamT="alignments/{sample}.tumor.dd.rec.bam",
		baiT="alignments/{sample}.tumor.dd.rec.bai",
		bamC="alignments/{sample}.germline.dd.rec.bam",
		baiC="alignments/{sample}.germline.dd.rec.bai",
	output:
		vcf="results/{sample}.mutect2.paired.vcf.gz",
	threads: config["threads"]
	log:
		"logs/{sample}.mutect.paired.log"
	params:
		ref=config["genome"],
		interval=config["intervals"],
		gnomAD=config["gnomAD"],
		pon=config["gatk_pon"]
	shell:
		"""
		gatk --java-options "-Xmx4G -XX:ParallelGCThreads={threads}" Mutect2 -R {params.ref} -I {input.bamT} -I {input.bamC} -pon {params.pon} --max-reads-per-alignment-start 0 --max-mnp-distance 0 -normal {wildcards.sample} -L {params.interval} --native-pair-hmm-threads {threads} --af-of-alleles-not-in-resource 0.0000025 --germline-resource {params.gnomAD} -O {output.vcf} 2>{log}
		"""

rule FilterMutectCallsPaired:
	input:
		vcf="results/{sample}.mutect2.paired.vcf.gz",
		ref=config["genome"],
		#bam="data/bam/{sample}.dd.rec.bam",
		intervals=config["intervals"],
		contamination="data/{sample}.paired.contamination.table", # from gatk CalculateContamination
		segmentation="data/{sample}.paired.tumorseg.txt", # from gatk CalculateContamination
		#f1r2="data/{sample}.artifacts_prior.tar.gz" # from gatk LearnReadOrientationBias
	output:
		vcf="results/{sample}.mutect2.paired.filtered.vcf.gz"
	log:
		"logs/{sample}.FilterMutectCallsPaired.log",
	params:
		#extra="--tumor-segmentation data/{wildcard.sample}.tumseg.txt",  # optional arguments, see GATK docs
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"])  # optional
	resources:
		mem_mb=4096,
	wrapper:
		"v3.3.3/bio/gatk/filtermutectcalls"

#######################################################################################   VARSCAN   #######################################################################################

rule GermlineMpileup:
    input:
        # single or list of bam files
        bam="alignments/{sample}.germline.dd.rec.bam",
        reference_genome=config["genome"],
    output:
        "data/mpileup/{sample}.germline.mpileup.gz",
    log:
        "logs/{sample}.GermlineMpileup.log",
    params:
        extra="-d 10000",  # optional
    wrapper:
        "v3.3.3/bio/samtools/mpileup"

rule TumorMpileup:
    input:
        # single or list of bam files
        bam="alignments/{sample}.tumor.dd.rec.bam",
        reference_genome=config["genome"],
    output:
        "data/mpileup/{sample}.tumor.mpileup.gz",
    log:
        "logs/{sample}.TumorMpileup.log",
    params:
        extra="-d 10000",  # optional
    wrapper:
        "v3.3.3/bio/samtools/mpileup"

rule BgzipGermlineMpileup:
	input:
		"data/mpileup/{sample}.germline.mpileup.gz"
	output:
		"data/mpileup/{sample}.germline.mpileup"
	log:
		"logs/{sample}.BgzipGermlineMpileup.log"
	threads:1
	shell:
		"bgzip -d {input}"

rule BgzipTumorMpileup:
	input:
		"data/mpileup/{sample}.tumor.mpileup.gz"
	output:
		"data/mpileup/{sample}.tumor.mpileup"
	log:
		"logs/{sample}.BgzipTumorMpileup.log"
	threads:1
	shell:
		"bgzip -d {input}"

rule VarscanSomatic:
	input:
		normal_pileup = "data/mpileup/{sample}.germline.mpileup",
		tumor_pileup = "data/mpileup/{sample}.tumor.mpileup"
	output:
		snp=temp("results/{sample}.varscan.paired.snp.vcf"),
		indel=temp("results/{sample}.varscan.paired.indel.vcf")
	log:
		"logs/{sample}.VarscanSomatic.log"
	conda:
		"../envs/varscan.yaml"
	message:
		"Calling somatic variants {wildcards.sample}"
	threads:1
	shell:
		"varscan somatic {input.normal_pileup} {input.tumor_pileup} --output-vcf --min-avg-qual 15 --p-value 0.05 --min-var-freq 0.03  --output-snp {output.snp} --output-indel {output.indel} 2>{log}"

rule VarscanBgzipIndex:
	input:
		snp="results/{sample}.varscan.paired.snp.vcf",
		indel="results/{sample}.varscan.paired.indel.vcf"
	output:
		snpbg=temp("results/{sample}.varscan.paired.snp.vcf.gz"),
		indelbg=temp("results/{sample}.varscan.paired.indel.vcf.gz"),
		snpix=temp("results/{sample}.varscan.paired.snp.vcf.gz.tbi"),
		indelix=temp("results/{sample}.varscan.paired.indel.vcf.gz.tbi")
	log:
		"logs/{sample}.VarscanBgzipIndex.log"
	threads:1
	shell:
		"bgzip {input.snp} && tabix {output.snpbg} 2>{log} && bgzip {input.indel} && tabix {output.indelbg} 2>{log}"

rule MergeVarscanOutput:
    input:
        calls=["results/{sample}.varscan.paired.snp.vcf.gz", "results/{sample}.varscan.paired.indel.vcf.gz"],
    output:
        temp("results/{sample}.varscan.paired.tmp.vcf.gz")
    log:
        "logs/{sample}.MergeVarscanOutput.log",
    params:
        uncompressed_bcf=True,
        extra="-a",  # optional parameters for bcftools concat (except -o)
    threads: 1
    resources:
        mem_mb=10,
    wrapper:
        "v3.3.3/bio/bcftools/concat"

rule ReheaderVarscanOutput:
	input:
		"results/{sample}.varscan.paired.tmp.vcf.gz"
	output:
		"results/{sample}.varscan.paired.vcf.gz"
	threads: 1
	params:
		txt="data/{sample}.rh.txt"
	log:
		"logs/{sample}.varscan.rh.log"
	conda:
		"../envs/bcftools.yaml"
	shell: 
		"echo {wildcards.sample}_germline > {params.txt} && "
		"echo {wildcards.sample}_tumor >> {params.txt} && "
		"bcftools reheader -s {params.txt} -o {output} {input} && "
		"rm {params.txt}"
