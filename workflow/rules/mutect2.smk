rule mutect2_tumor: 
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

rule mutect2_control: 
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

rule orientation_model_tumor:
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

rule orientation_model_control:
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

rule get_pileup_summaries_tumor:
	input:
		bam="alignments/{sample}.tumor.dd.rec.bam",
		intervals=config["intervals"],
		variants=config["gnomAD"]
	output:
		"data/{sample}.tumor.pileupTable"
	threads: 1
	resources:
		mem_mb=4096,
	params:
		extra="",
	log:
		"logs/{sample}.GetPileupSummariesTumor.log",
	wrapper:
		"v3.3.3/bio/gatk/getpileupsummaries"

rule get_pileup_summaries_control:
	input:
		bam="alignments/{sample}.control.dd.rec.bam",
		intervals=config["intervals"],
		variants=config["gnomAD"]
	output:
		"data/{sample}.control.pileupTable"
	threads: 1
	resources:
		mem_mb=4096,
	params:
		extra="",
	log:
		"logs/{sample}.GetPileupSummariesControl.log",
	wrapper:
		"v3.3.3/bio/gatk/getpileupsummaries"

rule calculate_contamination_tumor:
	input:
		"data/{sample}.tumor.pileupTable"
	output:
		contamination="data/{sample}.tumor.contaminationTable",
		segmentation="data/{sample}.tumor.tumorSeg.txt"
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

rule calculate_contamination_control:
	input:
		"data/{sample}.control.pileupTable"
	output:
		contamination="data/{sample}.control.contaminationTable",
		segmentation="data/{sample}.control.tumorSeg.txt"
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

rule filter_mutect_calls_tumor:
	input:
		vcf="results/{sample}_tumor/{sample}.mutect2.vcf.gz",
		ref=config["genome"],
		bam="alignments/{sample}.tumor.dd.rec.bam",
		intervals=config["intervals"],
		contamination="data/{sample}.tumor.contaminationTable", # from gatk CalculateContamination
		segmentation="data/{sample}.tumor.tumorSeg.txt", # from gatk CalculateContamination
		f1r2="data/{sample}.tumor.artifacts_prior.tar.gz" # from gatk LearnReadOrientationBias
	output:
		vcf=temp("results/{sample}_tumor/{sample}.mutect2.filt.vcf.gz"),
	log:
		"logs/{sample}.FilterMutectCallsTumor.log",
	params:
		#extra="--tumor-segmentation data/{wildcard.sample}.tumseg.txt",  # optional arguments, see GATK docs
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"])  # optional
	resources:
		mem_mb=4096,
	wrapper:
		"v3.3.3/bio/gatk/filtermutectcalls"

rule filter_mutect_calls_control:
	input:
		vcf="results/{sample}_control/{sample}.mutect2.vcf.gz",
		ref=config["genome"],
		bam="alignments/{sample}.control.dd.rec.bam",
		intervals=config["intervals"],
		contamination="data/{sample}.control.contaminationTable", # from gatk CalculateContamination
		segmentation="data/{sample}.control.tumorSeg.txt", # from gatk CalculateContamination
		f1r2="data/{sample}.control.artifacts_prior.tar.gz" # from gatk LearnReadOrientationBias
	output:
		vcf=temp("results/{sample}_control/{sample}.mutect2.filt.vcf.gz")
	log:
		"logs/{sample}.FilterMutectCallsControl.log",
	params:
		#extra="--tumor-segmentation data/{wildcard.sample}.tumseg.txt",  # optional arguments, see GATK docs
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"])  # optional
	resources:
		mem_mb=4096,
	wrapper:
		"v3.3.3/bio/gatk/filtermutectcalls"


rule clean_filter_mutect_tumor_output:
	input:
		"results/{sample}_tumor/{sample}.mutect2.filt.vcf.gz"
	output:
		vcf="results/{sample}_tumor/{sample}.mutect2.filtered.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.mutect2.filtered.vcf.gz.tbi"
	threads: 1
	log:
		"logs/{sample}.clean_filter_mutect_tumor_output.log"
	params:	
		excl=config["chr_to_exclude"],
		depth=config['filtering_tumors']['min_depth'],
		vaf=config['filtering_tumors']['vaf'],
		alt=config['filtering_tumors']['alt_depth'],
		ref=config["genome"],
	conda:
		"../envs/bcftools.yaml"
	shell:	
		"""
		bcftools view -Ov -i'FILTER == "PASS" || FILTER == "germline" & POPAF > 2' {input} |
		grep -v -f {params.excl} | 
		bcftools norm -m - -f {params.ref} | 
		bcftools view -i "FORMAT/DP[0] >= {params.depth} & FORMAT/AD[0:1] >= {params.alt} & FORMAT/AF >= {params.vaf}" -Oz -o {output.vcf} &&
		tabix -p vcf {output.vcf}
		"""

rule clean_filter_mutect_control_output:
	input:
		"results/{sample}_control/{sample}.mutect2.filt.vcf.gz"
	output:
		vcf="results/{sample}_control/{sample}.mutect2.filtered.vcf.gz",
		tbi="results/{sample}_control/{sample}.mutect2.filtered.vcf.gz.tbi"
	threads: 1
	log:
		"logs/{sample}.clean_filter_mutect_control_output.log"
	params:
		excl=config["chr_to_exclude"],
		depth=config['filtering_controls']['min_depth'],
		vaf=config['filtering_controls']['vaf'],
		alt=config['filtering_controls']['alt_depth'],
		ref=config["genome"],
	conda:
		"../envs/bcftools.yaml"
	shell:	
		"""
		bcftools view -Ov -i'FILTER == "PASS" || FILTER == "germline" & POPAF > 2' {input} |
		grep -v -f {params.excl} | 
		bcftools norm -m - -f {params.ref} | 
		bcftools view -i "FORMAT/DP[0] >= {params.depth} & FORMAT/AD[0:1] >= {params.alt} & FORMAT/AF >= {params.vaf}" -Oz -o {output.vcf} &&
		tabix -p vcf {output.vcf}
		"""

rule merge_mutect2_control_vars:
	input:
		vcf=expand(f"results/{{sample}}_control/{{sample}}.mutect2.filtered.vcf.gz", sample=config["controls"].values()),
		tbi=expand(f"results/{{sample}}_control/{{sample}}.mutect2.filtered.vcf.gz.tbi", sample=config["controls"].values())
	output:
		vcf="results/custom_pon.mutect2.vcf.gz",
		tbi="results/custom_pon.mutect2.vcf.gz.tbi" 
	threads: 1
	log:
		"logs/merge_mutect2_control_vars.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} && " 
		"tabix -p vcf {output.vcf}"

rule remove_mutect2_control_vars:
	input:
		vcf="results/{sample}_tumor/{sample}.mutect2.filtered.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.mutect2.filtered.vcf.gz.tbi",
		PoN_vcf="results/custom_pon.mutect2.vcf.gz",
		PoN_tbi="results/custom_pon.mutect2.vcf.gz.tbi"
	output:
		vcf=temp("results/{sample}_tumor/{sample}.mutect2.filtOnCtr.vcf.gz"),
		tbi=temp("results/{sample}_tumor/{sample}.mutect2.filtOnCtr.vcf.gz.tbi")
	threads: 1
	log:
		"logs/{sample}.remove_mutect2_control_vars.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools isec -w1 -Oz -c none -n~10 -o {output.vcf} {input.vcf} {input.PoN_vcf} 2> {log} && "
		"tabix -p vcf {output.vcf}"


rule merge_mutect2_tumor_output:
	input:
		vcf=expand(f"results/{{sample}}_tumor/{{sample}}.mutect2.filtOnCtr.vcf.gz", sample=config["samples"].values()),
		tbi=expand(f"results/{{sample}}_tumor/{{sample}}.mutect2.filtOnCtr.vcf.gz.tbi", sample=config["samples"].values())
	output:
		vcf="results/multisample.mutect2.vcf.gz",
		tbi="results/multisample.mutect2.vcf.gz.tbi"
	threads: 1
	log:
		"logs/merge_mutect2_tumor_output.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2>{log} && "
		"tabix -p vcf {output.vcf}"

