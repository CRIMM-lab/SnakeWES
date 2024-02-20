########################################################################### Mutect2 ###########################################################################
rule Mutect2TumorPoN: 
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
		"logs/{sample}.Mutect2TumorPoN.log",
	params:
		extra="--max-reads-per-alignment-start 0" 
	wrapper:
		"v3.3.3/bio/gatk/mutect"

rule Mutect2ControlPoN: 
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
		"logs/{sample}.Mutect2ControlPoN.log",
	params:
		extra="--max-reads-per-alignment-start 0" 
	wrapper:
		"v3.3.3/bio/gatk/mutect"

rule OrientationModelTumorPoN:
	input:
		f1r2="data/{sample}.tumor.f1r2.tar.gz",
	output:
		"data/{sample}.tumor.artifacts_prior.tar.gz",
	resources:
		mem_mb=4096,
	log:
		"logs/{sample}.OrientationModelTumorPoN.log",
	wrapper:
		"v3.3.3/bio/gatk/learnreadorientationmodel"

rule OrientationModelControlPoN:
	input:
		f1r2="data/{sample}.control.f1r2.tar.gz",
	output:
		"data/{sample}.control.artifacts_prior.tar.gz",
	resources:
		mem_mb=4096,
	log:
		"logs/{sample}.OrientationModelControlPoN.log",
	wrapper:
		"v3.3.3/bio/gatk/learnreadorientationmodel"

rule GetPileupSummariesTumorPoN:
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
		"logs/{sample}.GetPileupSummariesTumorPoN.log",
	wrapper:
		"v3.3.3/bio/gatk/getpileupsummaries"

rule GetPileupSummariesControlPoN:
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
		"logs/{sample}.GetPileupSummariesControlPoN.log"
	wrapper:
		"v3.3.3/bio/gatk/getpileupsummaries"

rule CalculateContaminationTumorPoN:
	input:
		"data/{sample}.tumor.pileupTable"
	output:
		contamination="data/{sample}.tumor.contaminationTable",
		segmentation="data/{sample}.tumor.tumorSeg.txt"
	threads: 1
	log:
		"logs/{sample}.CalculateContaminationTumorPoN.log"
	conda:
		"../envs/gatk4.yaml"
	params:
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"]),
		mem_mb="-Xmx4G"
		#extra="--tumor-segmentation data/{wildcards.sample}.tumseg.txt"
	shell:
		"gatk --java-options {params.mem_mb} CalculateContamination -I {input} -O {output.contamination} --tumor-segmentation {output.segmentation} > {log} 2>&1"

rule CalculateContaminationControlPoN:
	input:
		"data/{sample}.control.pileupTable"
	output:
		contamination="data/{sample}.control.contaminationTable",
		segmentation="data/{sample}.control.tumorSeg.txt"
	threads: 1
	log:
		"logs/{sample}.CalculateContaminationControlPoN.log"
	conda:
		"../envs/gatk4.yaml"
	params:
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"]),
		mem_mb="-Xmx4G"
		#extra="--tumor-segmentation data/{wildcards.sample}.tumseg.txt"
	shell:
		"gatk --java-options {params.mem_mb} CalculateContamination -I {input} -O {output.contamination} --tumor-segmentation {output.segmentation} > {log} 2>&1"

rule FilterMutectCallsTumorPoN:
	input:
		vcf="results/{sample}_tumor/{sample}.mutect2.vcf.gz",
		ref=config["genome"],
		bam="alignments/{sample}.tumor.dd.rec.bam",
		intervals=config["intervals"],
		contamination="data/{sample}.tumor.contaminationTable", # from gatk CalculateContamination
		segmentation="data/{sample}.tumor.tumorSeg.txt", # from gatk CalculateContamination
		f1r2="data/{sample}.tumor.artifacts_prior.tar.gz" # from gatk LearnReadOrientationBias
	output:
		vcf=temp("results/{sample}_tumor/{sample}.mutect2.filt.vcf.gz")
	log:
		"logs/{sample}.FilterMutectCallsTumorPoN.log"
	params:
		#extra="--tumor-segmentation data/{wildcard.sample}.tumseg.txt",  # optional arguments, see GATK docs
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"])  # optional
	resources:
		mem_mb=4096,
	wrapper:
		"v3.3.3/bio/gatk/filtermutectcalls"

rule FilterMutectCallsControlPoN:
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
		"logs/{sample}.FilterMutectCallsControlPoN.log"
	params:
		#extra="--tumor-segmentation data/{wildcard.sample}.tumseg.txt",  # optional arguments, see GATK docs
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"])  # optional
	resources:
		mem_mb=4096,
	wrapper:
		"v3.3.3/bio/gatk/filtermutectcalls"

rule CleanFilterMutectTumorOutputPoN:
	input:
		"results/{sample}_tumor/{sample}.mutect2.filt.vcf.gz"
	output:
		vcf="results/{sample}_tumor/{sample}.mutect2.filtered.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.mutect2.filtered.vcf.gz.tbi"
	threads: 1
	log:
		"logs/{sample}.CleanFilterMutectTumorOutputPoN.log"
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

rule CleanFilterMutectControlOutputPoN:
	input:
		"results/{sample}_control/{sample}.mutect2.filt.vcf.gz"
	output:
		vcf="results/{sample}_control/{sample}.mutect2.filtered.vcf.gz",
		tbi="results/{sample}_control/{sample}.mutect2.filtered.vcf.gz.tbi"
	threads: 1
	log:
		"logs/{sample}.CleanFilterMutectControlOutputPoN.log"
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

rule MergeMutect2ControlVarsPoN:
	input:
		vcf=expand(f"results/{{sample}}_control/{{sample}}.mutect2.filtered.vcf.gz", sample=config["controls"].values()),
		tbi=expand(f"results/{{sample}}_control/{{sample}}.mutect2.filtered.vcf.gz.tbi", sample=config["controls"].values())
	output:
		vcf="results/custom_pon.mutect2.vcf.gz",
		tbi="results/custom_pon.mutect2.vcf.gz.tbi" 
	threads: 1
	log:
		"logs/MergeMutect2ControlVarsPoN.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} && " 
		"tabix -p vcf {output.vcf}"

rule RemoveMutect2ControlVarsPoN:
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
		"logs/{sample}.RemoveMutect2ControlVarsPoN.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools isec -w1 -Oz -c none -n~10 -o {output.vcf} {input.vcf} {input.PoN_vcf} 2> {log} && "
		"tabix -p vcf {output.vcf}"


rule MergeMutect2TumorOutputPoN:
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

########################################################################### Bcftools ###########################################################################

rule BcftoolsCallControlsPoN:  
	input:
		bam="alignments/{sample}.control.dd.rec.bam",
		bai="alignments/{sample}.control.dd.rec.bai"
	output:
		"results/{sample}_control/{sample}.bcftools.vcf.gz"
	threads: 1
	params:
		ref=config["genome"],
		intervals=config["intervals"]
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.BcftoolsCallControlsPoN.log"
	shell:
		"bcftools mpileup -Ou -d 10000 -R {params.intervals} -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR -f {params.ref} {input.bam} | "
		"bcftools call -mv | "
		"bcftools norm -m - -f {params.ref} -Oz -o {output} 2> {log}"


rule BcftoolsCallTumorsPoN:  
	input:
		bam="alignments/{sample}.tumor.dd.rec.bam",
		bai="alignments/{sample}.tumor.dd.rec.bai"
	output:
		"results/{sample}_tumor/{sample}.bcftools.vcf.gz"
	threads: 1
	params:
		ref=config["genome"],
		intervals=config["intervals"]
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.BcftoolsCallTumorsPoN.log"
	shell:
		"bcftools mpileup -Ou -d 10000 -R {params.intervals} -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR -f {params.ref} {input.bam} | "
		"bcftools call -mv | "
		"bcftools norm -m - -f {params.ref} -Oz -o {output} 2> {log}"

rule BcftoolsAddAFFieldControlsPoN:  
	input:
		"results/{sample}_control/{sample}.bcftools.vcf.gz"
	output:
		vcf=temp("results/{sample}_control/{sample}.bcftools.addaf.vcf.gz"),
		tsv=temp("results/{sample}_control/{sample}.annot.tsv.gz"),
		tsv_tbi=temp("results/{sample}_control/{sample}.annot.tsv.gz.tbi")
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.BcftoolsAddAFFieldControlsPoN.log"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%AD{{0}}\t%AD{{1}}]' {input} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {output.tsv} && 
		tabix -b2 -e2 {output.tsv} && 
		bcftools annotate -a {output.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF {input} -Oz -o {output.vcf}
		"""

rule BcftoolsAddAFFieldTumorsPoN:  
	input:
		"results/{sample}_tumor/{sample}.bcftools.vcf.gz"
	output:
		vcf=temp("results/{sample}_tumor/{sample}.bcftools.addaf.vcf.gz"),
		tsv=temp("results/{sample}_tumor/{sample}.annot.tsv.gz"),
		tsv_tbi=temp("results/{sample}_tumor/{sample}.annot.tsv.gz.tbi")
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.BcftoolsAddAFFieldTumorsPoN.log"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%AD{{0}}\t%AD{{1}}]' {input} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {output.tsv} && 
		tabix -b2 -e2 {output.tsv} && 
		bcftools annotate -a {output.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF {input} -Oz -o {output.vcf}
		"""

rule BcftoolsFilterControlsPoN: 
	input:
		"results/{sample}_control/{sample}.bcftools.addaf.vcf.gz"
	output:
		vcf=temp("results/{sample}_control/{sample}.bcftools.filt.vcf.gz"),
		tbi=temp("results/{sample}_control/{sample}.bcftools.filt.vcf.gz.tbi")
	threads: 1
	params:
		excl=config["chr_to_exclude"],
		depth=config['filtering_controls']['min_depth'],
		vaf=config['filtering_controls']['vaf'],
		alt=config['filtering_controls']['alt_depth']		
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.BcftoolsFilterControlsPoN.log"
	shell:
		"bcftools view -i'FORMAT/DP >= {params.depth} & FORMAT/AF >= {params.vaf} & FORMAT/AD[0:1] >= {params.alt}' {input} | "
		"grep -v -f {params.excl} | "
		"bcftools sort -Oz -o {output.vcf} 2> {log} && "
		"tabix -p vcf {output.vcf}"

rule BcftoolsFilterTumorsPoN: 
	input:
		"results/{sample}_tumor/{sample}.bcftools.addaf.vcf.gz"
	output:
		vcf=temp("results/{sample}_tumor/{sample}.bcftools.filt.vcf.gz"),
		tbi=temp("results/{sample}_tumor/{sample}.bcftools.filt.vcf.gz.tbi")
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
		"logs/{sample}.BcftoolsFilterTumorsPoN.log"
	shell:
		"bcftools view -i 'FORMAT/DP >= {params.depth} & FORMAT/AF >= {params.vaf} & FORMAT/AD[0:1] >= {params.alt}' {input} | "
		"grep -v -f {params.excl} | "
		"bcftools sort -Oz -o {output.vcf} 2> {log} && "
		"tabix -p vcf {output.vcf}"

rule BcftoolsMergeBcftoolsControlsPoN:
	input:
		vcf=expand(f"results/{{sample}}_control/{{sample}}.bcftools.filt.vcf.gz", sample=config['controls'].values()),
		tbi=expand(f"results/{{sample}}_control/{{sample}}.bcftools.filt.vcf.gz.tbi", sample=config['controls'].values())
	output:
		vcf="results/custom_pon.bcftools.vcf.gz",
		tbi="results/custom_pon.bcftools.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/BcftoolsMergeBcftoolsControlsPoN.log"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2> {log} && "
		"tabix -p vcf {output.vcf}"


rule BcftoolsRemoveControlsPoN: 
	input:
		vcf="results/{sample}_tumor/{sample}.bcftools.filt.vcf.gz",
		vcf_tbi="results/{sample}_tumor/{sample}.bcftools.filt.vcf.gz.tbi",
		custom_pon="results/custom_pon.bcftools.vcf.gz",
		pon_tbi="results/custom_pon.bcftools.vcf.gz.tbi"
	output:
		vcf="results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/{sample}.BcftoolsRemoveControlsPoN.log"
	shell:
		"bcftools isec -w1 -Oz -c none -n~10 -o {output.vcf} {input.vcf} {input.custom_pon} 2> {log} && "
		"tabix -p vcf {output.vcf}"	


rule MergeBcftoolsTumorVariantsPoN: 
	input:
		vcf=expand("results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vcf.gz", sample=config['samples'].values()),
		tbi=expand("results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vcf.gz.tbi", sample=config['samples'].values())
	output:
		vcf="results/multisample.bcftools.vcf.gz",
		tbi="results/multisample.bcftools.vcf.gz.tbi"
	threads: 1
	log: 
		"logs/multisample.bcftools.log"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2> {log} && "
		"tabix -p vcf {output.vcf}"

########################################################################### Varscan ###########################################################################

rule VarscanCallSnvControlsPoN: 
	input:
		bam="alignments/{sample}.control.dd.rec.bam",
		bai="alignments/{sample}.control.dd.rec.bai"
	output:
		vcf_snv=temp("results/{sample}_control/{sample}.varscan.snv.vcf.gz"),
		tbi_snv=temp("results/{sample}_control/{sample}.varscan.snv.vcf.gz.tbi")
	message:
		"varscan call snv - control {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.varscan.snv.callandNorm.txt"
	threads: 1
	params:
		intervals=config["intervals"],
		ref=config["genome"],
	conda:
		"../envs/varscan.yaml"
	log:
		"logs/{sample}.VarscanCallSnvControlsPoN.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} | "
		"varscan mpileup2snp --output-vcf 1 --min-avg-qual 1 --p-value 1 --min-var-freq 0 --min-coverage 1 --min-reads2 1 | "
		"bgzip -c > {output.vcf_snv} 2> {log} &&"
		"tabix -p vcf {output.vcf_snv}"

rule VarscanCallSnvTumorsPoN: 
	input:
		bam="alignments/{sample}.tumor.dd.rec.bam",
		bai="alignments/{sample}.tumor.dd.rec.bai"
	output:
		vcf_snv=temp("results/{sample}_tumor/{sample}.varscan.snv.vcf.gz"),
		tbi_snv=temp("results/{sample}_tumor/{sample}.varscan.snv.vcf.gz.tbi")
	threads: 1
	params:
		#varscan="/home/simone/programs/VarScan.v2.3.9.jar",
		intervals=config["intervals"],
		ref=config["genome"]
	conda:
		"../envs/varscan.yaml"
	log:
		"logs/{sample}.VarscanCallSnvTumorsPoN.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} |"
		"varscan mpileup2snp --output-vcf --min-avg-qual 1 --p-value 1 --min-var-freq 0 --min-coverage 1 --min-reads2 1 |"
		"bgzip -c > {output.vcf_snv} &&"
		"tabix -p vcf {output.vcf_snv}"


rule VarscanCallIndelControlsPoN: 
	input:
		bam="alignments/{sample}.control.dd.rec.bam",
		bai="alignments/{sample}.control.dd.rec.bai"
	output:
		vcf_indel=temp("results/{sample}_control/{sample}.varscan.indel.vcf.gz"),
		tbi_indel=temp("results/{sample}_control/{sample}.varscan.indel.vcf.gz.tbi")
	message:
		"varscan call indel - control {wildcards.sample}"
	benchmark:
		"benchmarks/{sample}.varscan.indel.callandNorm.txt"
	threads: 1
	params:
		intervals=config["intervals"],
		ref=config["genome"]
	conda:
		"../envs/varscan.yaml"
	log:
		"logs/{sample}.varscan.indel.callandNorm.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} |"
		"varscan mpileup2indel --output-vcf --min-avg-qual 1 --p-value 1 --min-var-freq 0 --min-coverage 1 --min-reads2 1 |"
		"bgzip -c > {output.vcf_indel} 2> {log} &&"
		"tabix -p vcf {output.vcf_indel}"
		
rule VarscanCallIndelTumorsPoN: 
	input:
		bam="alignments/{sample}.tumor.dd.rec.bam",
		bai="alignments/{sample}.tumor.dd.rec.bai"
	output:
		vcf_indel=temp("results/{sample}_tumor/{sample}.varscan.indel.vcf.gz"),
		tbi_indel=temp("results/{sample}_tumor/{sample}.varscan.indel.vcf.gz.tbi")
	threads: 1
	params:
		intervals=config["intervals"],
		ref=config["genome"]
	conda:
		"../envs/varscan.yaml"
	log:
		"logs/{sample}.VarscanCallIndelTumorsPoN.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} | "
		"varscan mpileup2indel --output-vcf --min-avg-qual 1 --p-value 1 --min-var-freq 0 --min-coverage 1 --min-reads2 1 | "
		"bgzip -c > {output.vcf_indel} && "
		"tabix -p vcf {output.vcf_indel} "

rule VarscanConcatControlsPoN: 
	input:
		vcf_snv="results/{sample}_control/{sample}.varscan.snv.vcf.gz",
		vcf_indel="results/{sample}_control/{sample}.varscan.indel.vcf.gz",
		tbi_snv="results/{sample}_control/{sample}.varscan.snv.vcf.gz.tbi",
		tbi_indel="results/{sample}_control/{sample}.varscan.indel.vcf.gz.tbi"
	output:
		temp("results/{sample}_control/{sample}.varscan.concat.tmp.vcf.gz")
	params: 
		config["genome"]
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.VarscanConcatControlsPoN.log"
	shell: 
		"bcftools concat -a -Ov {input.vcf_snv} {input.vcf_indel} | "
		"bcftools norm -Oz -m - -f {params} -Oz -o {output} 2> {log}"				

rule VarscanConcatTumorsPoN: 
	input:
		vcf_snv="results/{sample}_tumor/{sample}.varscan.snv.vcf.gz",
		vcf_indel="results/{sample}_tumor/{sample}.varscan.indel.vcf.gz",
		tbi_snv="results/{sample}_tumor/{sample}.varscan.snv.vcf.gz.tbi",
		tbi_indel="results/{sample}_tumor/{sample}.varscan.indel.vcf.gz.tbi"
	output:
		temp("results/{sample}_tumor/{sample}.varscan.concat.tmp.vcf.gz")
	params: 
		config["genome"]
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.VarscanConcatTumorsPoN.log"
	shell: 
		"bcftools concat -a -Ov {input.vcf_snv} {input.vcf_indel} | "
		"bcftools norm -Oz -m - -f {params} -Oz -o {output} 2> {log}"

rule VarscanModAFControlsPoN:
	input:
		"results/{sample}_control/{sample}.varscan.concat.tmp.vcf.gz"
	output:
		vcf="results/{sample}_control/{sample}.varscan.concat.vcf.gz",
		tsv="results/{sample}_control/{sample}.varscan.annot.tsv.gz",
		tsv_tbi=temp("results/{sample}_control/{sample}.varscan.annot.tsv.gz.tbi")
	log:
		"logs/{sample}.VarscanModAFControlsPoN.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%FREQ]' {input} | 
		awk -F'\t' 'BEGIN {{OFS="\t"}} {{ $5 = $5 / 100; print }}'| bgzip -c > {output.tsv} && 
		tabix -b2 -e2 {output.tsv} &&
		bcftools annotate -x FORMAT/FREQ -a {output.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF  {input} -Oz -o {output.vcf}
		"""

rule VarscanModAFTumorsPoN:
	input:
		"results/{sample}_tumor/{sample}.varscan.concat.tmp.vcf.gz"
	output:
		vcf="results/{sample}_tumor/{sample}.varscan.concat.vcf.gz",
		tsv="results/{sample}_tumor/{sample}.varscan.annot.tsv.gz",
		tsv_tbi="results/{sample}_tumor/{sample}.varscan.annot.tsv.gz.tbi"
	log:
		"logs/{sample}.VarscanModAFTumorsPoN.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%FREQ]' {input} | awk -F'\t' 'BEGIN {{OFS="\t"}} {{ $5 = $5 / 100; print }}'| bgzip -c > {output.tsv} && 
		tabix -b2 -e2 {output.tsv} &&
		bcftools annotate -x FORMAT/FREQ -a {output.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF  {input} -Oz -o {output.vcf}
		"""

rule VarscanReheaderControlsPoN:
	input:
		"results/{sample}_control/{sample}.varscan.concat.vcf.gz"
	output:
		"results/{sample}_control/{sample}.varscan.vcf.gz"
	threads: 1
	params:
		txt="results/{sample}_control/{sample}.rh.txt"
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.VarscanReheaderControlsPoN.log"
	shell: 
		"echo {wildcards.sample} > {params.txt} && "
		"bcftools reheader -s {params.txt} -o {output} {input} && "
		"rm {params.txt}"


rule VarscanReheaderTumorsPoN:
	input:
		"results/{sample}_tumor/{sample}.varscan.concat.vcf.gz"
	output:
		"results/{sample}_tumor/{sample}.varscan.vcf.gz"
	threads: 1
	params:
		txt="results/{sample}_tumor/{sample}.rh.txt"
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.VarscanReheaderTumorsPoN.log"
	shell: 
		"echo {wildcards.sample} > {params.txt} && "
		"bcftools reheader -s {params.txt} -o {output} {input} && "
		"rm {params.txt}"

rule FilterVarscanControlsVars: 
	input:
		"results/{sample}_control/{sample}.varscan.vcf.gz"
	output:
		vcf="results/{sample}_control/{sample}.varscan.filt.vcf.gz",
		tbi="results/{sample}_control/{sample}.varscan.filt.vcf.gz.tbi"
	threads: 1
	params:
		excl=config["chr_to_exclude"],
		alt=config["filtering_controls"]["alt_depth"],
		min_depth=config["filtering_controls"]["min_depth"],
		vaf=config["filtering_controls"]["vaf"],
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.FilterVarscanControlsVars.log"
	shell:
		"bcftools view -Ov -i 'FORMAT/DP >= {params.min_depth} & FORMAT/AD >= {params.alt} & FORMAT/AF >= {params.vaf}' {input} |"
		"grep -v -f {params.excl} |"
		"bcftools sort -Oz -o {output.vcf} 2> {log} &&"
		"tabix -f -p vcf {output.vcf}"

rule FilterVarscanTumorsVarsPoN: 
	input:
		"results/{sample}_tumor/{sample}.varscan.vcf.gz"
	output:
		vcf="results/{sample}_tumor/{sample}.varscan.filt.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.varscan.filt.vcf.gz.tbi"
	threads: 1
	params:
		excl=config["chr_to_exclude"],
		alt=config["filtering_tumors"]["alt_depth"],
		min_depth=config["filtering_tumors"]["min_depth"],
		vaf=config["filtering_tumors"]["vaf"],
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.FilterVarscanTumorsVarsPoN.log"
	shell:
		"bcftools view -Ov -i'FORMAT/DP >= {params.min_depth} & FORMAT/AD >= {params.alt} & FORMAT/AF >= {params.vaf}' {input} | "
		"grep -v -f {params.excl} | "
		"bcftools sort -Oz -o {output.vcf} && "
		"tabix -f -p vcf {output.vcf}"

rule MergeVarscanControlsVariantsPoN:
	input:
		vcf=expand(f"results/{{sample}}_control/{{sample}}.varscan.filt.vcf.gz", sample=config['controls'].values()),
		tbi=expand(f"results/{{sample}}_control/{{sample}}.varscan.filt.vcf.gz.tbi", sample=config['controls'].values())
	output:
		vcf="results/custom_pon.varscan.vcf.gz",
		tbi="results/custom_pon.varscan.vcf.gz.tbi"
	message:
		"merge varscan controls vars"
	benchmark:
		"benchmarks/varscan.merge.ctr.vars.txt"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/MergeVarscanControlsVariantsPoN.log"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2> {log} &&"
		"tabix -p vcf {output.vcf} "


rule RemoveVarscanControlsVarsPoN: 
	input:
		vcf="results/{sample}_tumor/{sample}.varscan.filt.vcf.gz",
		custom_pon="results/custom_pon.varscan.vcf.gz",
		tbi="results/custom_pon.varscan.vcf.gz.tbi"
	output:
		vcf="results/{sample}_tumor/{sample}.varscan.filtOnCtr.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.varscan.filtOnCtr.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/{sample}.RemoveVarscanControlsVarsPoN.log"
	shell:
		"bcftools isec -w1 -Oz -c none -n~10 -o {output.vcf} {input.vcf} {input.custom_pon} 2> {log} && "
		"tabix -p vcf {output.vcf}"


rule MergeVarscanTumorsVariantsPoN: ## difference between single curly and doucle curly brackets ##
	input:
		vcf=expand(f"results/{{sample}}_tumor/{{sample}}.varscan.filtOnCtr.vcf.gz", sample=config['samples'].values()),
		tbi=expand(f"results/{{sample}}_tumor/{{sample}}.varscan.filtOnCtr.vcf.gz.tbi", sample=config['samples'].values())
	output:
		vcf="results/multisample.varscan.vcf.gz",
		tbi="results/multisample.varscan.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/multisample.MergeVarscanTumorsVariantsPoN.log"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2> {log} &&"
		"tabix -p vcf {output.vcf}"

########################################################################### Freebayes ###########################################################################

rule FreebayesCallNormTumorPoN: 
	input:
		bam="alignments/{sample}.tumor.dd.rec.bam",
		bai="alignments/{sample}.tumor.dd.rec.bai"
	output:
		"results/{sample}_tumor/{sample}.freebayes.vcf.gz"
	threads: 1
	conda:
		"../envs/freebayes.yaml"
	params:
		ref=config["genome"],
		intervals=config["intervals"]
	log:
		"logs/{sample}.FreebayesCallNormTumorPoN.log"
	shell:
		"freebayes -f {params.ref} -F 0.01 -C 2 -t {params.intervals} --pooled-continuous {input.bam} | "
		"bcftools norm -m - -f {params.ref} -Oz -o {output} 2>{log}"


rule FreebayesCallNormControlPoN: 
	input:
		bam="alignments/{sample}.control.dd.rec.bam",
		bai="alignments/{sample}.control.dd.rec.bai"
	output:
		"results/{sample}_control/{sample}.freebayes.vcf.gz"
	threads: 1
	conda:
		"../envs/freebayes.yaml"
	params:
		ref=config["genome"],
		intervals=config["intervals"]
	log:
		"logs/{sample}.FreebayesCallNormControlPoN.log"
	shell:
		"freebayes -f {params.ref} -F 0.01 -C 2 -t {params.intervals} --pooled-continuous {input.bam} | "
		"bcftools norm -m - -f {params.ref} -Oz -o {output} 2>{log}"

rule FreebayesAddAFTumorPoN:
	input:
		"results/{sample}_tumor/{sample}.freebayes.vcf.gz"
	output:
		vcf="results/{sample}_tumor/{sample}.freebayes.annotAF.vcf.gz",
		tsv="results/{sample}_tumor/{sample}.freebayes.annot.tsv.gz",
		tsv_tbi="results/{sample}_tumor/{sample}.freebayes.annot.tsv.gz.tbi"
	log:
		"logs/{sample}.FreebayesAddAFTumorPoN.log",
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%RO\t%AO\n' {input} 2> {log} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {output.tsv} 2>> {log} && 
		tabix -b2 -e2 {output.tsv} && 
		bcftools annotate -x INFO/AF -a {output.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF -Oz -o {output.vcf} {input} 2>> {log}  
		"""

rule FreebayesAddAFControlsPoN:
	input:
		"results/{sample}_control/{sample}.freebayes.vcf.gz"
	output:
		vcf="results/{sample}_control/{sample}.freebayes.annotAF.vcf.gz",
		tsv="results/{sample}_control/{sample}.freebayes.annot.tsv.gz",
		tsv_tbi="results/{sample}_control/{sample}.freebayes.annot.tsv.gz.tbi"
	log:
		"logs/{sample}.FreebayesAddAFControlsPoN.log",
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%RO\t%AO\n' {input} 2> {log} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {output.tsv} 2>> {log} && 
		tabix -b2 -e2 {output.tsv} && 
		bcftools annotate -x INFO/AF -a {output.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF -Oz -o {output.vcf} {input} 2>> {log}  
		"""

rule FreebayesFilterTumorsPoN: 
	input:
		"results/{sample}_tumor/{sample}.freebayes.annotAF.vcf.gz"
	output:
		"results/{sample}_tumor/{sample}.freebayes.filt.vcf.gz"
	log:
		"logs/{sample}.FreebayesFilterTumorsPoN.log",
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	params:
		excl=config["chr_to_exclude"],
		depth=config['filtering_tumors']['min_depth'],
		vaf=config['filtering_tumors']['vaf'],
		alt=config['filtering_tumors']['alt_depth']
	shell:
		"bcftools view -i 'QUAL > 1 & INFO/DP >= {params.depth} & FORMAT/AF >= {params.vaf} & INFO/AC >= {params.alt}' {input} 2> {log}| "
		"grep -v -f {params.excl} 2>> {log} | "
		"bcftools sort -Oz -o {output} 2>> {log} && "
		"tabix -p vcf {output} 2>> {log}"

rule FreebayesFilterControlsPoN: 
	input:
		"results/{sample}_control/{sample}.freebayes.annotAF.vcf.gz"
	output:
		vcf="results/{sample}_control/{sample}.freebayes.filt.vcf.gz",
		tbi="results/{sample}_control/{sample}.freebayes.filt.vcf.gz.tbi"
	log:
		"logs/{sample}.FreebayesFilterControlsPoN.log",
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	params:
		excl=config["chr_to_exclude"],
		depth=config['filtering_controls']['min_depth'],
		vaf=config['filtering_controls']['vaf'],
		alt=config['filtering_controls']['alt_depth']
	shell:
		"bcftools view -i 'QUAL > 1 & INFO/DP >= {params.depth} & FORMAT/AF >= {params.vaf} & INFO/AC >= {params.alt}' {input} 2> {log}| "
		"grep -v -f {params.excl} 2>> {log}| "
		"bcftools sort -Oz -o {output.vcf} 2>> {log} && "
		"tabix -p vcf {output.vcf}"

rule FreebayesMergeControlVarsPoN:
	input:
		vcf=expand(f"results/{{sample}}_control/{{sample}}.freebayes.filt.vcf.gz", sample=config['controls'].values()),
		tbi=expand(f"results/{{sample}}_control/{{sample}}.freebayes.filt.vcf.gz.tbi", sample=config['controls'].values())
	output:
		vcf="results/customPoN.freebayes.vcf.gz",
		tbi="results/customPoN.freebayes.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/FreebayesMergeControlVarsPoN.log"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2> {log} && "
		"tabix -p vcf {output.vcf} 2>> {log}"

rule FreebayesRemoveControlVarsPoN: 
	input:
		vcf="results/{sample}_tumor/{sample}.freebayes.filt.vcf.gz",
		PoN="results/customPoN.freebayes.vcf.gz",
		tbi="results/customPoN.freebayes.vcf.gz.tbi"
	output:
		vcf="results/{sample}_tumor/{sample}.freebayes.filtOnCtr.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.freebayes.filtOnCtr.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/{sample}.FreebayesRemoveControlVarsPoN.log"
	shell:
		"bcftools isec -w1 -Oz -c none -n~10 -o {output.vcf} {input.vcf} {input.PoN} 2> {log} && "
		"tabix -p vcf {output.vcf}"

rule MergeFreebayesTumorOutputPoN:
	input:
		vcf=expand(f"results/{{sample}}_tumor/{{sample}}.freebayes.filtOnCtr.vcf.gz", sample=config["samples"].values()),
		tbi=expand(f"results/{{sample}}_tumor/{{sample}}.freebayes.filtOnCtr.vcf.gz.tbi", sample=config["samples"].values())
	output:
		vcf="results/multisample.freebayes.vcf.gz",
		tbi="results/multisample.freebayes.vcf.gz.tbi"
	threads: 1
	log:
		"logs/multisample.MergeFreebayesTumorOutputPoN.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2>{log} && "
		"tabix -p vcf {output.vcf} 2>>{log}"

########################################################################### VEP ###########################################################################

rule VepSingleSamplePoN:
	input:
		vcf="results/{sample}_tumor/{sample}.{caller}.filtOnCtr.vcf.gz",
		cache=config["vepcache"]
	output:
		vcf="results/{sample}_tumor/{sample}.{caller}.filtOnCtr.vep.vcf.gz",  
		tbi="results/{sample}_tumor/{sample}.{caller}.filtOnCtr.vep.vcf.gz.tbi",
		html="results/{sample}_tumor/{sample}.{caller}.filtOnCtr.vep.html"
	threads: 
		config["threads"]
	log:
		"logs/{sample}.Vep{caller}SingleSamplePoN.log"
	conda:
		"../envs/vep.yaml"
	params:
		ref=config["genome"],
		clinvar=config["clinvar"],
		dbNSFP=config["dbNSFP"],
	shell:
		"vep -i {input.vcf} -o {output.vcf} --fork {threads} --compress_output bgzip --everything --offline --species homo_sapiens --stats_file {output.html} --assembly GRCh38 --cache --dir_cache {input.cache} --cache_version 110 --merged --fasta {params.ref} --format vcf --symbol --no_intergenic --merged --cache --pick --pick_allele --vcf --plugin dbNSFP,{params.dbNSFP},SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,FATHMM_converted_rankscore,FATHMM_pred --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN 2>{log} && tabix {output.vcf} 2>>{log}"

rule VepMultisamplePoN:
	input:
		vcf="results/multisample.{caller}.vcf.gz",
		cache=config["vepcache"]
	output:
		vcf="results/multisample.{caller}.vep.vcf.gz",
		tbi="results/multisample.{caller}.vep.vcf.gz.tbi",
		html="results/multisample.{caller}.vep.html"
	threads:
		config["threads"]
	log:
		"logs/Vep{caller}MultisamplePoN.log"
	conda:
		"../envs/vep.yaml"
	params:
		ref=config["genome"],
		clinvar=config["clinvar"],
		dbNSFP=config["dbNSFP"],
	shell:
		"vep -i {input.vcf} -o {output.vcf} --fork {threads} --compress_output bgzip --everything --offline --species homo_sapiens --stats_file {output.html} --assembly GRCh38 --cache --dir_cache {input.cache} --cache_version 110 --merged --fasta {params.ref} --format vcf --symbol --no_intergenic --merged --cache --pick --pick_allele --vcf --plugin dbNSFP,{params.dbNSFP},SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,FATHMM_converted_rankscore,FATHMM_pred --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN 2>{log} && tabix {output.vcf} 2>>{log}"

########################################################################### SPLIT VEP ###########################################################################

#MUTECT2

rule SplitVepMutect2SingleSamplePoN:
	input:
		vcf="results/{sample}_tumor/{sample}.mutect2.filtOnCtr.vep.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.mutect2.filtOnCtr.vep.vcf.gz.tbi",
		html="results/{sample}_tumor/{sample}.mutect2.filtOnCtr.vep.html",
	output:
		"results/{sample}_tumor/{sample}.mutect2.filtOnCtr.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/{sample}.SplitVepMutect2SingleSamplePoN.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule SplitVepMutect2MultisamplePoN:
	input:
		vcf="results/multisample.mutect2.vep.vcf.gz",
		tbi="results/multisample.mutect2.vep.vcf.gz.tbi",
		html="results/multisample.mutect2.vep.html"
	output:
		"results/multisample.mutect2.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/SplitVepMutect2MultisamplePoN.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

# FREEBAYES

rule SplitVepFreebayesSingleSamplePoN:
	input:
		vcf="results/{sample}_tumor/{sample}.freebayes.filtOnCtr.vep.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.freebayes.filtOnCtr.vep.vcf.gz.tbi",
		html="results/{sample}_tumor/{sample}.freebayes.filtOnCtr.vep.html",
	output:
		"results/{sample}_tumor/{sample}.freebayes.filtOnCtr.vep.tmp01.tsv"
	threads:1
	log:
		"logs/{sample}.SplitVepFreebayesSingleSamplePoN.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%RO\t%AO\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule SplitVepFreebayesMultisamplePoN:
	input:
		vcf="results/multisample.freebayes.vep.vcf.gz",
		tbi="results/multisample.freebayes.vep.vcf.gz.tbi",
		html="results/multisample.freebayes.vep.html"
	output:
		"results/multisample.freebayes.vep.tmp01.tsv"
	threads: 1 
	log:
		"logs/SplitVepFreebayesMultisamplePoN.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%RO\t%AO\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

# BCFTOOLS 

rule SplitVepBcftoolsSingleSamplePoN:
	input:
		vcf="results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vep.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vep.vcf.gz.tbi",
		html="results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vep.html",
	output:
		"results/{sample}_tumor/{sample}.bcftools.filtOnCtr.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/{sample}.SplitVepBcftoolsSingleSamplePoN.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule SplitVepBcftoolsMultisamplePoN:
	input:
		vcf="results/multisample.bcftools.vep.vcf.gz",
		tbi="results/multisample.bcftools.vep.vcf.gz.tbi",
		html="results/multisample.bcftools.vep.html"
	output:
		"results/multisample.bcftools.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/SplitVepBcftoolsMultisamplePoN.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

# VARSCAN 

rule SplitVepVarscanSingleSamplePoN:
	input:
		vcf="results/{sample}_tumor/{sample}.varscan.filtOnCtr.vep.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.varscan.filtOnCtr.vep.vcf.gz.tbi",
		html="results/{sample}_tumor/{sample}.varscan.filtOnCtr.vep.html",
	output:
		"results/{sample}_tumor/{sample}.varscan.filtOnCtr.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/{sample}.SplitVepVarscanSingleSamplePoN.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%RD\t%AD\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule SplitVepVarscanMultisamplePoN:
	input:
		vcf="results/multisample.varscan.vep.vcf.gz",
		tbi="results/multisample.varscan.vep.vcf.gz.tbi",
		html="results/multisample.varscan.vep.html"
	output:
		"results/multisample.varscan.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/SplitVepVarscanMultisamplePoN.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%RD\t%AD\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""
