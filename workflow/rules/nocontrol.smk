########################################################################### Mutect2 ###########################################################################

rule Mutect2TumorNocontrol: 
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
		"logs/{sample}.Mutect2TumorNocontrol.log",
	params:
		extra="--max-reads-per-alignment-start 0" 
	wrapper:
		"v3.3.3/bio/gatk/mutect"

rule OrientationModelTumorNocontrol:
	input:
		f1r2="data/{sample}.tumor.f1r2.tar.gz",
	output:
		"data/{sample}.tumor.artifacts_prior.tar.gz",
	resources:
		mem_mb=4096,
	log:
		"logs/{sample}.OrientationModelTumorNocontrol.log",
	wrapper:
		"v3.3.3/bio/gatk/learnreadorientationmodel"

rule GetPileupSummariesTumorNocontrol:
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
		"logs/{sample}.GetPileupSummariesTumorNocontrol.log",
	wrapper:
		"v3.3.3/bio/gatk/getpileupsummaries"

rule CalculateContaminationTumorNocontrol:
	input:
		"data/{sample}.tumor.pileupTable"
	output:
		contamination="data/{sample}.tumor.contaminationTable",
		segmentation="data/{sample}.tumor.tumorSeg.txt"
	threads: 1
	log:
		"logs/{sample}.CalculateContaminationTumorNocontrol.log",
	conda:
		"../envs/gatk4.yaml"
	params:
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"]),
		mem_mb="-Xmx4G"
		#extra="--tumor-segmentation data/{wildcards.sample}.tumseg.txt"
	shell:
		"gatk --java-options {params.mem_mb} CalculateContamination -I {input} -O {output.contamination} --tumor-segmentation {output.segmentation} > {log} 2>&1"


rule FilterMutectCallsTumorNocontrol:
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
		"logs/{sample}.FilterMutectCallsTumorNocontrol.log",
	params:
		#extra="--tumor-segmentation data/{wildcard.sample}.tumseg.txt",  # optional arguments, see GATK docs
		java_opts="-XX:ParallelGCThreads=" + str(config["threads"])  # optional
	resources:
		mem_mb=4096,
	wrapper:
		"v3.3.3/bio/gatk/filtermutectcalls"

rule CleanFilterMutectTumorOutputNocontrol:
	input:
		"results/{sample}_tumor/{sample}.mutect2.filt.vcf.gz"
	output:
		vcf="results/{sample}_tumor/{sample}.mutect2.filt.vep.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.mutect2.filt.vep.vcf.gz.tbi"
	threads: 1
	log:
		"logs/{sample}.CleanFilterMutectTumorOutputNocontrol.log"
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

rule MergeMutect2TumorOutputNocontrol:
	input:
		vcf=expand(f"results/{{sample}}_tumor/{{sample}}.mutect2.filt.vep.vcf.gz", sample=config["samples"].values()),
		tbi=expand(f"results/{{sample}}_tumor/{{sample}}.mutect2.filt.vep.vcf.gz.tbi", sample=config["samples"].values())
	output:
		vcf="results/multisample.mutect2.vcf.gz",
		tbi="results/multisample.mutect2.vcf.gz.tbi"
	threads: 1
	log:
		"logs/multisample.MergeMutect2TumorOutputNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2>{log} && "
		"tabix -p vcf {output.vcf}"


########################################################################### Bcftools ###########################################################################

rule BcftoolsCallOnTumorsNocontrol:  
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
		"logs/{sample}.BcftoolsCallOnTumorsNocontrol.log"
	shell:
		"bcftools mpileup -Ou -d 10000 -R {params.intervals} -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR -f {params.ref} {input.bam} | "
		"bcftools call -mv | "
		"bcftools norm -m - -f {params.ref} -Oz -o {output} 2>{log}"

rule BcftoolsAddAfFieldOnTumorsNocontrol:  
	input:
		"results/{sample}_tumor/{sample}.bcftools.vcf.gz"
	output:
		vcf=temp("results/{sample}_tumor/{sample}.bcftools.addaf.vcf.gz")
#		tsv=temp("results/{sample}_tumor/{sample}.annot.tsv.gz"),
#		tsv_tbi=temp("results/{sample}_tumor/{sample}.annot.tsv.gz.tbi")
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.BcftoolsAddAfFieldOnTumorsNoncontrol.log"
	params:
		tsv="results/{sample}_tumor/{sample}.annot.tsv.gz",
		tsv_tbi="results/{sample}_tumor/{sample}.annot.tsv.gz.tbi"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%AD{{0}}\t%AD{{1}}]' {input} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {params.tsv} 2>{log} && 
		tabix -b2 -e2 {params.tsv} 2>>{log} && 
		bcftools annotate -a {params.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF {input} -Oz -o {output.vcf} 2>>{log} && 
		rm {params.tsv} 2>>{log} && 
		rm {params.tsv_tbi} 2>>{log}
		"""

rule BcftoolsFilterTumorsNocontrol: 
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
		"logs/{sample}.BcftoolsFilterTumorsNocontrol.log"
	shell:
		"bcftools view -i 'FORMAT/DP >= {params.depth} & FORMAT/AF >= {params.vaf} & FORMAT/AD[0:1] >= {params.alt}' {input} | "
		"grep -v -f {params.excl} | "
		"bcftools sort -Oz -o {output.vcf} 2>{log} && "
		"tabix -p vcf {output.vcf} 2>>{log}"

rule MergeBcftoolsTumorOutputNocontrol:
	input:
		vcf=expand(f"results/{{sample}}_tumor/{{sample}}.bcftools.filt.vcf.gz", sample=config["samples"].values()),
		tbi=expand(f"results/{{sample}}_tumor/{{sample}}.bcftools.filt.vcf.gz.tbi", sample=config["samples"].values())
	output:
		vcf="results/multisample.bcftools.vcf.gz",
		tbi="results/multisample.bcftools.vcf.gz.tbi"
	threads: 1
	log:
		"logs/multisample.MergeBcftoolsTumorOutputNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2>{log} && "
		"tabix -p vcf {output.vcf} 2>>{log}"


########################################################################### Varscan ###########################################################################

rule VarscanCallSnvTumorsNocontrol: 
	input:
		bam="alignments/{sample}.tumor.dd.rec.bam",
		bai="alignments/{sample}.tumor.dd.rec.bai"
	output:
		vcf_snv=temp("results/{sample}_tumor/{sample}.varscan.snv.vcf.gz"),
		tbi_snv=temp("results/{sample}_tumor/{sample}.varscan.snv.vcf.gz.tbi")
	threads: 1
	params:
		intervals=config["intervals"],
		ref=config["genome"]
	conda:
		"../envs/varscan.yaml"
	log:
		"logs/{sample}.VarscanCallSnvTumorsNocontrol.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} |"
		"varscan mpileup2snp --output-vcf --min-avg-qual 1 --p-value 1 --min-var-freq 0 --min-coverage 1 --min-reads2 1 |"
		"bgzip -c > {output.vcf_snv} 2>{log} &&"
		"tabix -p vcf {output.vcf_snv} 2>>{log}"

rule VarscanCallIndelTumorsNocontrol: 
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
		"logs/{sample}.VarscanCallIndelTumorsNocontrol.log"
	shell:
		"samtools mpileup -l {params.intervals} -f {params.ref} {input.bam} | "
		"varscan mpileup2indel --output-vcf --min-avg-qual 1 --p-value 1 --min-var-freq 0 --min-coverage 1 --min-reads2 1 | "
		"bgzip -c > {output.vcf_indel} 2>{log} && "
		"tabix -p vcf {output.vcf_indel} 2>>{log}"

rule VarscanConcatTumorNocontrol: 
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
		"logs/{sample}.VarscanConcatTumorNocontrol.log"
	shell: 
		"bcftools concat -a -Ov {input.vcf_snv} {input.vcf_indel} | "
		"bcftools norm -Oz -m - -f {params} -Oz -o {output} 2>{log}"

rule VarscanModAFTumorNoncontrol:
	input:
		"results/{sample}_tumor/{sample}.varscan.concat.tmp.vcf.gz"
	output:
		vcf="results/{sample}_tumor/{sample}.varscan.concat.vcf.gz",
		tsv="results/{sample}_tumor/{sample}.varscan.annot.tsv.gz",
		tsv_tbi="results/{sample}_tumor/{sample}.varscan.annot.tsv.gz.tbi"
	log:
		"logs/{sample}.VarscanModAFTumorNoncontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%FREQ]' {input} | awk -F'\t' 'BEGIN {{OFS="\t"}} {{ $5 = $5 / 100; print }}'| bgzip -c > {output.tsv} 2>{log} && 
		tabix -b2 -e2 {output.tsv} 2>>{log} &&
		bcftools annotate -x FORMAT/FREQ -a {output.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF  {input} -Oz -o {output.vcf} 2>>{log}
		"""

rule VarscanReheaderTumorNontrol:
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
		"logs/{sample}.VarscanReheaderTumorNontrol.log"
	shell: 
		"echo {wildcards.sample} > {params.txt} 2>{log} && "
		"bcftools reheader -s {params.txt} -o {output} {input} 2>>{log} && "
		"rm {params.txt} 2>>{log}"

rule FilterVarscanTumorsNocontrol: 
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
		vaf=config["filtering_tumors"]["vaf"]
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.FilterVarscanTumorsNocontrol.log"
	shell:
		"bcftools view -Ov -i'FORMAT/DP >= {params.min_depth} & FORMAT/AD >= {params.alt} & FORMAT/AF >= {params.vaf}' {input} | "
		"grep -v -f {params.excl} | "
		"bcftools sort -Oz -o {output.vcf} 2>{log} && "
		"tabix -f -p vcf {output.vcf} 2>>{log}"

rule MergeVarscanTumorsNocontrol: ## difference between single curly and doucle curly brackets ##
	input:
		vcf=expand(f"results/{{sample}}_tumor/{{sample}}.varscan.filt.vcf.gz", sample=config['samples'].values()),
		tbi=expand(f"results/{{sample}}_tumor/{{sample}}.varscan.filt.vcf.gz.tbi", sample=config['samples'].values())
	output:
		vcf="results/multisample.varscan.vcf.gz",
		tbi="results/multisample.varscan.vcf.gz.tbi"
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	log: 
		"logs/multisample.varscan.log"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2>{log} &&"
		"tabix -p vcf {output.vcf} 2>>{log}"

########################################################################### Freebayes ###########################################################################

rule FreebayesCallNormTumorNocontrol: 
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
		"logs/{sample}.freebayes_call_and_norm_on_samples.log"
	shell:
		"freebayes -f {params.ref} -F 0.01 -C 2 -t {params.intervals} --pooled-continuous {input.bam} | "
		"bcftools norm -m - -f {params.ref} -Oz -o {output} 2>{log}"

rule FreebayesAddAFTumorNocontrol:
	input:
		"results/{sample}_tumor/{sample}.freebayes.vcf.gz"
	output:
		vcf="results/{sample}_tumor/{sample}.freebayes.annotAF.vcf.gz",
		#tsv="results/{sample}_tumor/{sample}.freebayes.annot.tsv.gz",
		#tsv_tbi="results/{sample}_tumor/{sample}.freebayes.annot.tsv.gz.tbi"
	log:
		"logs/{sample}.FreebayesAddAFTumorNocontrol.log",
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	params:
		tsv="results/{sample}_tumor/{sample}.freebayes.annot.tsv.gz",
		tsv_tbi="results/{sample}_tumor/{sample}.freebayes.annot.tsv.gz.tbi"		
	shell:
		"""
		bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%RO\t%AO\n' {input} 2>{log} | 
		awk 'OFS=FS="\t"''{{print $1,$2,$3,$4,$6/($5 + $6)}}' | bgzip -c > {params.tsv} 2>>{log} && 
		tabix -b2 -e2 {params.tsv} 2>>{log} && 
		bcftools annotate -x INFO/AF -a {params.tsv} -h <(echo '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">') --columns CHROM,POS,REF,ALT,FORMAT/AF -Oz -o {output.vcf} {input} 2>>{log} &&
		rm {params.tsv} 
		rm {params.tsv_tbi}
		"""

rule FreebayesFilterTumorNocontrol: 
	input:
		"results/{sample}_tumor/{sample}.freebayes.annotAF.vcf.gz"
	output:
		vcf="results/{sample}_tumor/{sample}.freebayes.filt.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.freebayes.filt.vcf.gz.tbi"
	log:
		"logs/{sample}.FreebayesFilterTumorNocontrol.log",
	threads: 1
	conda:
		"../envs/bcftools.yaml"
	params:
		excl=config["chr_to_exclude"],
		depth=config['filtering_tumors']['min_depth'],
		vaf=config['filtering_tumors']['vaf'],
		alt=config['filtering_tumors']['alt_depth']
	shell:
		"bcftools view -i 'QUAL > 1 & INFO/DP >= {params.depth} & FORMAT/AF >= {params.vaf} & INFO/AC >= {params.alt}' {input} | "
		"grep -v -f {params.excl} | "
		"bcftools sort -Oz -o {output.vcf} 2> {log} && "
		"tabix -p vcf {output.vcf} 2>>{log}"

rule MergeFreebayesTumorOutputNocontrol:
	input:
		vcf=expand(f"results/{{sample}}_tumor/{{sample}}.freebayes.filt.vcf.gz", sample=config["samples"].values()),
		tbi=expand(f"results/{{sample}}_tumor/{{sample}}.freebayes.filt.vcf.gz.tbi", sample=config["samples"].values())
	output:
		vcf="results/multisample.freebayes.vcf.gz",
		tbi="results/multisample.freebayes.vcf.gz.tbi"
	threads: 1
	log:
		"logs/multisample.MergeFreebayesTumorOutputNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools merge -m none -Oz -o {output.vcf} {input.vcf} 2>{log} && "
		"tabix -p vcf {output.vcf} 2>>{log}"


########################################################################### VEP ###########################################################################

rule VepSingleSampleNocontrol:
	input:
		vcf="results/{sample}_tumor/{sample}.{caller}.filt.vcf.gz",
		cache=config["vepcache"]
	output:
		vcf="results/{sample}_tumor/{sample}.{caller}.filt.vep.vcf.gz",  
		tbi="results/{sample}_tumor/{sample}.{caller}.filt.vep.vcf.gz.tbi",
		html="results/{sample}_tumor/{sample}.{caller}.filt.vep.html"
	threads: 
		config["threads"]
	log:
		"logs/{sample}.Vep{caller}SingleSampleNocontrol.log"
	conda:
		"../envs/vep.yaml"
	params:
		ref=config["genome"],
		clinvar=config["clinvar"],
		dbNSFP=config["dbNSFP"],
	shell:
		"vep -i {input.vcf} -o {output.vcf} --fork {threads} --compress_output bgzip --everything --offline --species homo_sapiens --stats_file {output.html} --assembly GRCh38 --cache --dir_cache {input.cache} --cache_version 110 --merged --fasta {params.ref} --format vcf --symbol --no_intergenic --merged --cache --pick --pick_allele --vcf --plugin dbNSFP,{params.dbNSFP},SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,FATHMM_converted_rankscore,FATHMM_pred --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN 2>{log} && tabix {output.vcf} 2>>{log}"

rule VepMultisampleNocontrol:
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

rule SplitVepMutect2SingleSampleNocontrol:
	input:
		vcf="results/{sample}_tumor/{sample}.mutect2.filt.vep.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.mutect2.filt.vep.vcf.gz.tbi",
		html="results/{sample}_tumor/{sample}.mutect2.filt.vep.html",
	output:
		"results/{sample}_tumor/{sample}.mutect2.filt.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/{sample}.SplitVepMutect2SingleSampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule SplitVepMutect2MultisampleNocontrol:
	input:
		vcf="results/multisample.mutect2.vep.vcf.gz",
		tbi="results/multisample.mutect2.vep.vcf.gz.tbi",
		html="results/multisample.mutect2.vep.html"
	output:
		"results/multisample.mutect2.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/SplitVepMutect2MultisampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

# FREEBAYES

rule SplitVepFreebayesSingleSampleNocontrol:
	input:
		vcf="results/{sample}_tumor/{sample}.freebayes.filt.vep.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.freebayes.filt.vep.vcf.gz.tbi",
		html="results/{sample}_tumor/{sample}.freebayes.filt.vep.html",
	output:
		"results/{sample}_tumor/{sample}.freebayes.filt.vep.tmp01.tsv"
	threads:1
	log:
		"logs/{sample}.SplitVepFreebayesSingleSampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%RO\t%AO\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule SplitVepFreebayesMultisampleNocontrol:
	input:
		vcf="results/multisample.freebayes.vep.vcf.gz",
		tbi="results/multisample.freebayes.vep.vcf.gz.tbi",
		html="results/multisample.freebayes.vep.html"
	output:
		"results/multisample.freebayes.vep.tmp01.tsv"
	threads: 1 
	log:
		"logs/SplitVepFreebayesMultisampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%RO\t%AO\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

# BCFTOOLS 

rule SplitVepBcftoolsSingleSampleNocontrol:
	input:
		vcf="results/{sample}_tumor/{sample}.bcftools.filt.vep.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.bcftools.filt.vep.vcf.gz.tbi",
		html="results/{sample}_tumor/{sample}.bcftools.filt.vep.html",
	output:
		"results/{sample}_tumor/{sample}.bcftools.filt.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/{sample}.SplitVepBcftoolsSingleSampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule SplitVepBcftoolsMultisampleNocontrol:
	input:
		vcf="results/multisample.bcftools.vep.vcf.gz",
		tbi="results/multisample.bcftools.vep.vcf.gz.tbi",
		html="results/multisample.bcftools.vep.html"
	output:
		"results/multisample.bcftools.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/SplitVepBcftoolsMultisampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%AD{{0}}\t%AD{{1}}\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

# VARSCAN 

rule SplitVepVarscanSingleSampleNocontrol:
	input:
		vcf="results/{sample}_tumor/{sample}.varscan.filt.vep.vcf.gz",
		tbi="results/{sample}_tumor/{sample}.varscan.filt.vep.vcf.gz.tbi",
		html="results/{sample}_tumor/{sample}.varscan.filt.vep.html",
	output:
		"results/{sample}_tumor/{sample}.varscan.filt.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/{sample}.SplitVepVarscanSingleSampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%RD\t%AD\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""

rule SplitVepVarscanMultisampleNocontrol:
	input:
		vcf="results/multisample.varscan.vep.vcf.gz",
		tbi="results/multisample.varscan.vep.vcf.gz.tbi",
		html="results/multisample.varscan.vep.html"
	output:
		"results/multisample.varscan.vep.tmp01.tsv"
	threads: 1
	log:
		"logs/SplitVepVarscanMultisampleNocontrol.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%CSQ\t[%RD\t%AD\t%AF][\t%GT]\n" -d -A tab > {output} 2>{log}
		"""
