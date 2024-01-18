configfile: "config/conf.yaml"

include: "rules/fastqc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/alfred.smk"
include: "rules/paired.smk"
include: "rules/mutect2.smk"
#include: "rules/fastqc.smk"
#include: "rules/depth.smk"
#include: "rules/annotation.smk"
include: "rules/freebayes.smk"
#include: "rules/varscan.smk"
#include: "rules/vep_annotation.smk"
#include: "rules/bcftools.smk"

rule fastqc:
	input:
		expand(f"qc/{{sample}}_tumor_{{strand}}_fastqc.html", sample=config["tumors"].values(), strand=config["strand"].values()),
		expand(f"qc/{{sample}}_control_{{strand}}_fastqc.html", sample=config["controls"].values(), strand=config["strand"].values()),
		expand(f"qc/{{sample}}_tumor_{{strand}}_tr_fastqc.html", sample=config["tumors"].values(), strand=config["strand"].values()),
		expand(f"qc/{{sample}}_control_{{strand}}_tr_fastqc.html", sample=config["controls"].values(), strand=config["strand"].values()),

rule align:
	input:
		expand(f"alignments/{{sample}}.tumor.dd.rec.bam", sample=config["tumors"].values()),
		expand(f"alignments/{{sample}}.control.dd.rec.bam", sample=config["controls"].values())

rule bamqc:
	input:
		expand(f"qc/{{sample}}.tumor.qc.tsv.gz.pdf", sample=config["tumors"].values()),
		expand(f"qc/{{sample}}.control.qc.tsv.gz.pdf", sample=config["controls"].values())

rule mutect2:
	input:
		expand(f"results/{{sample}}_tumor/{{sample}}.mutect2.filtered.vcf.gz", sample=config["tumors"].values()),
		expand(f"results/{{sample}}_control/{{sample}}.mutect2.filtered.vcf.gz", sample=config["controls"].values())

rule freebayes:
	input:
		expand(f"results/{{sample}}_tumor/{{sample}}.freebayes.filtOnCtr.vcf.gz", sample=config["tumors"].values()),
		expand(f"results/{{sample}}_control/{{sample}}.freebayes.filt.vcf.gz", sample=config["controls"].values())
		#expand(f"results/freebayes/custom_pon/{{control}}.ctr.freebayes.filt.vcf.gz", control=config["controls"].values()),
		#expand(f"results/freebayes/{{sample}}.freebayes.filtOnCtr.vep.tsv", sample={**config["samples"],**config["controls"]}.values())

#rule varscan: 
#	input:
#		expand(f"results/varscan/custom_pon/{{control}}.ctr.varscan.filt.vcf.gz", control=config["controls"].values()),
#		expand(f"results/varscan/{{sample}}.sample.varscan.filtOnCtr.vcf.gz", sample=config["samples"].values()),
		#f"results/varscan/multisample.varscan.vcf.gz"

#rule bcftools:
#	input:
		#expand(f"results/variants/{{sample}}.bcftools.addaf.vcf.gz", sample={**config["tumors"],**config["normals"]}.values()),
#		expand(f"results/controls/{{ctr}}.bcftools.filt.vcf.gz", ctr=config["controls"].values()),
#		expand(f"results/{{tum}}/{{tum}}.bcftools.filtOnCtr.vcf.gz", tum=config["tumors"].values()),
		#expand(f"results/bcftools/custom_pon/{{control}}.ctr.bcftools.filt.vcf.gz", control=config["controls"].values()),
		#expand(f"results/bcftools/{{sample}}.sample.bcftools.filtOnCtr.vcf.gz", sample=config["samples"].values()),
		#f"results/bcftools/multisample.bcftools.vcf.gz"


rule paired:
	input:
		expand(f"qc/{{sample}}_{{type}}_{{strand}}_fastqc.html", sample=config["tumors"].values(), type=config["paired"].values(), strand=config["strand"].values()),
		expand(f"qc/{{sample}}_{{type}}_{{strand}}_tr_fastqc.html", sample=config["tumors"].values(), type=config["paired"].values(), strand=config["strand"].values()),
		expand(f"alignments/{{sample}}.{{type}}.dd.rec.bam", sample=config["tumors"].values(),type=config["paired"].values()),
		expand(f"qc/{{sample}}.{{type}}.qc.tsv.gz.pdf", sample=config["tumors"].values(),type=config["paired"].values()),
		expand(f"results/{{sample}}.mutect2.paired.filtered.vcf.gz",sample=config["tumors"].values()),
		expand(f"results/{{sample}}.varscan.paired.vcf.gz",sample=config["tumors"].values())


#rule vep_annotation: 
#	input:
#		f"results/varscan/multisample.varscan.vep.tsv"

#---------------------------
#rule qc:
#	input:
#		expand(f"results/fastqc/{{sample}}.{{strand}}_fastqc.html",sample=config["samples"].values(),strand=config["strand"].values())

#rule trimming:
#	input:
#		expand(f"data/fastp/{{sample}}.R1.fastq.gz", sample=config["samples"].values()),
#		expand(f"data/fastp/{{sample}}.R2.fastq.gz", sample=config["samples"].values())
		#expand(f"/home/simone/mnt/part1/wesMastobis/data/fastq/{{sample}}.R2.fastq.gz", sample=config["samples"].values())

#rule align:
#	input:
#		expand(f"data/bam/{{sample}}.markdup.bqsr.bam",sample=config["samples"].values())

#rule snpcalling:
#	input:
#		expand(f"results/mutect2/{{sample}}.mutect2.filtered.vcf.gz",sample=config["samples"].values())

#rule depth:
#	input:
#		expand(f"data/bam/{{sample}}.mosdepth.summary.txt",sample=config["samples"].values())

#rule annotation:
#	input: 
#		f"{mutectdir}/multisample.mutect2.annot.vcf"
