rule AlfredTumorStats:
	input:
		bam="alignments/{sample}.tumor.dd.rec.bam",
		bai="alignments/{sample}.tumor.dd.rec.bai"
	output:
		"qc/{sample}.tumor.qc.tsv.gz"
	threads: 1
	log:
		"logs/{sample}.tumor.alfred.log"
	conda:
		"../envs/alfred.yaml"
	params:
		ref=config["genome"],
		interval=config["intervals"]	
	shell:
		"alfred qc -r {params.ref} -b {params.interval} -o {output} {input.bam} 2>{log}"

rule AlfredControlStats:
	input:
		bam="alignments/{sample}.control.dd.rec.bam",
		bai="alignments/{sample}.control.dd.rec.bai"
	output:
		"qc/{sample}.control.qc.tsv.gz"
	threads: 1
	log:
		"logs/{sample}.control.alfred.log"
	conda:
		"../envs/alfred.yaml"
	params:
		ref=config["genome"],
		interval=config["intervals"]	
	shell:
		"alfred qc -r {params.ref} -b {params.interval} -o {output} {input.bam} 2>{log}"

rule AlfredGermlineStats:
	input:
		bam="alignments/{sample}.germline.dd.rec.bam",
		bai="alignments/{sample}.germline.dd.rec.bai"
	output:
		"qc/{sample}.germline.qc.tsv.gz"
	threads: 1
	log:
		"logs/{sample}.germline.alfred.log"
	conda:
		"../envs/alfred.yaml"
	params:
		ref=config["genome"],
		interval=config["intervals"]	
	shell:
		"alfred qc -r {params.ref} -b {params.interval} -o {output} {input.bam} 2>{log}"

rule PlotAlfredTumorStats:
	input:
		"qc/{sample}.tumor.qc.tsv.gz"
	output:
		"qc/{sample}.tumor.qc.tsv.gz.pdf"
	threads: 1
	log:
		"logs/{sample}.tumor.alfred.plot.log"
	conda:
		"../envs/alfred.yaml"
	params:
		script="workflow/scripts/stats.R"
	shell:
		"Rscript {params.script} {input} 2>{log}"

rule PlotAlfredControlStats:
	input:
		"qc/{sample}.control.qc.tsv.gz"
	output:
		"qc/{sample}.control.qc.tsv.gz.pdf"
	threads: 1
	log:
		"logs/{sample}.control.alfred.plot.log"
	conda:
		"../envs/alfred.yaml"
	params:
		script="workflow/scripts/stats.R"
	shell:
		"Rscript {params.script} {input} 2>{log}"

rule PlotAlfredGermlineStats:
	input:
		"qc/{sample}.germline.qc.tsv.gz"
	output:
		"qc/{sample}.germline.qc.tsv.gz.pdf"
	threads: 1
	log:
		"logs/{sample}.germline.alfred.plot.log"
	conda:
		"../envs/alfred.yaml"
	params:
		script="workflow/scripts/stats.R"
	shell:
		"Rscript {params.script} {input} 2>{log}"