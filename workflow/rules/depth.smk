rule Mosdepth:
    input:
        bam="alignments/{sample}.tumor.dd.rec.bam",
        bai="alignments/{sample}.tumor.dd.rec.bai"
    output:
        "alignments/{sample}.regions.bed.gz"
    log:
        "logs/{sample}.mosdepth.log"
    conda: "../envs/mosdepth.yaml"
    threads: 10
    params:
        prefix = "alignments/{sample}",
        target= config["intervals"]
    shell:
        "mosdepth --fast-mode -t {threads} --by {params.target} {params.prefix} {input.bam} 2> {log}"
