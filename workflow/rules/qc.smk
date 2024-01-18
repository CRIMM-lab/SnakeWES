rule fastqc:
    input:
        "data/fastq/{sample}.R1.tr.fastq.gz"
    output:
        html="qc/fastqc/{sample}.R1.tr.html",
        zip="qc/fastqc/{sample}.R1.tr.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}.log"
    threads: 1
    resources:
        mem_mb=4096 
    wrapper:
        "v3.3.3/bio/fastqc"


rule samtools_stats:
    input:
        bam="data/bam/{sample}.dd.rec.bam",
        bed=config['interval']
    output:
        "qc/samtools_stats/{sample}.txt"
    log:
        "logs/stats/{sample}.samtools.log",
    wrapper:
        "v3.3.3/bio/samtools/stats"
