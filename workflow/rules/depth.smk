rule mosdepth: 
        input: 
                "data/bam/{sample}.markdup.bqsr.bam" 
        output:
        	"data/bam/{sample}.mosdepth.summary.txt"
        threads: 4
        params:
        	mosdepth="/home/simone/programs/mosdepth-0.3.1/mosdepth",
        	intervals=config["intervals"],
        log:
                "logs/{sample}.mosdepth.log"
        shell:
        	"{params.mosdepth} --by {params.intervals} {wildcards.sample} -t {threads} {input} 2> {log}"

rule plot_cov:
        input:
