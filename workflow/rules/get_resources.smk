
## RUN samtools faidx referece.fasta

rule create_dict:
    input:
        config['genome'],
    output:
        "/home/simone/mnt/part1/smk_wes/resouces/GRCh38_full_analysis_set_plus_decoy_hla.dict"
    log:
        "logs/dict/genome.dict.log",
    params:
         #java_opts="-XX:ParallelGCThreads=10"
    resources:
        mem_mb=4096
    wrapper:
        "v3.3.3/bio/picard/createsequencedictionary"
