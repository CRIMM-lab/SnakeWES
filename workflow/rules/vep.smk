rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species="homo_sapiens_merged",
        build="GRCh38",
        release="110",
    log:
        "logs/vep_cache.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v3.3.5/bio/vep/cache"


rule download_vep_plugins:
    output:
        directory("resources/vep/plugins")
    params:
        release=110
    wrapper:
        "v3.3.5/bio/vep/plugins"


rule vep_single_sample:
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
        "logs/{sample}.vep_{caller}_single_sample.log"
    conda:
        "../envs/vep.yaml"
    params:
        ref=config["genome"],
        clinvar=config["clinvar"],
        dbNSFP=config["dbNSFP"],
    shell:
        "vep -i {input.vcf} -o {output.vcf} --fork {threads} --compress_output bgzip --everything --offline --species homo_sapiens --stats_file {output.html} --assembly GRCh38 --cache --dir_cache {input.cache} --cache_version 110 --merged --fasta {params.ref} --format vcf --symbol --no_intergenic --merged --cache --pick --pick_allele --vcf --plugin dbNSFP,{params.dbNSFP},SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,FATHMM_converted_rankscore,FATHMM_pred --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN 2>{log} && tabix {output.vcf} 2>>{log}"

rule vep_multisample:
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
        "logs/vep_{caller}_multisample.log"
    conda:
        "../envs/vep.yaml"
    params:
        ref=config["genome"],
        clinvar=config["clinvar"],
        dbNSFP=config["dbNSFP"],
    shell:
        "vep -i {input.vcf} -o {output.vcf} --fork {threads} --compress_output bgzip --everything --offline --species homo_sapiens --stats_file {output.html} --assembly GRCh38 --cache --dir_cache {input.cache} --cache_version 110 --merged --fasta {params.ref} --format vcf --symbol --no_intergenic --merged --cache --pick --pick_allele --vcf --plugin dbNSFP,{params.dbNSFP},SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,FATHMM_converted_rankscore,FATHMM_pred --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN 2>{log} && tabix {output.vcf} 2>>{log}"


#rule split_vep_single_sample:
#    input:
#        vcf="results/{sample}_tumor/{sample}.{caller}.filtOnCtr.vep.vcf.gz",
#    output:
#        tsv="results/{sample}_tumor/{sample}.{caller}.filtOnCtr.vep.tmp.tsv"
#    threads:
#        config["threads"]
#    log:
#        "logs/{sample}.split_vep_single_sample.log"
#    conda:
#        "../envs/bcftools.yaml"
#    shell:
