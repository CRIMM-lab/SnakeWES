rule snpEff: 
    input: 
        f"{mutectdir}/multisample.mutect2.vcf"
    output:
        temp(f"{mutectdir}/multisample.snpEff.vcf") 
    #threads: 
    params: 
        stats=f"{mutectdir}/multisample.snpEff_summary.html",
        intervals=config["intervals"],
        snpEff="/home/alessio/programs/snpEff"
    log: 
        "logs/snpEff.log"
    shell: 
        "java -Xmx20G -jar {params.snpEff}/snpEff.jar eff -v -stats {params.stats} -fi {params.intervals} GRCh38.105 {input} > {output}"

rule snpEff_norm: 
    input: 
        f"{mutectdir}/multisample.snpEff.vcf"
    output: 
        temp(f"{mutectdir}/multisample.snpEff.norm.vcf")
    #threads: 
    params: 
        snpEff="/home/alessio/programs/snpEff",
        header=f"{mutectdir}/h.txt"
    log: 
        "logs/norm.log"
    shell:
        """
        bcftools view -h {input} > {params.header} && cat {params.header} <(bcftools view -H {input} | awk 'OFS=FS="\t"''{split($8,info,"ANN=");split(info[2],ann,","); for (i in ${{ann}}) {print $1,$2,$3,$4,$5,$6,$7,info[1]"ANN="ann[${{i}}],$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19}}') > {output}
        """

rule dbsnp: 
    input: 
        f"{mutectdir}/multisample.snpEff.norm.vcf"
    output: 
        temp(f"{mutectdir}/multisample.snpEff.dbsnp.vcf")
    #threads: 
    params:
        snpEff="/home/alessio/programs/snpEff",
        dbsnp156="/home/simone/mnt/part1/resources/snpSift/hg38/dbSNP.b156.filt.vcf.gz"
    log: 
        "logs/dbsnp.log"
    shell: 
        "java -Xmx20G -jar {params.snpEff}/SnpSift.jar annotate -tabix -v {params.dbsnp156} {input} > {output}"


rule clnvar: 
    input:
        f"{mutectdir}/multisample.snpEff.dbsnp.vcf" 
    output:
        temp(f"{mutectdir}/multisample.snpEff.dbsnp.clnvar.vcf") 
    #threads: 
    params:
        snpEff="/home/alessio/programs/snpEff",
        clinvar="/home/simone/mnt/part1/resources/snpSift/hg38/clinvar_20230410.vcf.gz"
    log:
        "logs/clnvar.log" 
    shell: 
        "java -Xmx20G -jar {params.snpEff}/SnpSift.jar annotate -v {params.clinvar} {input} > {output}"

rule gnomAD: 
    input:
        f"{mutectdir}/multisample.snpEff.dbsnp.clnvar.vcf" 
    output:
        temp(f"{mutectdir}/multisample.snpEff.dbsnp.clnvar.gnomAD.vcf") 
    #threads: 
    params: 
        snpEff="/home/alessio/programs/snpEff",
        gnomAD="/home/simone/mnt/part1/resources/snpSift/hg38/af-only-gnomad.hg38.vcf.gz"
    log:
        "logs/gnomAD.log" 
    shell: 
        "java -Xmx20G -jar {params.snpEff}/SnpSift.jar annotate -v {params.gnomAD} {input} > {output}"

rule cosmic: 
    input:
        f"{mutectdir}/multisample.snpEff.dbsnp.clnvar.gnomAD.vcf" 
    output:
        temp(f"{mutectdir}/multisample.snpEff.dbsnp.clnvar.gnomAD.cosmic.vcf") 
    #threads: 
    params:
        snpEff="/home/alessio/programs/snpEff", 
        cosmic="/home/simone/mnt/part1/resources/snpSift/hg38/CosmicCodingMuts.vcf.gz"
    log:
        "logs/cosmic.log" 
    shell: 
        "java -Xmx20G -jar {params.snpEff}/SnpSift.jar annotate -v {params.cosmic} {input} > {output}" 

rule dbNSFP: 
    input:
        f"{mutectdir}/multisample.snpEff.dbsnp.clnvar.gnomAD.cosmic.vcf" 
    output:
        f"{mutectdir}/multisample.mutect2.annot.vcf" 
    #threads: 
    params:
        snpEff="/home/alessio/programs/snpEff",
        dbNSFP="/home/simone/mnt/qnap/dbNSFP/dbNSFP4.3a_grch38.gz"
    log: 
    shell: 
        "java -Xmx20G -jar {params.snpEff}/SnpSift.jar DbNsfp -v -db {params.dbNSFP} {input} > {output}"
