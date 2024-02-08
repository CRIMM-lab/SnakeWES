#!/bin/bash

echo "performing calling with freebayes"

# project
project="/home/simone/mnt/part1/masto_wes"
threads=4

# Resources
refgen="/home/simone/mnt/part1/resources/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
intervals="/home/simone/mnt/part1/resources/snpSift/hg38/SureSelectV8.wes.mod.bed"
dbsnp156="/home/simone/mnt/part1/resources/snpSift/hg38/dbSNP.b156.filt.vcf.gz"
known_indels="/home/simone/mnt/part1/resources/snpSift/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
oneKG="/home/simone/mnt/part1/resources/snpSift/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
gnomAD="/home/simone/mnt/part1/resources/snpSift/hg38/af-only-gnomad.hg38.vcf.gz"

# Programs
gatk="/home/alessio/programs/gatk-4.4.0.0/gatk"
picard="/home/alessio/programs/picard/build/libs/picard.jar"
freebayes="/home/simone/programs/freebayes/build/freebayes"

# Data
fastqdir="/home/simone/mnt/qnap/Danilo/Biodiversa/CD3+/trimming_fastp"
bamdir="${project}/data/bam" && mkdir -p ${bamdir}
dupMetricsdir="${project}/data/dupMetrics" && mkdir -p ${dupMetricsdir}
logs="${project}/logs" && mkdir -p ${logs}
freebayesdir="/home/simone/mnt/part1/masto_wes/results/freebayes_with_ctr"

#samples=$(find ${bamdir} -type f -iname '*.markdup.bqsr.bam' -exec basename {} \; | grep -E '^(A|C)[0-9]{4}')

#for bam in ${samples}; do

#	name=$(echo ${bam} | sed -e 's/.markdup.bqsr.bam//')

#	if [ -s ${freebayes_with_ctr}/${name}.freebayes.vcf ]; then

#	    echo "${name} vcf already exist"

#	elif [ -s ${bamdir}/${bam} ]; then 

#		echo "calling ${name} with freebayes"

#	    ${freebayes} -f ${refgen} -F 0.01 -C 1 --pooled-continuous ${bamdir}/${bam} > ${freebayesdir}/${name}.freebayes.vcf

#	else 

#		echo "bam for ${name}"

#	fi 

#done

controls=$(find ${bamdir} -type f -iname '*.markdup.bqsr.bam' -exec basename {} \; | grep -E '^[C][D]')

run_freebayes() {
    refgen="/home/simone/mnt/part1/resources/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    freebayes="/home/simone/programs/freebayes/build/freebayes"
	bamdir="/home/simone/mnt/part1/masto_wes/data/bam"
	freebayesdir="/home/simone/mnt/part1/masto_wes/results/freebayes_with_ctr"

    bam="$1"
    name=$(echo ${bam} | sed -e 's/.markdup.bqsr.bam//')
    
    ${freebayes} -f ${refgen} -F 0.01 -C 1 --pooled-continuous ${bamdir}/${bam} > ${freebayesdir}/${name}.freebayes.vcf
}

# Export the function for use with parallel
export -f run_freebayes

# Read sample names from a file and run FreeBayes in parallel
echo $controls | tr ' ' '\n' | parallel -j $threads run_freebayes

