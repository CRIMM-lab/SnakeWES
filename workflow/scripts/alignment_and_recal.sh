#!/bin/bash

echo "performing alignment and recalibration - WES"

# project
project="/home/simone/mnt/part1/masto_wes"
threads=20

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


# Data
fastqdir="/home/simone/mnt/qnap/Danilo/Biodiversa/CD3+/trimming_fastp"
bamdir="${project}/data/bam" && mkdir -p ${bamdir}
dupMetricsdir="${project}/data/dupMetrics" && mkdir -p ${dupMetricsdir}
logs="${project}/logs" && mkdir -p ${logs}


for fastq in $(find ${fastqdir} -type f -iname '*.R1.fastq.gz'); do 

	echo $fastq

	name=$(basename "$fastq" .R1.fastq.gz)

	echo $name

	R1="${fastqdir}/${name}.R1.fastq.gz"	 

	R2="${fastqdir}/${name}.R2.fastq.gz"

	echo "STEP 1: Map to reference using BWA-MEM and Sort Bam - ${name}"

	if [ -e ${bamdir}/${name}.markdup.bqsr.bam ] && [ -e ${bamdir}/${name}.markdup.bqsr.bai ]; then 	

		echo "BAM file already exists for $name, skipping alignment."

	else 

		echo "Aligning $name..."
	
		bwa mem -t ${threads} -M -R "@RG\tID:${name}_wesMasto\tSM:${name}\tDT:20230515\tPI:150\tLB:SureSelectV8\tPL:ILLUMINA\tCN:BIODIVERSA" ${refgen} ${R1} ${R2} | \
		samtools sort -@ ${threads} -n -T ${name} -O bam -o ${bamdir}/${name}.bam - 2> ${logs}/${name}_alignment.log

		if [ $? -eq 0 ]; then

			echo "STEP 2: Mark Duplicates and Sort - ${name}"

			${gatk} --java-options "-Xmx30G" MarkDuplicatesSpark --spark-master local[${threads}] -I ${bamdir}/${name}.bam -O ${bamdir}/${name}.markdup.bam \
			-M ${dupMetricsdir}/${name}.marked.metrics.txt 2> ${logs}/${name}_MarkDuplicatesSpark.log
		
		else

			echo "alignment or sorting failed for ${name}"
			exit 1 
		fi

		if [ $? -eq 0 ]; then
		
			echo "STEP 4: Base quality recalibration - ${name}"

			${gatk} --java-options "-Xmx30G -XX:ParallelGCThreads=${threads}" BaseRecalibrator -I ${bamdir}/${name}.markdup.bam -R ${refgen} -L ${intervals} --known-sites ${dbsnp156} \
			--known-sites ${oneKG} --known-sites ${known_indels} --known-sites ${gnomAD} -O ${dupMetricsdir}/${name}.recal.data.table 2> ${logs}/${name}_BaseRecalibrator.log && \

			${gatk} --java-options "-Xmx30G -XX:ParallelGCThreads=${threads}" ApplyBQSR -I ${bamdir}/${name}.markdup.bam -R ${refgen} -L ${intervals} \
			--bqsr-recal-file ${dupMetricsdir}/${name}.recal.data.table -O ${bamdir}/${name}.markdup.bqsr.bam 2> ${logs}/${name}_ApplyBQSR.log

		else
			
			echo "recalibration failed for ${name}"
			exit 1 
		fi

		[ -s ${bamdir}/${name}.markdup.bqsr.bam ] && rm ${bamdir}/${name}.bam ${bamdir}/${name}.markdup.bam ${bamdir}/${name}.markdup.bam.sbi ${bamdir}/${name}.markdup.bam.bai
	
	fi

done
