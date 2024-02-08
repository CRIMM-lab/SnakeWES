#!/bin/bash

# Create somatic panel of normal 
project="/home/simone/mnt/part1/masto_wes"

# Set Variables
threads=10

# Results
mutectdir="${project}/results/mutect2" 

# Resources
refgen="/home/simone/mnt/part1/resources/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
intervals="/home/simone/mnt/part1/resources/snpSift/hg38/SureSelectV8.wes.mod.bed"
gnomAD="/home/simone/mnt/part1/resources/snpSift/hg38/af-only-gnomad.hg38.vcf.gz"

# Data
bamdir="${project}/data/bam" 
logs="${project}/data/logs" 

# Programs
gatk="/home/alessio/programs/gatk-4.4.0.0/gatk"
picard="/home/alessio/programs/picard/build/libs/picard.jar"

#controls=$(find ${bamdir} -type f -iname '*.markdup.bqsr.bam' -exec basename {} \; | grep -E '^[C][D]')

#for bam in ${controls}; do  
			
#	echo $bam
#	name=$(echo ${bam} | sed -e 's/.markdup.bqsr.bam//')

#	echo $name 
	# calling mutect2 in tumor-only mode
#	 ${gatk} --java-options "-Xmx30G" Mutect2 -R ${refgen} \ ## RUN Mutect2 WITH -max-mnp-distance 0  
#		-I ${bamdir}/${bam} \
#		-tumor ${name} \
#		--native-pair-hmm-threads ${threads} \
#		--max-reads-per-alignment-start 0 \
#		--germline-resource ${gnomAD} \
#		-L ${intervals} \
#		-O ${mutectdir}/${name}_for_pon.vcf.gz 

#	echo -e "${name}\t${mutectdir}/${name}_for_pon.vcf.gz" >> ${mutectdir}/cohort.sample_map
#done 

# Create a file with path to all vcf files 
#find ${mutectdir} -type f -iname '*_for_pon.vcf.gz' > ${mutectdir}/all_vcf_for_pon.args

#${gatk} GenomicsDBImport -R ${refgen} -L ${intervals} \
#	--genomicsdb-workspace-path ${mutectdir}/pon_db \
#	--sample-name-map ${mutectdir}/cohort.sample_map \
#	--tmp-dir ${mutectdir}\
#	--merge-input-intervals true \
#	--reader-threads ${threads}\

#Combine the normal calls using CreateSomaticPanelOfNormals.
#${gatk} CreateSomaticPanelOfNormals \
#   -vcfs ${mutectdir}/all_vcf_for_pon.args \
#   -O ${mutectdir}/pon.vcf.gz
