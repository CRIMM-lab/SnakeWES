########## WHOLE GENOME SEQUENCING ANALYSIS WITH GATK #######################

# Define project

project="/home/simone/mnt/part1/masto_wes"
mkdir -p ${project}/data
mkdir -p ${project}/docs
mkdir -p ${project}/results
mkdir -p ${project}/resources
mkdir -p ${project}/scripts

# Set Variables

# Results
mutectdir="${project}/results/mutect2" && mkdir -p ${mutectdir}
freebayesdir="${project}/results/freebayes" && mkdir -p ${freebayesdir}

# Resources
refgen="/home/simone/mnt/part1/resources/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
intervals="/home/simone/mnt/part1/resources/snpSift/hg38/SureSelectV8.wes.mod.bed"
dbsnp156="/home/simone/mnt/part1/resources/snpSift/hg38/dbSNP.b156.filt.vcf.gz"
#dbsnp138="/home/simone/mnt/part1/resources/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
known_indels="/home/simone/mnt/part1/resources/snpSift/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
oneKG="/home/simone/mnt/part1/resources/snpSift/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
gnomAD="/home/simone/mnt/part1/resources/snpSift/hg38/af-only-gnomad.hg38.vcf.gz"
clinvar="/home/simone/mnt/part1/resources/snpSift/hg38/clinvar_20230410.vcf.gz"
dbNSFP="/home/simone/mnt/qnap/dbNSFP/dbNSFP4.3a_grch38.gz"
PoN="/home/simone/mnt/part1/resources/snpSift/hg38/1000g_pon.hg38.vcf.gz"
cosmic="/home/simone/mnt/part1/resources/snpSift/hg38/CosmicCodingMuts.vcf.gz"
chrToExclude="${project}/resources/GRCh38_full_analysis_set_plus_decoy_hla.exclude.txt"

# Scripts
countSamples="${project}/scripts/countSamples.sh"

# Data
fastqdir="${project}/data/fastq" && mkdir -p ${fastqdir}/untrimmed
bamdir="${project}/data/bam" && mkdir -p ${bamdir}
dupMetricsdir="${project}/data/dupMetrics" && mkdir -p ${dupMetricsdir}
fastpdir="${project}/data/fastp" && mkdir -p ${fastpdir}
#fastQCdir="${project}/data/fastQC" && mkdir -p ${fastQCdir}/{trimmed,untrimmed}
stderrdir="${project}/data/stderr" && mkdir -p ${stderrdir}
coverage="${project}/data/coverage"

# Programs
gatk="/home/alessio/programs/gatk-4.4.0.0/gatk"
fastp="/home/simone/programs/fastp/fastp"
picard="/home/alessio/programs/picard/build/libs/picard.jar"
mosdepth="/home/simone/programs/mosdepth-0.3.1/mosdepth"
snpEff="/home/alessio/programs/snpEff"
#fastQC="/home/alessio/programs/FastQC/fastqc"
freebayes="/home/simone/programs/freebayes/build/freebayes"

# Trimming with fastp

#R1=("${fastqdir}/untrimmed/A0934_SA_L001_R1_001.fastq.gz") #"${fastqdir}/untrimmed/A0030_SA_L001_R1_001.fastq.gz" "${fastqdir}/untrimmed/C3068_SA_L001_R1_001.fastq.gz" "${fastqdir}/untrimmed/A0789_SA_L001_R1_001.fastq.gz" "${fastqdir}/untrimmed/C3142_SA_L001_R1_001.fastq.gz" "${fastqdir}/untrimmed/A1428_SA_L001_R1_001.fastq.gz" "${fastqdir}/untrimmed/C1294_SA_L001_R1_001.fastq.gz" "${fastqdir}/untrimmed/A0734_SA_L001_R1_001.fastq.gz" "${fastqdir}/untrimmed/A0340_SA_L001_R1_001.fastq.gz" "${fastqdir}/untrimmed/A1018_SA_L001_R1_001.fastq.gz") 
#R2=("${fastqdir}/untrimmed/A0934_SA_L001_R2_001.fastq.gz") #"${fastqdir}/untrimmed/A0030_SA_L001_R2_001.fastq.gz" "${fastqdir}/untrimmed/C3068_SA_L001_R2_001.fastq.gz" "${fastqdir}/untrimmed/A0789_SA_L001_R2_001.fastq.gz" "${fastqdir}/untrimmed/C3142_SA_L001_R2_001.fastq.gz" "${fastqdir}/untrimmed/A1428_SA_L001_R2_001.fastq.gz" "${fastqdir}/untrimmed/C1294_SA_L001_R2_001.fastq.gz" "${fastqdir}/untrimmed/A0734_SA_L001_R2_001.fastq.gz" "${fastqdir}/untrimmed/A0340_SA_L001_R2_001.fastq.gz" "${fastqdir}/untrimmed/A1018_SA_L001_R2_001.fastq.gz")

#for ((i=0; i<${#R1[@]}; i++))
#do

#	echo ${R1[$i]}

#	echo "Extract the sample name from the R1 file-path"
#	name=$(basename ${R1[$i]} | sed -e "s/_SA_L001_R1_001.fastq.gz//")

#	echo ${name}

#	${fastp} -i ${R1[$i]} -I ${R2[$i]} -o ${fastqdir}/${name}.R1.fastq.gz -O ${fastqdir}/${name}.R2.fastq.gz -j ${fastpdir}/${name}.json -h ${fastpdir}/${name}.html -w 20 --failed_out ${fastpdir}/${name}.failedreads.txt --length_required 30 --detect_adapter_for_pe --disable_quality_filtering 2> ${stderrdir}/${name}_fastp_stderr.log

#	echo "STEP 1: Map to reference using BWA-MEM and Sort Bam - ${name}"

#	bwa mem -t 40 -M -R "@RG\tID:${name}_wesMasto\tSM:${name}\tDT:20230515\tPI:150\tLB:SureSelectV8\tPL:ILLUMINA\tCN:BIODIVERSA" ${refgen} ${fastqdir}/${name}.R1.fastq.gz ${fastqdir}/${name}.R2.fastq.gz \
#	| samtools sort -@ 40 -n -T ${name} -O bam -o ${bamdir}/${name}.bam - 2> ${stderrdir}/${name}_alignment_stderr.log

#	echo "STEP 2: Mark Duplicates and Sort - ${name}"

#	${gatk} --java-options "-Xmx30G" MarkDuplicatesSpark --spark-master local[40] -I ${bamdir}/${name}.bam -O ${bamdir}/${name}.markdup.bam -M ${dupMetricsdir}/${name}.marked.metrics.txt 2> ${stderrdir}/${name}_MarkDuplicatesSpark_stderr.log

#	echo "STEP 4: Base quality recalibration - ${name}"

#	${gatk} --java-options "-Xmx30G -XX:ParallelGCThreads=40" BaseRecalibrator -I ${bamdir}/${name}.markdup.bam -R ${refgen} -L ${intervals} --known-sites ${dbsnp156} --known-sites ${oneKG} --known-sites ${known_indels} --known-sites ${gnomAD} -O ${dupMetricsdir}/${name}.recal.data.table 2> ${stderrdir}/${name}_BaseRecalibrator_stderr.log

#	${gatk} --java-options "-Xmx30G -XX:ParallelGCThreads=40" ApplyBQSR -I ${bamdir}/${name}.markdup.bam -R ${refgen} -L ${intervals} --bqsr-recal-file ${dupMetricsdir}/${name}.recal.data.table -O ${bamdir}/${name}.markdup.bqsr.bam 2> ${stderrdir}/${name}_ApplyBQSR_stderr.log

#	echo "STEP 5: Call Variants with Mutect2 in tumor only mode - ${name}"

#	${gatk} --java-options "-Xmx30G" Mutect2 -R ${refgen} -I ${bamdir}/${name}.markdup.bqsr.bam --native-pair-hmm-threads 40 -L ${intervals} -O ${mutectdir}/${name}.mutect2.vcf.gz --panel-of-normals ${PoN} --germline-resource ${gnomAD} --bam-output ${bamdir}/${name}.mutect2.bam --max-reads-per-alignment-start 0 --f1r2-tar-gz ${mutectdir}/${name}.f1r2.tar.gz 2> ${stderrdir}/${name}_Mutect2__stderr.log

#	${gatk} --java-options "-Xmx30G" LearnReadOrientationModel -I ${mutectdir}/${name}.f1r2.tar.gz -O ${mutectdir}/${name}.read.orientation.model.tar.gz 2> ${stderrdir}/${name}_OrModel_stderr.log

#	${gatk} --java-options "-Xmx30G" GetPileupSummaries -I ${bamdir}/${name}.markdup.bqsr.bam -L ${intervals} -V ${gnomAD} -O ${mutectdir}/${name}.pileups.table 2> ${stderrdir}/${name}_Pileup_stderr.log

#	${gatk} --java-options "-Xmx30G" CalculateContamination -I ${mutectdir}/${name}.pileups.table -O ${mutectdir}/${name}.contamination.table --tumor-segmentation ${mutectdir}/${name}.tumseg.txt 2> ${stderrdir}/${name}_calcCont_stderr.log

#	${gatk} --java-options "-Xmx30G" FilterMutectCalls -R ${refgen} -V ${mutectdir}/${name}.mutect2.vcf.gz --contamination-table ${mutectdir}/${name}.contamination.table -O ${mutectdir}/${name}.mutect2.filtered.vcf.gz --tumor-segmentation ${mutectdir}/${name}.tumseg.txt --orientation-bias-artifact-priors ${mutectdir}/${name}.read.orientation.model.tar.gz 2> ${stderrdir}/${name}_filtMutectCall_stderr.log
	
#	rm ${bamdir}/${name}.bam ${bamdir}/${name}.markdup.bam
#done

# coverage calculation

#for file in ${bamdir}/*.markdup.bqsr.bam
#do
#    name=$(basename $file .markdup.bqsr.bam)
#    ${mosdepth} --by ${intervals} ${coverage}/${name}.coverage -t 4 ${bamdir}/${name}.markdup.bqsr.bam

#done


# filtering vcf file 

#for file in ${mutectdir}/*.mutect2.filtered.vcf.gz; do 

#    echo ${file}
    
#    name=$(basename ${file} .mutect2.filtered.vcf.gz)
    
#    echo ${name}
    
#    bcftools view -O v -i'FILTER == "PASS" || FILTER == "germline" & POPAF > 2' ${mutectdir}/${name}.mutect2.filtered.vcf.gz | \
#    bcftools view -i 'FORMAT/DP[0] >= 50' | grep -v -f ${chrToExclude} | grep -v -E '^chrX|^chrY' > ${mutectdir}/${name}.mutect2.pass.vcf  

#    bgzip ${mutectdir}/${name}.mutect2.pass.vcf
#    tabix -p vcf ${mutectdir}/${name}.mutect2.pass.vcf.gz


#done


#bcftools merge -m none $(ls ${mutectdir}/*.mutect2.pass.vcf.gz) | bcftools view -i 'count(FORMAT/GT[*]!="./.")>=2' > ${mutectdir}/multisample.mutect2.vcf


#vep annotation

#perl /home/alessio/programs/ensembl-vep/vep -i ${mutectdir}/multisample.mutect2.vcf -o ${mutectdir}/multisample.mutect2.vep.vcf.gz --everything --species homo_sapiens --vcf --compress_output gzip --stats_file ${mutectdir}/multisample.vep.summary.html --assembly GRCh38 --no_progress --buffer_size 5000 --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --cache --cache_version 109 --merged --tsl --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --dir_cache /home/alessio/.vep --fasta /home/simone/mnt/part1/resources/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa --format vcf --offline --pubmed --fork 10 --sift b --polyphen b --af --af_1kg --af_gnomade --af_gnomadg --regulatory

#perl /home/alessio/programs/ensembl-vep/filter_vep -i ${mutectdir}/multisample.mutect2.vep.vcf.gz --filter "EUR_AF < 0.01 or not EUR_AF" | \

#perl /home/alessio/programs/ensembl-vep/filter_vep --filter "gnomADe_NFE_AF < 0.01 or not gnomADe_NFE_AF" | \

#perl /home/alessio/programs/ensembl-vep/filter_vep --filter "gnomADg_NFE_AF < 0.01 or not gnomADg_NFE_AF" | \

#perl /home/alessio/programs/ensembl-vep/filter_vep -o ${mutectdir}/multisample.mutect2.vep.filt.vcf --filter "MAX_AF < 0.01 or not MAX_AF"


bcftools +split-vep -o ${mutectdir}/multisample.mutect2.vep.filt.tmp.bed ${mutectdir}/multisample.mutect2.vep.filt.vcf -f "%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%IMPACT\t%SYMBOL\t%Gene\t%Feature_type\t%Feature\t%BIOTYPE\t%EXON\t%INTRON\t%HGVSc\t%HGVSp\t%CDS_position\t%Amino_acids\t%Codons\t%Existing_variation\t%DISTANCE\t%PICK\t%VARIANT_CLASS\t%CANONICAL\t%MANE_SELECT\t%MANE_PLUS_CLINICAL\t%TSL\t%ENSP\t%UNIPROT_ISOFORM\t%REFSEQ_MATCH\t%GENE_PHENO\t%SIFT\t%PolyPhen\t%DOMAINS\t%AF\t%EUR_AF\t%gnomADe_AF\t%gnomADe_NFE_AF\t%gnomADg_AF\t%gnomADg_NFE_AF\t%MAX_AF\t%CLIN_SIG\t%SOMATIC\t%PHENO\t%PUBMED\t%MOTIF_NAME\t%HIGH_INF_POS\t%MOTIF_SCORE_CHANGE\t%TRANSCRIPTION_FACTORS\t[%GT;%AD;%AF\t]\n" -d -A tab 

#bash ${countSamples} ${mutectdir}/multisample.mutect2.vep.filt.vcf ${mutectdir}/multisample.mutect2.vep.filt.tmp.bed ${mutectdir}
bash ${countSamples} ${mutectdir}/multisample.mutect2.vep.filt.tmp.bed ${mutectdir}


echo -e "CHROM\tPOS\tREF\tALT\tsamplesID\tSUPP\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tCDS_position\tAmino_acids\tCodons\tExisting_variation\tDISTANCE\tPICK\tVARIANT_CLASS\tCANONICAL\tMANE_SELECT\tMANE_PLUS_CLINICAL\tTSL\tENSP\tUNIPROT_ISOFORM\tREFSEQ_MATCH\tGENE_PHENO\tSIFT\tPolyPhen\tDOMAINS\tAF\tEUR_AF\tgnomADe_AF\tgnomADe_NFE_AF\tgnomADg_AF\tgnomADg_NFE_AF\tMAX_AF\tCLIN_SIG\tSOMATIC\tPHENO\tPUBMED\tMOTIF_NAME\tHIGH_INF_POS\tMOTIF_SCORE_CHANGE\tTRANSCRIPTION_FACTORS\tA0030\tA0340\tA0734\tA0789\tA0934\tA1018\tA1428\tC1294\tC3068\tC3142" > ${mutectdir}/header.txt

paste <(cut -f1-4 ${mutectdir}/multisample.mutect2.vep.filt.tmp.bed) ${mutectdir}/countSamples.tsv <(cut -f5- ${mutectdir}/multisample.mutect2.vep.filt.tmp.bed) > ${mutectdir}/multisample.mutect2.vep.filt.bed







