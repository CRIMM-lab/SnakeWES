
# Define project

project="/home/simone/mnt/part1/masto_wes"

# Set Variables
threads=20

# Results
mutectdir="${project}/results/mutect2" && mkdir -p ${mutectdir}
freebayesdir="${project}/results/freebayes" && mkdir -p ${freebayesdir}
biodiversadir="${project}/results/biodiversa" && mkdir -p ${biodiversadir}

# Script
makebed="${project}/scripts/makebed_wes.py"

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
#freebayes="/home/simone/programs/freebayes/build/freebayes"

# Filtering vcf biodiversa

#for file in ${biodiversadir}/*_SA_L001_R1_001.vcf; do

#	echo "$file"
	 
#	name=$(basename ${file} _SA_L001_R1_001.vcf)
#	echo "${name}" > ${name}.rh.txt
#
#	bcftools reheader -s ${name}.rh.txt ${biodiversadir}/${name}_SA_L001_R1_001.vcf | bcftools view -O v -i'FILTER == "PASS" || FILTER == "germline"' | bcftools norm -m - | \
#	bcftools view -i 'FORMAT/DP[0] >= 50' | grep -v -f ${chrToExclude} | grep -v -E '^chrX|^chrY' > ${biodiversadir}/${name}.biodiversa.pass.vcf && rm ${name}.rh.txt 
#
#	bgzip -f ${biodiversadir}/${name}.biodiversa.pass.vcf
#	tabix -p vcf ${biodiversadir}/${name}.biodiversa.pass.vcf.gz

#done

#bcftools merge -m none $(ls ${biodiversadir}/*.biodiversa.pass.vcf.gz) > ${biodiversadir}/multisample.biodiversa.vcf #bcftools view -i 'count(FORMAT/GT[*]!="./.")>=2'


#vep annotation

perl /home/alessio/programs/ensembl-vep/vep -i ${biodiversadir}/multisample.biodiversa.vcf.gz -o ${biodiversadir}/multisample.biodiversa.vep.vcf.gz --everything --species homo_sapiens --vcf --compress_output gzip \
	--stats_file ${biodiversadir}/multisample.vep.summary.html --assembly GRCh38 --no_progress --buffer_size 5000 --ccds --uniprot --hgvs --symbol --numbers --domains \
	--gene_phenotype --canonical --protein --biotype --uniprot --cache --cache_version 109 --merged --tsl --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number \
	--no_escape --xref_refseq --failed 1 --pick --pick_order canonical,tsl,biotype,rank,ccds,length --dir_cache /home/alessio/.vep \
	--fasta /home/simone/mnt/part1/resources/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa --format vcf --offline --pubmed --fork ${threads} \
	--sift b --polyphen b --af --af_1kg --af_gnomade --af_gnomadg --regulatory

perl /home/alessio/programs/ensembl-vep/filter_vep -i ${biodiversadir}/multisample.biodiversa.vep.vcf.gz --filter "EUR_AF < 0.01 or not EUR_AF" | \
perl /home/alessio/programs/ensembl-vep/filter_vep --filter "gnomADe_NFE_AF < 0.01 or not gnomADe_NFE_AF" | \
perl /home/alessio/programs/ensembl-vep/filter_vep --filter "gnomADg_NFE_AF < 0.01 or not gnomADg_NFE_AF" | \
perl /home/alessio/programs/ensembl-vep/filter_vep -o ${biodiversadir}/multisample.biodiversa.vep.filt.vcf --filter "MAX_AF < 0.01 or not MAX_AF"


bcftools +split-vep -o ${biodiversadir}/multisample.biodiversa.vep.filt.tmp.bed ${biodiversadir}/multisample.biodiversa.vep.filt.vcf -f "%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%IMPACT\t%SYMBOL\t%Gene\t%Feature_type\t%Feature\t%BIOTYPE\t%EXON\t%INTRON\t%HGVSc\t%HGVSp\t%CDS_position\t%Amino_acids\t%Codons\t%Existing_variation\t%DISTANCE\t%VARIANT_CLASS\t%CANONICAL\t%MANE_SELECT\t%MANE_PLUS_CLINICAL\t%TSL\t%ENSP\t%UNIPROT_ISOFORM\t%REFSEQ_MATCH\t%GENE_PHENO\t%SIFT\t%PolyPhen\t%DOMAINS\t%AF\t%EUR_AF\t%gnomADe_AF\t%gnomADe_NFE_AF\t%gnomADg_AF\t%gnomADg_NFE_AF\t%MAX_AF\t%CLIN_SIG\t%SOMATIC\t%PHENO\t%PUBMED\t%MOTIF_NAME\t%HIGH_INF_POS\t%MOTIF_SCORE_CHANGE\t%TRANSCRIPTION_FACTORS[\t%AD][\t%AF][\t%GT]\n" -d -A tab 

bcftools query -l ${biodiversadir}/multisample.biodiversa.vep.filt.vcf > samples_list.txt 

python ${makebed} multisample.biodiversa.vep.filt.tmp.bed multisample.biodiversa.vep.filt.bed samples_list.txt 

#echo -e "CHROM\tPOS\tREF\tALT\tsamplesID\tSUPP\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tCDS_position\tAmino_acids\tCodons\tExisting_variation\tDISTANCE\tPICK\tVARIANT_CLASS\tCANONICAL\tMANE_SELECT\tMANE_PLUS_CLINICAL\tTSL\tENSP\tUNIPROT_ISOFORM\tREFSEQ_MATCH\tGENE_PHENO\tSIFT\tPolyPhen\tDOMAINS\tAF\tEUR_AF\tgnomADe_AF\tgnomADe_NFE_AF\tgnomADg_AF\tgnomADg_NFE_AF\tMAX_AF\tCLIN_SIG\tSOMATIC\tPHENO\tPUBMED\tMOTIF_NAME\tHIGH_INF_POS\tMOTIF_SCORE_CHANGE\tTRANSCRIPTION_FACTORS\tA0030\tA0340\tA0734\tA0789\tA0934\tA1018\tA1428\tC1294\tC3068\tC3142" > header.txt

# plot AF for each variants for all samples 
mkdir ${biodiversadir}/plots

bash ${project}/scripts/vaf_parser.sh multisample.biodiversa.vep.filt.vcf.gz plot/af_biodiversa_filt.txt

Rscript ${biodiversadir}/plots/densityplot.R ${biodiversadir}/plots/af_biodiversa_filt.txt ${biodiversadir}/plots/af_biodiversa_filt.pdf
