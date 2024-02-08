# perform vcf to maf convertion 

VCF2MAF="/home/alessio/programs/vcf2maf-1.6.21/vcf2maf.pl"
FREEBAYES_DIR="/home/simone/mnt/part1/wesMastobis/results/freebayes/lowfreq"
REF_GEN="/home/simone/mnt/part1/resources/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
threads=20


# Define a function to process a single sample
process_sample() {
    local s="$1"    

    # Pipe the VCF data to the temporary file
    bcftools view -s ${s} multisample.freebayesLowFreq.vep.filt.vcf | bcftools view -e 'FORMAT/GT == "./."' > ${s}.tmp.vcf

    # Run the Perl script with the temporary file as input
    perl /home/alessio/programs/vcf2maf-1.6.21/vcf2maf.pl --input-vcf ${s}.tmp.vcf --output-maf ${s}.maf --tumor-id ${s} --ncbi-build GRCh38 --inhibit-vep --ref-fasta /home/simone/mnt/part1/resources/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa

    # Remove the temporary file
    rm ${s}.tmp.vcf
}

# Export the function so it's available to parallel processes
export -f process_sample

# List of samples to process (replace with your list)
sample_list=( $(bcftools query -l ${FREEBAYES_DIR}/multisample.freebayesLowFreq.vep.vcf.gz) )

# Use 'parallel' to process samples in parallel
echo "${sample_list[@]}" | tr ' ' '\n' | parallel -j ${threads} process_sample {}


# Unify all singlesamples maf files
head -n2 A0340.maf > header.txt 
parallel -j 6 'awk "FNR > 2" {}' ::: *.maf > multisample.tmp.maf 
cat header.txt multisample.tmp.maf > multisample.maf 
rm multisample.tmp.maf


#find ${FREEBAYES_DIR} -type f -iname '*.freebayes.filt.vep.vcf' | xargs -I {} -n 1 -P ${threads} \

#bash -c 'perl /home/alessio/programs/vcf2maf-1.6.21/vcf2maf.pl --input-vcf {} --output-maf ./$(basename {} .vcf).maf --tumor-id $(basename {} .freebayes.filt.vep.vcf) --ncbi-build GRCh38 --inhibit-vep --ref-fasta /home/simone/mnt/part1/resources/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa'

#head -n2 ${FREEBAYES_DIR}/A0030.freebayes.filt.vep.maf > ${FREEBAYES_DIR}/header.txt && awk 'FNR > 2' ${FREEBAYES_DIR}/*.freebayes.filt.vep.maf > ${FREEBAYES_DIR}/multisample.tmp.maf \

#&& cat ${FREEBAYES_DIR}/header.txt ${FREEBAYES_DIR}/multisample.tmp.maf > ${FREEBAYES_DIR}/multisample.maf && rm ${FREEBAYES_DIR}/multisample.tmp.maf

#bgzip ${FREEBAYES_DIR}/multisample.maf


