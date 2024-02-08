#	vep annotation 

mutectdir="/home/simone/mnt/part1/wesMastobis/results/mutect2"

perl /home/alessio/programs/ensembl-vep/vep -i ${mutectdir}/multisample.mutect2.vcf -o ${mutectdir}/multisample.mutect2.vep.vcf.gz --everything --species homo_sapiens --vcf --compress_output gzip --stats_file multisample.vep.summary.html --assembly GRCh38 --no_progress --buffer_size 5000 --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --cache --cache_version 109 --merged --tsl --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --dir_cache /home/alessio/.vep --fasta /home/simone/mnt/part1/resources/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa --format vcf --offline --pubmed --fork 10 --sift b --polyphen b --af --af_1kg --af_gnomade --af_gnomadg --regulatory

perl /home/alessio/programs/ensembl-vep/filter_vep -i ${mutectdir}/multisample.mutect2.vep.vcf.gz --filter "EUR_AF < 0.01 or not EUR_AF" | \

perl /home/alessio/programs/ensembl-vep/filter_vep --filter "gnomADe_NFE_AF < 0.01 or not gnomADe_NFE_AF" | \

perl /home/alessio/programs/ensembl-vep/filter_vep --filter "gnomADg_NFE_AF < 0.01 or not gnomADg_NFE_AF" | \

perl /home/alessio/programs/ensembl-vep/filter_vep -o ${mutectdir}/multisample.mutect2.vep.filt.vcf --filter "MAX_AF < 0.01 or not MAX_AF"


bcftools +split-vep -o ${mutectdir}/multisample.mutect2.vep.filt.bed ${mutectdir}/multisample.mutect2.vep.filt.vcf -f "%CHROM\t%POS\tALT\t%REF\t%Consequence\t%IMPACT\t%SYMBOL\t%Gene\t%Feature_type\t%Feature\t%BIOTYPE\t%EXON\t%INTRON\t%HGVSc\t%HGVSp\t%CDS_position\t%Amino_acids\t%Codons\t%Existing_variation\t%DISTANCE\t%PICK\t%VARIANT_CLASS\t%CANONICAL\t%MANE_SELECT\t%MANE_PLUS_CLINICAL\t%TSL\t%ENSP\t%UNIPROT_ISOFORM\t%REFSEQ_MATCH\t%GENE_PHENO\t%SIFT\t%PolyPhen\t%DOMAINS\t%AF\t%EUR_AF\t%gnomADe_AF\t%gnomADe_NFE_AF\t%gnomADg_AF\t%gnomADg_NFE_AF\t%MAX_AF\t%CLIN_SIG\t%SOMATIC\t%PHENO\t%PUBMED\t%MOTIF_NAME\t%HIGH_INF_POS\t%MOTIF_SCORE_CHANGE\t%TRANSCRIPTION_FACTORS\t[%GT\t%AD\t%AF]\n" -d -A tab 




