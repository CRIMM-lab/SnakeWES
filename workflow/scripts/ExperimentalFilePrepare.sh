#!/bin/bash

### EXCAVATOR2 

# Prepare ExperimentalFilePrepare.txt for RC and downstream analysis 
bamdir="/home/simone/mnt/part1/masto_wes/data/bam"
outdir="/home/alessio"

# path to tumor samples bam files
find ${bamdir} -type f -iname '[A|C][0-9]*.dedup.recal.bam' > tum_bam_list.txt

# path to control samples bam files 
find ${bamdir} -type f -iname '[C][D][3]*.dedup.recal.bam' > ctr_bam_list.txt

# Unify tumor and control files
[ -s tum_bam_list.txt ] && [ -s ctr_bam_list.txt ] && cat tum_bam_list.txt ctr_bam_list.txt > all_bam_list.txt && rm tum_bam_list.txt ctr_bam_list.txt

# tumor samples names
find ${bamdir} -type f -iname '[A|C][0-9]*.dedup.recal.bam' -exec basename {} \; | sed -e's/.dedup.recal.bam//g' > tum_names_list.txt

# control samples names 
find ${bamdir} -type f -iname '[C][D][3]*.dedup.recal.bam' -exec basename {} \; | sed -e's/.dedup.recal.bam//g' > ctr_names_list.txt  

[ -s tum_names_list.txt ] && [ -s ctr_names_list.txt ] && cat tum_names_list.txt ctr_names_list.txt > all_names_list.txt 

# tumor samples path 
sed -e "s|^|${outdir}/|" tum_names_list.txt > tum_outdir_list.txt && rm tum_names_list.txt

# control samples path 
sed -e "s|^|${outdir}/|" ctr_names_list.txt > ctr_outdir_list.txt && rm ctr_names_list.txt

[ -s tum_outdir_list.txt ] && [ -s ctr_outdir_list.txt ] && cat tum_outdir_list.txt ctr_outdir_list.txt > all_outdir_list.txt && rm tum_outdir_list.txt ctr_outdir_list.txt

paste -d' ' all_bam_list.txt all_outdir_list.txt all_names_list.txt > ExperimentalFilePrepare.w100000.txt && rm all_bam_list.txt all_outdir_list.txt all_names_list.txt
