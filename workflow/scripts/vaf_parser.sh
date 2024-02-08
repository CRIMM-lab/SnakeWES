#!/bin/bash

# Extract AF from multisample.vcf (for each variants for all samples) 

INPUT_VCF=$1
OUTPUT_TXT=$2
#REMOTE_MACHINE_DIR="/home/alessio/testdir/wesMasto/."



# Check if the file exists
if [ -e "$INPUT_VCF" ]; then

	# Check if the file is in gzip format
	
	file_type=$(file -b "$INPUT_VCF")

	if [[ "$file_type" == *gzip* ]]; then

		echo "File is in gzip format."
	
	else
				
		gzip "$INPUT_VCF"
		
		
		if [ $? -eq 0 ]; then
		
			echo "File successfully gzipped."
		
			
			if [ "${INPUT_VCF: -3}" == ".gz" ]; then
			
		    	INPUT_VCF=${INPUT_VCF}.gz 
    			
    			echo "Renamed gzipped file"
			
    		else

    			echo "fail to rename gzipped file"
				
				exit 1 
			
			fi


		else
		
			echo "Error gziping the file."
			

		fi
	
	fi

else

	echo "File does not exist."

	exit 1
fi


echo "$INPUT_VCF"

if zgrep -qi "freebayes" "$INPUT_VCF"; then
	
	echo "The file contains 'freebayes' (case-insensitive). Continue with your script."
	
	## parsing ALLELE frequency in freebayesLowfreq

	bcftools query -f'[%AO\t%RO\n]' "$INPUT_VCF" | awk 'OFS=FS="\t"''{if($1 != ".")print}' | \

	while read -r AO RO; do
			
		if [[ $AO =~ ^[0-9]+$ && $RO =~ ^[0-9]+$ ]]; then

			result=$(echo "scale=2; $AO / ($AO + $RO)" | bc)
		
			echo -e "0$result" >> "$OUTPUT_TXT"
		
		else

			echo -e "Invalid input: $AO\t$RO"

			exit 1
		
		fi
	
	done 	

elif zgrep -qi "mutect2" "$INPUT_VCF"; then  
	
	echo "The file contains 'mutect2' (case-insensitive). Continue with your script."
	
	bcftools query -f'[%AF\t]\n' "$INPUT_VCF" | awk '{ for (i=1; i<=NF; i++) if ($i != ".") print $i }' > "$OUTPUT_TXT"

	
elif zgrep -qi "biodiversa" "$INPUT_VCF"; then  

		echo "The file contains 'biodiversa' (case-insensitive). Continue with your script."

		bcftools query -f'[%AF\t]\n' "$INPUT_VCF" | awk '{ for (i=1; i<=NF; i++) if ($i != ".") print $i }' > "$OUTPUT_TXT"

	
else

	echo "Error: check the name of the file"

	exit 1

fi  



## Copy file to remote host to remote
#scp "$OUTPUT_TXT" "alessio@192.168.1.40:${REMOTE_MACHINE_DIR}"


# Check the exit status of the scp command
#if [ $? -eq 0 ]; then

#	echo "File copied successfully."

#else
	
#	echo "Error: File copy failed."
	
#	exit 1

#fi

#rm "$OUTPUT_TXT"

echo "Done"
