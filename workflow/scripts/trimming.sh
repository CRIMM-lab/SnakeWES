########## WHOLE GENOME SEQUENCING ANALYSIS WITH GATK #######################

# Define project

project="/home/simone/mnt/part1/masto_wes"

# Set Variables

threads=20

fastqdir="/home/simone/mnt/qnap/Danilo/Biodiversa/CD3+" && mkdir -p ${fastqdir}/trimming_fastp
logs="${project}/logs" && mkdir -p ${logs}

# Programs
fastp="/home/simone/programs/fastp/fastp"


for fastq in $(find ${fastqdir} -type f -iname '*_SA_L001_R1_001.fastq.gz'); do 

	echo $fastq

	name=$(basename "$fastq" _SA_L001_R1_001.fastq.gz)

	echo $name

	R1="${fastqdir}/${name}_SA_L001_R1_001.fastq.gz"

	R2="${fastqdir}/${name}_SA_L001_R2_001.fastq.gz"

	echo -e "$R1\t$R2"

	${fastp} -i ${R1} -I ${R2} -o ${fastqdir}/trimming_fastp/${name}.R1.fastq.gz -O ${fastqdir}/trimming_fastp/${name}.R2.fastq.gz 	\
	-j ${fastqdir}/trimming_fastp/${name}.json -h ${fastqdir}/trimming_fastp/${name}.html -w ${threads} \
	--failed_out ${fastqdir}/trimming_fastp/${name}.failedreads.txt --length_required 30 --detect_adapter_for_pe --disable_quality_filtering 2> ${logs}/${name}_fastp.log

done
