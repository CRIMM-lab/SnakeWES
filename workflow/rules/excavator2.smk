rule excavator2_exp_file_prep_w100000:
	input:
		bam_ctr=expand(f"data/bam/{{control}}.dedup.recal.bam", control=config['controls'].values()),
		bam_sample=expand(f"data/bam/{{sample}}.dedup.recal.bam", sample=config['samples'].values()),
		bai_ctr=expand(f"data/bam/{{control}}.dedup.recal.bam", control=config['controls'].values()),
		bai_sample=expand(f"data/bam/{{sample}}.dedup.recal.bam", sample=config['samples'].values())
	output:
		"results/excavator2/ExperimentalFilePrepare.w100000.txt"
	params:
		bamlist="results/excavator2/bam_list.txt"
		namelist="results/excavator2/names_list.txt"
	shell:
		"cat <(find {input.bam_dir} -type f -iname '[A|C][0-9]*.dedup.recal.bam') <(find ${bamdir} -type f -iname '[C][D][3]*.dedup.recal.bam') > {params.bamlist} && "
		"cat <(find {input.bam_dir} -type f -iname '[A|C][0-9]*.dedup.recal.bam' -exec basename {} \\; | sed -e's/.dedup.recal.bam//g' <(find ${bamdir} -type f -iname '[C][D][3]*.dedup.recal.bam' -exec basename {} \\; | sed -e's/.dedup.recal.bam//g' > {params.namelist} "
