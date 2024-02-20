rule get_vep_cache:
	output:
		directory("resources/vep/cache"),
	params:
		species="homo_sapiens_merged",
		build="GRCh38",
		release="110",
	log:
		"logs/vep_cache.log",
	cache: "omit-software"  # save space and time with between workflow caching (see docs)
	wrapper:
		"v3.3.5/bio/vep/cache"


rule download_vep_plugins:
	output:
		directory("resources/vep/plugins")
	params:
		release=110
	wrapper:
		"v3.3.5/bio/vep/plugins"