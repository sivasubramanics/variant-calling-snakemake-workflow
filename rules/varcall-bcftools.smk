rule list_bams:
"""
	Rule to create a list of bam files.
"""
	input:
		expand(OUTPUTDIR + "bamfiles/" + "{sample}.dedup.bam", sample=SAMPLES),
		expand(OUTPUTDIR + "bamfiles/" + "{sample}.dedup.bam.bai", sample=SAMPLES)
	output:
		temp(OUTPUTDIR + "SortedBams")
	params:
		bam_dir = OUTPUTDIR + "bamfiles/*.dedup.bam"
	message:
		"Creating bam list file..!!"
	shell:
		r"ls -c1 {params.bam_dir} > {output}"

rule bcftools_varcalling:
"""
	Rule to call variations using bcftools.
"""
	input:
		genome_fasta = WORKINGDIR + "reference/genome.fasta",
		sortbamlist = OUTPUTDIR + "SortedBams"
	output:
		BCFTOOLS_VCFFILE
	params:
		"" # optional params string
	message:
		"Calling variants using bcftools..!!"
	log:
		"logs/bcftools/mpileup/stdout.log"
	threads:
		workflow.cores * 0.5
	shell:
		"{bcftools} mpileup -a DP,AD -O u --threads {threads} -f {input.genome_fasta} --bam-list {input.sortbamlist} 2> {log} | {bcftools} call --threads {threads} -mv -O z -o {output} 2>> {log}"

