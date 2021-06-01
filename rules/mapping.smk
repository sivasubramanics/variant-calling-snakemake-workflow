rule minimap_mapping:
"""
	Rule to map the trimmed fastq reads on the indexed genome and the ouput is piped to samtools, in order sort the bam entries
"""
	input:
		genome_index = WORKINGDIR + "reference/genome.mmi",
		genome_fasta = WORKINGDIR + "reference/genome.fasta",
		read1 = WORKINGDIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
		read2 = WORKINGDIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz"
	output:
		temp(OUTPUTDIR + "bamfiles/" + "{sample}.sort.bam")
	params:
		readgroup = "'@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA'"
	threads:
		10
	message:
		"Miminap2 mapping of {wildcards.sample} on the reference genome..!!"
	log:
		"logs/minimap2/mapping/{sample}.log"
	run:
		if is_single_end(wildcards.sample):
			shell(
			"{minimap2} -ax sr \
			-R {params.readgroup} \
			-t {threads} \
			-2 \
			--secondary=no \
			{input.genome_index} \
			{input.read1} 2> {log}| {samtools} sort \
			-@ {threads} \
			--reference {input.genome_fasta} \
			-o {output} - 2> {log}")
		else:
			shell(
			"{minimap2} -ax sr \
			-R {params.readgroup} \
			-t {threads} \
			-2 \
			--secondary=no \
			{input.genome_index} \
			{input.read1} \
			{input.read2} 2> {log} | {samtools} sort \
			-@ {threads} \
			--reference {input.genome_fasta} \
			-o {output} - 2> {log}")

rule bam_index:
"""
	Rule to index the sorted bam files
"""
	input:
		OUTPUTDIR + "bamfiles/{sample}.sort.bam"
	output:
		temp(OUTPUTDIR + "bamfiles/{sample}.sort.bam.bai")
	threads:
		8
	message:
		"Building bam index for the sample {wildcards.sample}..!!"
	log:
		"logs/samtools/bamindex/{sample}.log"
	wrapper:
		"0.74.0/bio/samtools/index"

rule mark_duplicates:
"""
	Rule to mark PCR duplicated reads on the bam file, using picard
"""
    input:
        OUTPUTDIR + "bamfiles/{sample}.sort.bam"
    output:
        bam = OUTPUTDIR + "bamfiles/{sample}.dedup.bam",
        metrics = OUTPUTDIR + "bamfiles/{sample}.dedup.metrics.txt"
    params:
        "REMOVE_DUPLICATES=true"
    resources:
        mem_mb=100000
    message:
    	"Marking duplicates for the sample {wildcards.sample}..!!"
    log:
        "logs/picard/dedup/{sample}.log"
    wrapper:
        "0.74.0/bio/picard/markduplicates"

rule dedup_bam_index:
"""
	Rule to index deduped bam files.
"""
	input:
		OUTPUTDIR + "bamfiles/{sample}.dedup.bam"
	output:
		OUTPUTDIR + "bamfiles/{sample}.dedup.bam.bai"
	threads:
		8
	message:
		"Building bam index for the sample {wildcards.sample}..!!"
	log:
		"logs/samtools/bamindex/{sample}.log"
	wrapper:
		"0.74.0/bio/samtools/index"
